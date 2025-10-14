# import duckdb as ddb
# import polars as pl
import zipfile
import geopandas as gpd
from shapely.ops import substring
from shapely.geometry import LineString, Point
import numpy as np
import pandas as pd

# import similaritymeasures

import functools
# from loguru import logger 

# loguru is an external logging library that QGIS’s Python environment does not include it by default.
# so we did the following:
import logging
logger = logging.getLogger("vehicle_passenger_flow") 
logging.basicConfig(level=logging.INFO) # This ensures that .info(), .debug(), .error() will actually output something.

# THIS IS A NEW HELPER FUNCTION TO TURN MULTILINESTRING TO LINESTRING IN FRECHET DISTANCE CALCULATION
from shapely.geometry import MultiLineString, GeometryCollection
from shapely.ops import linemerge
def to_single_linestring(geom):
    """
    Return a LineString for frechet distance:
    - LineString => itself
    - MultiLineString => longest component (or linemerge if it helps)
    - GeometryCollection => try merging/extracting the longest LineString
    - Anything else / empty => None
    """
    if geom is None or geom.is_empty:
        return None
    if isinstance(geom, LineString):
        return geom
    if isinstance(geom, MultiLineString):
        # pick the longest linestring component if any
        if len(geom.geoms) == 0:
            return None
        longest = max(geom.geoms, key=lambda g: g.length)
        # try to merge first if it helps, else use longest
        merged = linemerge(geom)
        if isinstance(merged, LineString):
            return merged
        return longest
    if isinstance(geom, GeometryCollection):
        # extract all lines, pick the longest
        lines = [g for g in geom.geoms if isinstance(g, (LineString, MultiLineString))]
        if not lines:
            return None
        merged = linemerge(lines)
        if isinstance(merged, LineString):
            return merged
        if isinstance(merged, MultiLineString) and len(merged.geoms) > 0:
            return max(merged.geoms, key=lambda g: g.length)
        # fallback to longest single linestring among originals
        flat = []
        for g in lines:
            if isinstance(g, LineString):
                flat.append(g)
            elif isinstance(g, MultiLineString):
                flat.extend(list(g.geoms))
        return max(flat, key=lambda g: g.length) if flat else None
    return None


import math
def discrete_frechet(P, Q):
    if len(P) == 0 or len(Q) == 0:  # WE ADDED THE IF CONDITION AS A DEFENSIVE CHECK TO AVOID AN ERROR
        return float("inf")  # or return None
    """
    Discrete Frechet distance between curves P and Q (each a list of [x, y] coordinates).
    """
    ca = [[-1 for _ in range(len(Q))] for _ in range(len(P))]

    def c(i, j):
        if ca[i][j] > -1:
            return ca[i][j]
        elif i == 0 and j == 0:
            ca[i][j] = math.dist(P[0], Q[0])
        elif i > 0 and j == 0:
            ca[i][j] = max(c(i - 1, 0), math.dist(P[i], Q[0]))
        elif i == 0 and j > 0:
            ca[i][j] = max(c(0, j - 1), math.dist(P[0], Q[j]))
        elif i > 0 and j > 0:
            ca[i][j] = max(
                min(c(i - 1, j), c(i - 1, j - 1), c(i, j - 1)),
                math.dist(P[i], Q[j])
            )
        else:
            ca[i][j] = float("inf")
        return ca[i][j]

    return c(len(P) - 1, len(Q) - 1)


# @functools.cache
def read_trips_gdf_from_gtfs(gtfs_path):
    logger.info("Extracting trips GDF from GTFS feed")
    gtfs_archive = zipfile.ZipFile(gtfs_path, "r")

    routes_df = pd.read_csv(gtfs_archive.open("routes.txt"))
    trips_df = pd.read_csv(gtfs_archive.open("trips.txt"))
    shapes_df = pd.read_csv(gtfs_archive.open("shapes.txt"))

    # Merge and sort shapes
    shapes_df = shapes_df.sort_values(by=["shape_id", "shape_pt_sequence"])
    shapes_gdf = shapes_df.groupby("shape_id").apply(lambda group: LineString(zip(group.shape_pt_lon, group.shape_pt_lat)))
    shapes_gdf = shapes_gdf.reset_index()
    shapes_gdf.columns = ["shape_id", "geometry"]
    shapes_gdf = gpd.GeoDataFrame(shapes_gdf, geometry=gpd.GeoSeries.from_wkt(shapes_gdf.geometry.astype(str)), crs="EPSG:4326")
    shapes_gdf = shapes_gdf.to_crs(3857)

    # Merge trips with routes and shapes
    trips_df = trips_df.merge(routes_df, on="route_id", how="left")
    trips_df = trips_df.merge(shapes_gdf, on="shape_id", how="left")

    trips = gpd.GeoDataFrame(trips_df, geometry="geometry", crs="EPSG:3857")
    
    return trips


# @functools.cache
def get_trips_segments_from_gtfs(gtfs_path):
    # def gtfs2gps(gtfs_feed_path):
    gtfs_archive = zipfile.ZipFile(gtfs_path, "r")

    stops_df = pd.read_csv(gtfs_archive.open("stops.txt"))
    stop_times_df = pd.read_csv(gtfs_archive.open("stop_times.txt"))

    # Convert to GeoDataFrame
    stops = gpd.GeoDataFrame(
        stops_df,
        geometry=gpd.points_from_xy(stops_df["stop_lon"], stops_df["stop_lat"]),
        crs="EPSG:4326"
    ).to_crs("EPSG:3857")

    stop_times_df["arrival_time"] = pd.to_timedelta(stop_times_df["arrival_time"])
    stop_times_df["departure_time"] = pd.to_timedelta(stop_times_df["departure_time"])

    stop_times_df = stop_times_df.sort_values(["trip_id", "stop_sequence"])
    stop_times_df["segment_order"] = stop_times_df.groupby("trip_id").cumcount()

    # Lag stop_id and time
    stop_times_df["from_id"] = stop_times_df.groupby("trip_id")["stop_id"].shift()
    stop_times_df["to_id"] = stop_times_df["stop_id"]
    stop_times_df["start_time"] = stop_times_df.groupby("trip_id")["arrival_time"].shift()
    stop_times_df["end_time"] = stop_times_df["arrival_time"]
    stop_times_df = stop_times_df.dropna(subset=["from_id", "start_time"])

    stop_times_df["start_time_secs"] = stop_times_df["start_time"].dt.total_seconds()
    stop_times_df["duration_secs"] = (stop_times_df["end_time"] - stop_times_df["start_time"]).dt.total_seconds()
    stop_times_df["trip_start_time"] = stop_times_df.groupby("trip_id")["start_time_secs"].transform("min")
    stop_times_df["start_time"] = stop_times_df["start_time_secs"] - stop_times_df["trip_start_time"]

    # Merge stops geometry
    trips = read_trips_gdf_from_gtfs(gtfs_path)

    # THIS SECTION IS A REPLACEMENT TO THE COMMENTED OUT PARTS ABOVE
    segments = stop_times_df.copy()
    trip_geometry_map = trips.set_index("trip_id")["geometry"].to_dict()
    segments["trip_geometry"] = segments["trip_id"].map(trip_geometry_map)
    # First, create a dictionary mapping stop_id to geometry
    stop_geometry_map = stops.set_index("stop_id")["geometry"].to_dict()
    # Assign geometry for from_id and to_id directly
    segments["geometry_from"] = segments["from_id"].map(stop_geometry_map)
    segments["geometry_to"] = segments["to_id"].map(stop_geometry_map)
    segments["geometry"] = segments.apply(
        lambda row: substring(row.trip_geometry, row.trip_geometry.project(row.geometry_from), row.trip_geometry.project(row.geometry_to))
        if isinstance(row.trip_geometry, LineString) else None,
        axis=1
    )


    trips_segments_gdf = gpd.GeoDataFrame(
        segments,
        geometry="geometry",
        crs="EPSG:3857"
    )

    # Compute speed
    trips_segments_gdf["speed_m_s"] = trips_segments_gdf.geometry.length / trips_segments_gdf.duration_secs
    trips_segments_gdf["length"] = trips_segments_gdf.geometry.length
    trips_segments_gdf = trips_segments_gdf[trips_segments_gdf["length"] > 0]

    return trips_segments_gdf, stops


def get_segments(curve):
    return list(map(LineString, zip(curve.coords[:-1], curve.coords[1:])))


# @functools.cache
def trips_segments_to_transport_model(
    trips_segments_gdf, trip_segments_segmentization_threshold_meters
):
    trips_segementized = (  # noqa: F841
    # TODO: Rename trips_segments_gdf -> gtfs_trips_segments_gdf
        trips_segments_gdf.assign(
            subsegments=lambda df: df.apply(
                lambda row: get_segments(
                    row.geometry.segmentize(
                        trip_segments_segmentization_threshold_meters
                    )
                ),
                axis=1,
            )
        )
        .explode("subsegments")
        .set_geometry("subsegments")
        .set_crs(3857)
        .drop(columns=["geometry"])
        .rename_geometry("geometry")
        .assign(
            duration_secs=lambda df: (df.geometry.length / df.speed_m_s),
            length_m=lambda df: df.geometry.length,
        )
    )


    transport_model_points = (
        trips_segementized
        .assign(
            point_order=lambda df: df.groupby("trip_id").cumcount(),
        )
        .assign(
            geometry=lambda df: df.geometry.apply(lambda geom: Point(geom.coords[0]))
        )
        .set_geometry("geometry")
        # .to_wkt()
    )
    
    transport_model = (
        transport_model_points
        .assign(
            timestamp=lambda df: df.groupby("trip_id")["duration_secs"].cumsum().shift(fill_value=0)
        )
    )

    return transport_model, trips_segementized



def get_vehicle_appearances_on_road_segments(
    trips_segments_gdf: gpd.GeoDataFrame,
    road_segments_gdf: gpd.GeoDataFrame,
    transport_model: gpd.GeoDataFrame,
    frequencies_df: pd.DataFrame,
    segment_matching_buffer_meters: float,
    frechet_dist_densify_param: float,  # kept for signature compatibility (not used by discrete version)
    frechet_dist_segment_length_ratio_max_thresh: float,
):
    """
    DuckDB-free rewrite of the original function.
    Uses GeoPandas/Shapely for spatial ops and a discrete Fréchet implementation.
    Returns the same tuple:
      (vehicle_appearances, candidates_prep, candidates_with_verdict, matched, trip_instances, temp2)
    """

    # 1) Build candidates_prep by spatial join between buffered trip segments and roads
    ts = trips_segments_gdf[["trip_id", "segment_order", "geometry"]].copy()
    rs = road_segments_gdf[["gid", "geometry"]].copy()

    # small speed-up: precompute buffers once
    ts_buf = ts.copy()
    ts_buf["buffer_geom"] = ts_buf.geometry.buffer(segment_matching_buffer_meters)

    # spatial join: roads intersect trip buffer
    # (GeoPandas >=0.10 uses predicate=; older versions use op=)
    rs_for_join = rs.rename(columns={"geometry": "rs_geom"})
    ts_for_join = ts_buf.rename(columns={"geometry": "ts_geom"})  # keep original line in ts_geom
    ts_for_join = ts_for_join.set_geometry("buffer_geom")

    joined_idx = gpd.sjoin(
        rs_for_join.set_geometry("rs_geom"),
        ts_for_join[["trip_id", "segment_order", "ts_geom", "buffer_geom"]].set_geometry("buffer_geom"),
        how="inner",
        predicate="intersects",
    ).reset_index(drop=True)

    # compute the intersections like original SQL:
    # NOTE: we need a buffer for RS too (match original logic)
    rs_buffered = rs_for_join.copy()
    rs_buffered["rs_buffer"] = rs_buffered["rs_geom"].buffer(segment_matching_buffer_meters)
    rs_buffer_map = rs_buffered.set_index("gid")["rs_buffer"].to_dict()

    def row_intersections(row):
        rs_geom = row["rs_geom"]              # road centerline
        ts_geom = row["ts_geom"]              # trip segment centerline (we carried this as a column)
        rs_buf  = rs_buffer_map.get(row["gid"])

        # buffer the trip geometry here (no need to read row["buffer_geom"])
        ts_buf = ts_geom.buffer(segment_matching_buffer_meters) if ts_geom is not None else None

        part_in_rs = rs_geom.intersection(ts_buf) if (rs_geom is not None and ts_buf is not None) else None
        part_in_ts = ts_geom.intersection(rs_buf) if (ts_geom is not None and rs_buf is not None) else None
        dist = rs_geom.distance(ts_geom) if (rs_geom is not None and ts_geom is not None) else np.nan
        return part_in_rs, part_in_ts, ts_geom, dist

    out = joined_idx.apply(lambda r: row_intersections(r), axis=1, result_type="expand")
    out.columns = ["intersected_part_in_rs", "intersected_part_in_ts", "trip_segment_geometry", "dist"]

    candidates_prep = pd.concat(
        [
            joined_idx[["trip_id", "segment_order", "gid"]].reset_index(drop=True),
            out.reset_index(drop=True),
        ],
        axis=1,
    )

    # 2) Compute Fréchet distance + lengths, robust to multipart/empty geometries
    def safe_frechet(rs_geom, ts_geom):
        rs = to_single_linestring(rs_geom)
        ts = to_single_linestring(ts_geom)
        if rs is None or ts is None:
            return np.inf
        P = list(rs.coords)
        Q = list(ts.coords)
        if not P or not Q:
            return np.inf
        return min(
            discrete_frechet(P, Q),
            discrete_frechet(list(reversed(P)), Q),
        )

    geom_rs = candidates_prep["intersected_part_in_rs"]
    geom_ts = candidates_prep["intersected_part_in_ts"]
    trip_geom = candidates_prep["trip_segment_geometry"]

    frechet_vals = [
        safe_frechet(rg, tg) for rg, tg in zip(geom_rs, geom_ts)
    ]

    intersection_rs_length = gpd.GeoSeries(geom_rs, crs=road_segments_gdf.crs).length
    intersection_ts_length = gpd.GeoSeries(geom_ts, crs=road_segments_gdf.crs).length
    trip_segment_length   = gpd.GeoSeries(trip_geom, crs=road_segments_gdf.crs).length

    candidates = gpd.GeoDataFrame(
        candidates_prep.assign(
            frechet_dist=frechet_vals,
            intersection_rs_length=intersection_rs_length,
            intersection_ts_length=intersection_ts_length,
            trip_segment_length=trip_segment_length,
        ),
        geometry="intersected_part_in_rs",
        crs=road_segments_gdf.crs,
    ).explode(index_parts=True).reset_index(drop=True)

    # 3) Verdict
    candidates_with_verdict = candidates.assign(
        accepted=(candidates["frechet_dist"] / candidates["intersection_rs_length"].replace(0, np.nan)) 
                  < frechet_dist_segment_length_ratio_max_thresh
    )
    candidates_with_verdict["accepted"] = candidates_with_verdict["accepted"].fillna(False)

    matched = candidates_with_verdict.query("accepted").copy()

    # 4) Map transport_model points to nearest matched road segment centroid for (trip_id, segment_order, gid)
    temp_df = (
        road_segments_gdf
        .assign(centroid=lambda df: df.geometry.centroid)
        # .rename_geometry("geometry")
        [["gid", "geometry", "centroid"]]
        .merge(
            transport_model.merge(
                matched[["trip_id", "segment_order", "gid"]],
                on=["trip_id", "segment_order"],
            ),
            on="gid",
        )
        .assign(
            dist=lambda df: gpd.GeoSeries(df["centroid"]).distance(
                gpd.GeoSeries(df["geometry_y"]).set_crs(road_segments_gdf.crs)
            )
        )
    )

    # Ensure both sides have CRS set, then compute distances
    temp_df = temp_df.copy()
    temp_df["dist"] = gpd.GeoSeries(temp_df["centroid"], crs=road_segments_gdf.crs).distance(
        gpd.GeoSeries(temp_df["geometry_y"]).set_crs(road_segments_gdf.crs)
    )

    temp2 = (
        temp_df.sort_values(["trip_id", "segment_order", "gid", "dist"])
               .drop_duplicates(["trip_id", "segment_order", "gid"], keep="first")
    )

    # 5) trip_instances from frequencies_df (no SQL)
    # expects columns: trip_id, start_time, end_time, headway_secs (all seconds)
    freq = frequencies_df.copy()
    # guard types
    for col in ["start_time", "end_time", "headway_secs"]:
        freq[col] = pd.to_numeric(freq[col], errors="coerce")
    freq = freq.dropna(subset=["trip_id", "start_time", "end_time", "headway_secs"])
    freq["n_instances"] = ((freq["end_time"] - freq["start_time"]) // freq["headway_secs"]).astype(int)
    freq = freq[freq["n_instances"] > 0]

    # expand instances
    freq["instance_no"] = freq["n_instances"].apply(lambda n: list(range(n)))
    trip_instances = (
        freq.explode("instance_no")
            .assign(
                instance_no=lambda df: df["instance_no"].astype(int),
                trip_start_time=lambda df: df["start_time"] + df["instance_no"] * df["headway_secs"]
            )[["trip_id", "headway_secs", "trip_start_time", "instance_no"]]
            .reset_index(drop=True)
    )

    # 6) vehicle_appearances
    vehicle_appearances = gpd.GeoDataFrame(
        trip_instances.merge(
            temp2[["trip_id", "segment_order", "duration_secs", "timestamp", "gid", "geometry_y"]]
            .rename(columns={"geometry_y": "geometry"}),
            on="trip_id",
        )
        .assign(timeofday_secs=lambda df: df["trip_start_time"] + df["timestamp"])
        .drop_duplicates(subset=["gid", "trip_id", "instance_no", "timeofday_secs"])
        [["gid", "timeofday_secs", "duration_secs", "trip_id", "segment_order", "geometry"]],
        geometry="geometry",
        crs=road_segments_gdf.crs,
    )

    return vehicle_appearances, candidates_prep, candidates_with_verdict, matched, trip_instances, temp2


# @functools.cache
def read_or_generate_frequencies(gtfs_path):
    trips_segments_gdf, stops = get_trips_segments_from_gtfs(gtfs_path)

    # Read frequencies.txt, if exists
    gtfs_archive = zipfile.ZipFile(gtfs_path, "r")
    try:
        # frequencies_txt = pl.read_csv(gtfs_archive.read("frequencies.txt")).to_pandas()
        frequencies_txt = pd.read_csv(gtfs_archive.open("frequencies.txt"))
    except BaseException:
        frequencies_txt = pd.DataFrame(  # noqa: F841
            {"trip_id": [], "start_time": [], "end_time": [], "headway_secs": []}
        )

    # Convert times to seconds since midnight
    for col in ["start_time", "end_time"]:
        frequencies_txt[col] = pd.to_timedelta(frequencies_txt[col], errors="coerce")
        frequencies_txt[col] = frequencies_txt[col].dt.total_seconds()

    frequencies_txt["headway_secs"] = pd.to_numeric(frequencies_txt["headway_secs"], errors="coerce")

    # Add synthetic frequencies where missing
    frequencies_augmented = (
        trips_segments_gdf[["trip_id", "trip_start_time"]]
        .drop_duplicates("trip_id")
        .merge(frequencies_txt, how="left", on="trip_id")
        .assign(
            start_time=lambda df: df["start_time"].fillna(df["trip_start_time"]),
            end_time=lambda df: df["end_time"].fillna(df["start_time"] + 10),
            headway_secs=lambda df: df["headway_secs"].fillna(10),
        )
    )
    return frequencies_augmented
    


def get_vehicle_appearances(
    gtfs_path,
    road_segments_gdf,
    agency_filter_list=None,
    config={
        "trip_segments_segmentization_threshold_meters": 300,
        "segment_matching_buffer_meters": 80,
        "frechet_dist_densify_param": 0.1,
        "frechet_dist_segment_length_ratio_max_thresh": 0.2,
    },
):
    # Process GTFS trips
    logger.info("Extracting trip segments from GTFS...")
    trips_segments_gdf, stops = get_trips_segments_from_gtfs(gtfs_path)

    if agency_filter_list:
        trips_segments_gdf = trips_segments_gdf[
            trips_segments_gdf["agency_id"].isin(agency_filter_list)
        ]

    logger.info("Synthesizing transport model")
    transport_model, trips_segmentized = trips_segments_to_transport_model(
        # TODO: Rewrite to passing all configuration using the reader Monad
        trips_segments_gdf, config["trip_segments_segmentization_threshold_meters"]
    )

    frequencies_df = read_or_generate_frequencies(gtfs_path)

    logger.info("Estimating number of vehicles per road segment per interval")
    vehicle_appearances, candidates_prep, candidates_with_verdict, matched, trip_instances, temp2 = (
        get_vehicle_appearances_on_road_segments(
            road_segments_gdf=road_segments_gdf,
            transport_model=transport_model,
            trips_segments_gdf=trips_segmentized,
            frequencies_df=frequencies_df,
            segment_matching_buffer_meters=config["segment_matching_buffer_meters"],
            frechet_dist_densify_param=config["frechet_dist_densify_param"],
            frechet_dist_segment_length_ratio_max_thresh=config[
                "frechet_dist_segment_length_ratio_max_thresh"
            ],
        )
    )

    # Keep legacy variable names for downstream compatibility
    candidates_base = candidates_prep
    candidate_trip_segments = candidates_with_verdict
    matched_trip_segments = matched

    # Return all these for debugging the intermediate steps
    return (vehicle_appearances,
        transport_model,
        candidates_base,
        candidate_trip_segments,
        matched_trip_segments,
        trips_segments_gdf,
        stops,
        trip_instances,
        temp2)


# @functools.cache
def get_avg_occupancy_per_segment_v2(
    trips_segments_gdf,
    stops_gdf,
    intervals_df,
    raw_stops_gdf,
    onboard_instances_gdf,
):
    trips_subsegments = (
        trips_segments_gdf.assign(
            subsegments=lambda df: df.apply(
                lambda row: get_segments(row.geometry.segmentize(300)),
                axis=1,
            )
        )
        .explode("subsegments")
        .set_geometry("subsegments")
        .set_crs(3857)
        .drop(columns=["geometry"])
        .rename_geometry("geometry")
        .assign(
            duration_secs=lambda df: (df.geometry.length / df.speed_m_s),
            length_m=lambda df: df.geometry.length,
            point_order=lambda df: df.groupby("trip_id").cumcount(),
            segment_id= lambda df: df['trip_id'] + "_" + df['point_order'].astype(str)
        )
    )
    valid_onboard_instances = (
        onboard_instances_gdf[["id", "trip_id", "status", "valid", "departed_at"]]
        # .rename(columns={"departed_a": "departed_at"})
        .query("status=='finished'")
        .query("valid==True")
    )

    valid_raw_stops = raw_stops_gdf.drop(columns="gid").merge(
        valid_onboard_instances, left_on="onboard_instance_observer_id", right_on="id"
    )

    filtered_stops = valid_raw_stops[
        valid_raw_stops["trip_id"].isin(trips_subsegments["trip_id"].unique())
    ].reset_index(drop=True)

    matched_stops = (
        filtered_stops.groupby("trip_id")
        .apply(
            lambda group: group.reset_index()
            .rename(columns={"index": "orig_index"})
            .sjoin_nearest(trips_subsegments.query(f"trip_id == '{group.name}'"))[
                ["segment_id", "orig_index"]
            ]
            .reset_index()
        )
        .rename(columns={"trip_id": "trip_id_temp"})
        .reset_index()
        .merge(
            filtered_stops.reset_index()
            .rename(columns={"index": "orig_index"})
            .drop(columns="trip_id"),
            on="orig_index",
        )
        .groupby(
            ["onboard_instance_observer_id", "departed_at", "trip_id", "segment_id"],
            as_index=False,
        )
        .agg({"board": "sum", "alight": "sum", "created_at": "first"})
    )

    onboard_segments_with_occupancy = (
        trips_subsegments.drop(columns="geometry")
        # Merge to filter trips not found in the matched_stops layer
        .merge(
            matched_stops[
                ["trip_id", "onboard_instance_observer_id"]
            ].drop_duplicates(),
            on="trip_id",
        )
        .merge(
            matched_stops.drop(columns="trip_id"),
            how="left",
            left_on=["onboard_instance_observer_id", "segment_id"],
            right_on=["onboard_instance_observer_id", "segment_id"],
        )
        .sort_values(["trip_id", "onboard_instance_observer_id", "segment_order"])
        .fillna({"board": 0, "alight": 0})
        .assign(
            survey_departed_at=lambda df: df.groupby(
                ["trip_id", "onboard_instance_observer_id"]
            )["departed_at"]
            .ffill()
            .bfill()
        )
        .assign(
            survey_departed_at=lambda df: pd.to_datetime(df.survey_departed_at)
            .dt.tz_convert("UTC")
            .dt.tz_localize(None)
        )
        .assign(
            survey_departed_at=lambda df: (
                df.survey_departed_at - df.survey_departed_at.dt.floor("d")
            ).dt.total_seconds()
        )
        .assign(
            vehicle_occupancy=lambda df: df.groupby(
                ["trip_id", "onboard_instance_observer_id"]
            )[["board", "alight"]]
            .cumsum()
            .assign(alight=lambda df: df.alight * -1)
            .sum(axis=1)
        )
    )

    avg_occupancy_per_trip_segment_per_interval = (
        onboard_segments_with_occupancy.merge(intervals_df, how="cross")
        .query(
            "survey_departed_at + start_time < interval_end_secs & survey_departed_at + start_time > interval_start_secs"
        )
        .groupby(
            [
                "trip_id",
                "segment_order",
                "interval_name",
                "interval_start_secs",
                "interval_end_secs",
            ],
            as_index=False,
        )
        .agg({"vehicle_occupancy": "median"})
    )

    return (
        avg_occupancy_per_trip_segment_per_interval,
        onboard_segments_with_occupancy,
        matched_stops,
        filtered_stops,
    )
