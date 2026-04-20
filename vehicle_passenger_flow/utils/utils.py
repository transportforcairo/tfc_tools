import geopandas as gpd
from shapely.ops import substring
from shapely.geometry import LineString, Point
from shapely.strtree import STRtree
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


def _agency_vehicle_name_map(agencies_df: pd.DataFrame, vehicles_df: pd.DataFrame) -> pd.DataFrame:
    """
    Returns a df: agency_id, vehicle_name.
    Assumes agencies has vehicle_id referencing vehicles.gid (common in your SDI).
    """
    a = agencies_df.copy()
    v = vehicles_df.copy()

    if "vehicle_id" in a.columns and "gid" in v.columns:
        out = a.merge(v[["gid", "name"]], left_on="vehicle_id", right_on="gid", how="left")
        out = out.rename(columns={"name": "vehicle_name"})
        return out[["agency_id", "vehicle_name"]].drop_duplicates()

    if "vehicle_name" in a.columns:
        return a[["agency_id", "vehicle_name"]].drop_duplicates()

    raise ValueError("Could not build agency->vehicle_name mapping (need agencies.vehicle_id + vehicles.gid or agencies.vehicle_name).")


def build_trip_stop_pair_segments_from_sdi(
    trips_view_gdf: gpd.GeoDataFrame,
    trip_stops_sequence_df: pd.DataFrame,
    agencies_df: pd.DataFrame,
    vehicles_df: pd.DataFrame,
) -> gpd.GeoDataFrame:
    """
    Builds trip stop-to-stop segments using:
      - transit.trips_view geometry (LineString)
      - transit.trip_stops_sequence distance_frac (0..1) and stop_sequence

    IMPORTANT: In RL2SDI GeoPackage exports, trip_stops_sequence often contains:
      - t_id (numeric internal trip gid)
      - observer_trip_id (string, the RouteLab/Observer trip id)
    While trips_view uses observer_id (string). So we prefer observer_trip_id when present.

    Output columns:
      trip_id, segment_order, from_id, to_id, vehicle_name, geometry, length_m
    """
    tv = trips_view_gdf.copy()
    tss = trip_stops_sequence_df.copy()

    # Attach vehicle_name to trips_view via agency->vehicle mapping (if needed)
    if "vehicle_name" not in tv.columns:
        vn = _agency_vehicle_name_map(agencies_df, vehicles_df)
        if "agency_id" in tv.columns:
            tv = tv.merge(vn, on="agency_id", how="left")

    # Expect canonical trip_id already normalized upstream
    if "trip_id" not in tss.columns:
        raise ValueError("trip_stops_sequence_df must include canonical trip_id (observer trip id).")

    tss["trip_id"] = tss["trip_id"].astype(str).str.strip()
    tv["trip_id"] = tv["trip_id"].astype(str).str.strip()

    if "stop_sequence" not in tss.columns:
        raise ValueError("trip_stops_sequence must contain stop_sequence.")
    if "stop_id" not in tss.columns:
        raise ValueError("trip_stops_sequence must contain stop_id.")
    if "distance_frac" not in tss.columns:
        raise ValueError("trip_stops_sequence must contain distance_frac (0..1).")

    # Ensure distance_frac is numeric
    tss["distance_frac"] = pd.to_numeric(tss["distance_frac"], errors="coerce")

    # --- Prep stop pairs ---
    tss = tss.sort_values(["trip_id", "stop_sequence"]).copy()
    tss["from_id"] = tss["stop_id"]
    tss["to_id"] = tss.groupby("trip_id")["stop_id"].shift(-1)
    tss["from_frac"] = tss["distance_frac"]
    tss["to_frac"] = tss.groupby("trip_id")["distance_frac"].shift(-1)
    tss["segment_order"] = tss.groupby("trip_id").cumcount()

    seg = tss.dropna(subset=["to_id", "to_frac", "from_frac"]).copy()

    # --- Map trip geometry & vehicle_name ---
    tv["trip_id"] = tv["trip_id"].astype(str).str.strip()
    geom_map = tv.set_index("trip_id")["geometry"].to_dict()
    veh_map = tv.set_index("trip_id")["vehicle_name"].to_dict() if "vehicle_name" in tv.columns else {}

    seg["trip_id"] = seg["trip_id"]
    seg["trip_geometry"] = seg["trip_id"].map(geom_map)
    seg["vehicle_name"] = seg["trip_id"].map(veh_map)

    # Defensive diagnostics: if nothing maps, fail loudly with hints
    if seg["trip_geometry"].isna().all():
        raise ValueError(
            "No trip geometries could be mapped from trips_view to trip_stops_sequence. "
            "This usually means trip IDs don't match. "
            "Check trips_view.trip_id values vs trip_stops_sequence.observer_trip_id (preferred) / trip_id / t_id."
        )

    # --- Build segment geometry using shapely substring with normalized fractions ---
    seg["geometry"] = seg.apply(
        lambda r: substring(
            r["trip_geometry"],
            float(min(r["from_frac"], r["to_frac"])),
            float(max(r["from_frac"], r["to_frac"])),
            normalized=True,
        )
        if isinstance(r["trip_geometry"], LineString)
        else None,
        axis=1,
    )

    out = gpd.GeoDataFrame(
        seg[["trip_id", "segment_order", "from_id", "to_id", "vehicle_name", "geometry"]],
        geometry="geometry",
        crs=tv.crs,
    )
    out = out[~out.geometry.isna()].copy()
    out["length_m"] = out.geometry.length
    out = out[out["length_m"] > 0].copy()

    return out


def build_flow_stop_pair_segments(trip_segments_gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Deduplicate stop-pair segments across trips (like the spirit of od_stats),
    keeping vehicle_name in the key.
    """
    cols = ["from_id", "to_id", "vehicle_name", "geometry"]
    base = trip_segments_gdf[cols].drop_duplicates().copy()
    base["flow_seg_id"] = (
        base["from_id"].astype(str) + "__" +
        base["to_id"].astype(str) + "__" +
        base["vehicle_name"].astype(str)
    )
    return gpd.GeoDataFrame(base, geometry="geometry", crs=trip_segments_gdf.crs)





def get_segments(curve):
    return list(map(LineString, zip(curve.coords[:-1], curve.coords[1:])))


def get_avg_occupancy_per_segment_v3_sdi_timeproxy(
    trip_segments_gdf,
    intervals_df,
    raw_stops_gdf,
    onboard_instances_gdf,
    od_stats_df,
    feedback=None
):
    """
    Same occupancy logic as v2, but:
      - trip segments come from SDI
      - start_time is computed per interval using od_stats.duration
      - duration join key: (o_id, d_id, interval_name, vehicle_name)
    Returns:
      avg_occupancy_per_trip_segment_per_interval, onboard_segments_with_occupancy, matched_stops, filtered_stops
        """
    
    # ------------------------------------------------------------
    # Option A: snap raw stop events directly to stop-pair segments
    # (no micro-segmentation)
    # ------------------------------------------------------------
    
    segments_targets = trip_segments_gdf[[
        "trip_id", "segment_order", "from_id", "to_id", "vehicle_name", "geometry"
    ]].copy()
    segments_targets = gpd.GeoDataFrame(segments_targets, geometry="geometry", crs=trip_segments_gdf.crs)
    segments_targets = segments_targets[~segments_targets.geometry.isna()].copy()
    
    if segments_targets.empty:
        raise ValueError("trip_segments_gdf produced no valid geometries (segments_targets is empty).")
    
    valid_onboard_instances = onboard_instances_gdf[["id", "trip_id", "status", "valid", "departed_at"]].copy()
    valid_onboard_instances = valid_onboard_instances.query("status == 'finished'")
    valid_onboard_instances["valid"] = (
    valid_onboard_instances["valid"]
    .astype(str)
    .str.strip()
    .str.lower()
    .map({"true": True, "false": False}))

    valid_onboard_instances = valid_onboard_instances[
        valid_onboard_instances["valid"] == True
    ]
    
    valid_raw_stops = raw_stops_gdf.drop(columns="gid", errors="ignore").merge(
        valid_onboard_instances,
        left_on="onboard_instance_observer_id",
        right_on="id",
        how="inner",
    )
    
    filtered_stops = (
        valid_raw_stops[valid_raw_stops["trip_id"].isin(segments_targets["trip_id"].unique())]
        .reset_index(drop=True)
        .copy()
    )
    filtered_stops["orig_index"] = filtered_stops.index
    
    def _snap_points_to_stop_pair_segments(group: gpd.GeoDataFrame) -> pd.DataFrame:
        """Return a DataFrame with segment_order/from_id/to_id/vehicle_name + orig_index for this trip."""
        segs = segments_targets[segments_targets["trip_id"] == group.name]
        if segs.empty:
            out = group[["orig_index"]].copy()
            out["segment_order"] = np.nan
            out["from_id"] = np.nan
            out["to_id"] = np.nan
            out["vehicle_name"] = np.nan
            return out[["segment_order", "from_id", "to_id", "vehicle_name", "orig_index"]]

        segs = segs.reset_index(drop=True)
        geoms = list(segs.geometry.values)
        tree = STRtree(geoms)
        geom_pos = {id(g): i for i, g in enumerate(geoms)}

        pos = []
        for pt in group.geometry.values:
            if pt is None or getattr(pt, "is_empty", False):
                pos.append(None)
                continue

            nearest_result = tree.nearest(pt)

            if nearest_result is None:
                pos.append(None)
            elif isinstance(nearest_result, (int, np.integer)):
                pos.append(int(nearest_result))
            else:
                pos.append(geom_pos.get(id(nearest_result)))

        # Build result
        rows = []
        for p in pos:
            if p is None:
                rows.append((np.nan, np.nan, np.nan, np.nan))
            else:
                r = segs.loc[p, ["segment_order", "from_id", "to_id", "vehicle_name"]]
                rows.append(tuple(r.values.tolist()))

        out = pd.DataFrame(rows, columns=["segment_order", "from_id", "to_id", "vehicle_name"])
        out["orig_index"] = group["orig_index"].values
        return out[["segment_order", "from_id", "to_id", "vehicle_name", "orig_index"]]

    matched_pairs = (
        filtered_stops.groupby("trip_id", group_keys=False)
        .apply(_snap_points_to_stop_pair_segments)
        .reset_index(drop=True)
    )

    # if feedback:
    #     feedback.pushInfo(f"raw stops rows: {len(valid_raw_stops)}")
    #     feedback.pushInfo(f"raw stops trip id: {valid_raw_stops["trip_id"]}")
    #     feedback.pushInfo(f"segments targets id: {segments_targets["trip_id"]}")
    #     feedback.pushInfo(f"segments columns: {matched_pairs.columns}")
    
    matched_stops = (
        matched_pairs.merge(
            filtered_stops,
            on="orig_index",
            how="left",
        )
        .groupby(
            ["onboard_instance_observer_id", "departed_at", "trip_id", "segment_order"],
            as_index=False,
        )
        .agg({"board": "sum", "alight": "sum", "created_at": "first"})
    )

    
    onboard_segments_with_occupancy = (
        segments_targets
        .merge(
            matched_stops[["trip_id", "onboard_instance_observer_id"]].drop_duplicates(),
            on="trip_id",
            how="inner",
        )
        .merge(
            matched_stops,
            on=["trip_id", "segment_order", "onboard_instance_observer_id"],
            how="left",
        )
        .sort_values(["trip_id", "onboard_instance_observer_id", "segment_order"])
        .fillna({"board": 0, "alight": 0})
        .assign(
            survey_departed_at=lambda df: df.groupby(
                ["trip_id", "onboard_instance_observer_id"]
            )["departed_at"].ffill().bfill()
        )
        .assign(
            survey_departed_at=lambda df: pd.to_datetime(df.survey_departed_at, utc=True)
            .dt.tz_localize(None)
        )
        .assign(
            survey_departed_at=lambda df: (
                df.survey_departed_at - df.survey_departed_at.dt.floor("d")
            ).dt.total_seconds()
        )
        .assign(
            vehicle_occupancy=lambda df: (
                df.groupby(["trip_id", "onboard_instance_observer_id"])
                .apply(lambda g: (g["board"] - g["alight"]).cumsum())
                .reset_index(level=[0, 1], drop=True)
            )
        )
    )
    
# ---- Build interval-specific start_time for each (trip_id, segment_order, interval_name) using od_stats.duration ----
    # Normalize od_stats columns: expect o_id, d_id, interval_name, vehicle_name, duration
    dur = od_stats_df.copy()
    if "o_id" not in dur.columns or "d_id" not in dur.columns:
        # od_stats in your RL2SDI is stop-pair based; if column names differ, adjust here
        raise ValueError("od_stats_df must include o_id and d_id columns.")
    if "interval_name" not in dur.columns:
        raise ValueError("od_stats_df must include interval_name (joinable to intervals_df).")
    if "vehicle_name" not in dur.columns:
        raise ValueError("od_stats_df must include vehicle_name for speed differentiation.")
    if "duration" not in dur.columns:
        raise ValueError("od_stats_df must include duration (seconds).")

    # Build a per-trip per-interval table of segment durations
    seg = trip_segments_gdf[["trip_id", "segment_order", "from_id", "to_id", "vehicle_name"]].copy()
    seg["key"] = 1
    i2 = intervals_df[["interval_name", "interval_start_secs", "interval_end_secs"]].copy()
    i2["key"] = 1
    seg_i = seg.merge(i2, on="key").drop(columns="key")

    dur2 = dur.rename(columns={"o_id": "from_id", "d_id": "to_id"})[
        ["from_id", "to_id", "interval_name", "vehicle_name", "duration"]
    ].copy()

    seg_i = seg_i.merge(
        dur2,
        on=["from_id", "to_id", "interval_name", "vehicle_name"],
        how="left"
    )

    # Fallback: median by vehicle_name, then global
    seg_i["duration"] = seg_i["duration"].fillna(
        seg_i.groupby("vehicle_name")["duration"].transform("median")
    )
    seg_i["duration"] = seg_i["duration"].fillna(seg_i["duration"].median())

    seg_i = seg_i.sort_values(["trip_id", "interval_name", "segment_order"])
    seg_i["start_time"] = seg_i.groupby(["trip_id", "interval_name"])["duration"].cumsum().shift(1).fillna(0)

    # Merge interval-specific start_time into occupancy rows
    occ = onboard_segments_with_occupancy.merge(
        seg_i[["trip_id", "segment_order", "interval_name", "start_time", "interval_start_secs", "interval_end_secs"]],
        on=["trip_id", "segment_order"],
        how="inner",
    )

    avg_occupancy_per_trip_segment_per_interval = (
        occ
        .query(
            "survey_departed_at + start_time < interval_end_secs & survey_departed_at + start_time > interval_start_secs"
        )
        .groupby(
            ["trip_id", "segment_order", "interval_name"],
            as_index=False,
        )
        .agg({"vehicle_occupancy": "median"})
        .rename(columns={"vehicle_occupancy": "vehicle_occupancy"})
    )

    return (
        avg_occupancy_per_trip_segment_per_interval,
        onboard_segments_with_occupancy,
        matched_stops,
        filtered_stops,
    )
