# %%
import geopandas as gpd
import numpy as np
import pandas as pd
import os
from .utils import utils
from shapely.geometry import LineString, box, Point, MultiLineString, GeometryCollection
import json
from shapely.geometry import mapping
import zipfile

# import sys
# plugin_dir = os.path.dirname(__file__)
# libs_path = os.path.join(plugin_dir, "libs")
# if libs_path not in sys.path:
#     sys.path.insert(0, libs_path)

# from loguru import logger
import logging
logger = logging.getLogger("vehicle_passenger_flow") 
logging.basicConfig(level=logging.INFO) # This ensures that .info(), .debug(), .error() will actually output something.

import math
# from shapely.geometry import LineString
def discrete_frechet(P, Q):
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

# --- QGIS writer: save GeoPandas GDF to disk without fiona/pyogrio ---
from qgis.core import (
    QgsVectorLayer,
    QgsField,
    QgsFeature,
    QgsGeometry,
    QgsFields,
    QgsCoordinateTransformContext,
    QgsVectorFileWriter,
    QgsWkbTypes,
    QgsProject
)
from qgis.PyQt.QtCore import QVariant


def _qvariant_type_for_series(s: pd.Series):
    if pd.api.types.is_integer_dtype(s):
        return QVariant.LongLong
    if pd.api.types.is_float_dtype(s):
        return QVariant.Double
    if pd.api.types.is_bool_dtype(s):
        return QVariant.Bool
    return QVariant.String

def _wkb_from_geom_type(geom_type_str: str):
    # Map common GeoPandas geom types to QGIS WKB
    g = geom_type_str.lower()
    if "multilinestring" in g: return QgsWkbTypes.MultiLineString
    if "linestring"      in g: return QgsWkbTypes.LineString
    if "multipolygon"    in g: return QgsWkbTypes.MultiPolygon
    if "polygon"         in g: return QgsWkbTypes.Polygon
    if "multipoint"      in g: return QgsWkbTypes.MultiPoint
    if "point"           in g: return QgsWkbTypes.Point
    return QgsWkbTypes.Unknown


from qgis.core import QgsCoordinateReferenceSystem
def gdf_to_qgis_layer(gdf, layer_name):
    """Convert a GeoDataFrame to a QGIS memory layer."""
    # Create memory layer with the same geometry type and CRS
    geom_type = gdf.geometry.iloc[0].geom_type if not gdf.empty else "Point"

    # crs = gdf.crs.to_wkt() if gdf.crs else "EPSG:4326"
    # vl = QgsVectorLayer(f"{geom_type}?crs={crs}", layer_name, "memory")

    crs_obj = QgsCoordinateReferenceSystem()
    if gdf.crs:
        crs_obj.createFromWkt(gdf.crs.to_wkt())
    else:
        crs_obj.createFromEpsgId(4326)
    
    vl = QgsVectorLayer(f"{geom_type}?crs={crs_obj.authid()}", layer_name, "memory")
    
    # Add fields
    pr = vl.dataProvider()
    pr.addAttributes([QgsField(str(col), QVariant.String) for col in gdf.columns if col != gdf.geometry.name])
    vl.updateFields()
    
    # Add features
    for _, row in gdf.iterrows():
        feat = QgsFeature()
        feat.setGeometry(QgsGeometry.fromWkt(row.geometry.wkt))
        attrs = [str(row[col]) for col in gdf.columns if col != gdf.geometry.name]
        feat.setAttributes(attrs)
        pr.addFeature(feat)
    
    vl.updateExtents()
    return vl


def _interpret_writer_result(res):
    """Return (err, new_file, new_layer) no matter the QGIS version."""
    if isinstance(res, tuple):
        # Newer PyQGIS often returns (err, newFileName, newLayerName)
        if len(res) == 3:
            return res
        # Some variants return 4-tuple; be defensive
        if len(res) >= 1:
            return (res[0], None, None)
        return (res, None, None)
    else:
        # Older style: just the error enum
        return (res, None, None)
    

def save_gdf_with_qgis_writer(gdf, out_path, layer_name, feedback=None):
    # quick exits
    if gdf is None or gdf.empty:
        if feedback: feedback.pushInfo(f"[INFO] '{layer_name}': nothing to write.")
        return

    # make sure folder exists
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    # build an in‑memory layer from the gdf (your existing helper)
    vl = gdf_to_qgis_layer(gdf, layer_name)

    # writer options
    options = QgsVectorFileWriter.SaveVectorOptions()
    options.driverName = "GPKG" if out_path.lower().endswith(".gpkg") else "GeoJSON"
    options.layerName = layer_name

    # choose action based on file existence
    if out_path.lower().endswith(".gpkg"):
        if os.path.exists(out_path):
            # overwrite/replace this layer inside the existing gpkg
            options.actionOnExistingFile = QgsVectorFileWriter.CreateOrOverwriteLayer
        else:
            # create a new gpkg file
            options.actionOnExistingFile = QgsVectorFileWriter.CreateOrOverwriteFile

    # try write
    transform_ctx = QgsProject.instance().transformContext() or QgsCoordinateTransformContext()
    res = QgsVectorFileWriter.writeAsVectorFormatV3(vl, out_path, transform_ctx, options)

    err, new_file, new_layer = _interpret_writer_result(res)

    if err != QgsVectorFileWriter.NoError:
        # retry logic for GPKG only
        if out_path.lower().endswith(".gpkg") and os.path.exists(out_path):
            try:
                if feedback: feedback.pushInfo(f"[WARN] Writer error {err}. Retrying after removing {out_path}…")
                os.remove(out_path)
                options.actionOnExistingFile = QgsVectorFileWriter.CreateOrOverwriteFile
                res2 = QgsVectorFileWriter.writeAsVectorFormatV3(vl, out_path, transform_ctx, options)
                err2, _, _ = _interpret_writer_result(res2)
                if err2 != QgsVectorFileWriter.NoError:
                    raise RuntimeError(f"Failed writing '{layer_name}' to {out_path}. Writer code: {res2}")
                if feedback: feedback.pushInfo(f"[OK] Wrote layer '{layer_name}' to {out_path} (fresh file).")
                return
            except Exception as e:
                raise RuntimeError(f"Failed writing '{layer_name}' to {out_path} after retry: {e}")
        else:
            raise RuntimeError(f"Failed writing '{layer_name}' to {out_path}. Writer code: {res}")
    else:
        # success
        if feedback:
            # prefer the returned file/layer names when present
            final_file = new_file or out_path
            final_layer = new_layer or layer_name
            feedback.pushInfo(f"[OK] Wrote layer '{final_layer}' to {final_file}.")
# --- end of QGIS writer: save GeoPandas GDF to disk without fiona/pyogrio ---

def get_boundary_from_gtfs(gtfs_zip_path, feedback=None):
    try:
        if feedback:
            feedback.pushInfo("Extracting boundary from GTFS zip...")

        with zipfile.ZipFile(gtfs_zip_path, 'r') as z:
            with z.open("stops.txt") as f:
                stops_df = pd.read_csv(f)

        # Ensure lat/lon exist
        if not {"stop_lat", "stop_lon"}.issubset(stops_df.columns):
            raise ValueError("GTFS stops.txt is missing stop_lat or stop_lon")

        # Convert to GeoDataFrame
        stops_gdf = gpd.GeoDataFrame(
            stops_df,
            geometry=[Point(xy) for xy in zip(stops_df["stop_lon"], stops_df["stop_lat"])],
            crs="EPSG:4326"
        )

        # Option 1: bounding box with small buffer
        bounds = stops_gdf.total_bounds
        boundary = box(*bounds).buffer(0.01)

        # Option 2 (optional): convex hull
        # boundary = stops_gdf.unary_union.convex_hull.buffer(0.005)

        if feedback:
            feedback.pushInfo("Boundary extracted from GTFS successfully.")

        return boundary

    except Exception as e:
        if feedback:
            feedback.pushInfo(f"Failed to extract boundary from GTFS: {e}")
        raise

# sys.path.insert(1, "./utils/utils.py")

# THIS PART IS FROM local_map_data_extractor PLUGIN FOR PIS
import requests
# Function to extract road data within the polygon layer's extent using Overpass API
def extract_osm_roads(self, map_boundary, feedback=None):
    if feedback:
            feedback.pushInfo("Loading map boundary...")

    # Ensure it’s in EPSG:4326 (WGS84 lat/lon)
    if map_boundary.crs != 'EPSG:4326':
        map_boundary = map_boundary.to_crs('EPSG:4326')

    if feedback:
            feedback.pushInfo("map boundary correct CRS is finished.")

    # Calculate the extent (bounding box) of the entire polygon layer
    bounds_polygon = map_boundary.total_bounds

    # Combine the bounds to create a new bounding box
    south = bounds_polygon[1]
    west = bounds_polygon[0]
    north = bounds_polygon[3]
    east = bounds_polygon[2]

    if feedback:
        feedback.pushInfo("Bounding box for OSM query calculated.")

    # --- THIS SECTION IS FOR DEBUGGING BY EXPORTING THE BOUNDING BOX IN 4326
    # if self.output_folder and isinstance(self.output_folder, str):
    #     try:
    #         debug_path = os.path.join(self.output_folder, "debug_boundary.geojson")
    #         # Write a minimal GeoJSON FeatureCollection
    #         with open(debug_path, "w", encoding="utf-8") as f:
    #             geojson_dict = {
    #                 "type": "FeatureCollection",
    #                 "features": [
    #                     {
    #                         "type": "Feature",
    #                         "geometry": mapping(map_boundary.geometry.iloc[0]),
    #                         "properties": {},
    #                     }
    #                 ],
    #             }
    #             json.dump(geojson_dict, f)
    #         if feedback:
    #             feedback.pushInfo(f"Saved debug boundary to {debug_path}")
    #     except Exception as e:
    #         if feedback:
    #             feedback.pushInfo(f"Failed to export debug boundary: {e}")


    # Construct Overpass API query
    overpass_url = "http://overpass-api.de/api/interpreter"
    overpass_query = f"""
    [out:json];
    (
    way["highway"]({south},{west},{north},{east});
    );
    out body;
    >;
    out skel qt;
    """

    # Fetch data from Overpass API
    response = requests.get(overpass_url, params={'data': overpass_query})
    # response.raise_for_status() # a line by ChatGPT that was not in local_map_data_extractor
    if response.status_code != 200: # this part is by ChatGPT similar to the line above
        if feedback:
            feedback.pushInfo(f"Error: Overpass returned status code {response.status_code}")
        return gpd.GeoDataFrame(columns=["road_type", "road_name", "geometry"], crs="EPSG:4326")

    data = response.json()

    if feedback:
            feedback.pushInfo("fetch data from Overpass API is finished.")

    # Parse the data to extract nodes and ways
    elements = data['elements']
    nodes = {el['id']: (el['lon'], el['lat']) for el in elements if el['type'] == 'node'}
    ways = [el for el in elements if el['type'] == 'way']

    if feedback:
            feedback.pushInfo("parse data to extract nodes and ways is finished.")

    edges = []

    # Iterate over ways to process edges
    for way in ways:
        if 'tags' in way and 'highway' in way['tags']:
            road_type = way['tags'].get('highway', 'Unknown')
            road_name = way['tags'].get('name', None) # A NEW LINE BY CHATGPT THAT WAS NOT IN local_map_data_extractor plugin

            coords = []
            for node_id in way['nodes']:
                if node_id in nodes:
                    coords.append(nodes[node_id])
            
            if len(coords) >= 2:
                geometry = LineString(coords)
                edges.append({
                            'gid': len(edges),  # unique ID for each row
                            'road_type': road_type,
                            'road_name': road_name,
                            'geometry': geometry
                        })  # NEW PART AS A REPLACEMENT TO THE LINES BELOW

    if feedback:
            feedback.pushInfo("Converted ways to LineString geometrie.")

    feedback.pushInfo(f"Number of edges collected: {len(edges)}")
    if len(edges) > 0:
        feedback.pushInfo(f"Sample edge: {edges[0]}")

    # Convert edges to a GeoDataFrame
    gdf_edges = gpd.GeoDataFrame(edges, crs="EPSG:4326", geometry='geometry')


    # Reproject polygon to match edges CRS (EPSG:4326)
    map_boundary = map_boundary.to_crs(gdf_edges.crs)

    # Clip roads to the actual polygon boundary (not just bounding box)
    gdf_edges_clipped = gpd.clip(gdf_edges, map_boundary)
    
    print(f"Clipped road data To boundary")

    return gdf_edges_clipped


class FlowEstimator:
    def __init__(self, gtfs_zip_path, connection, output_folder):
        self.gtfs_zip_path = gtfs_zip_path
        self.connection = connection
        self.output_folder = output_folder

        self.analysis_config = {
            "trip_segments_segmentization_threshold_meters": 300,
            "segment_matching_buffer_meters": 10,
            "frechet_dist_densify_param": 0.1,
            "frechet_dist_segment_length_ratio_max_thresh": 0.5, #revise this part
        }

        self.interval_boundaries = [
            ("06:00:00", "10:00:00"), #revise timezones later
            ("10:05:00", "15:45:00")
        ]        

    def download_required_layers(self, feedback=None):
        if feedback:
            feedback.pushInfo("downloading required layers is starting.")
        # os.makedirs(self.output_folder, exist_ok=True) # REMOVE this if not exporting files anymore

        # 1.1 Load raw.onboard_instances with WKT instead of native PostGIS geometry
        onboard_instances_df = pd.read_sql(
            "SELECT *, ST_AsText(geometry) AS geometry_wkt FROM raw.onboard_instances",
            con=self.connection
        )

        if feedback:
            feedback.pushInfo("Loading raw.onboard_instances is finished.")

        # Replace geometry column with one parsed from WKT, keep column name as 'geometry'
        onboard_instances_df["geometry"] = gpd.GeoSeries.from_wkt(onboard_instances_df["geometry_wkt"])
        onboard_instances = gpd.GeoDataFrame(
            onboard_instances_df.drop(columns=["geometry_wkt"]),
            geometry="geometry",
            crs="EPSG:4326"
        )

        if feedback:
            feedback.pushInfo("modifying raw.onboard_instances is finished.")

        # 1.2 Load raw.stops and join with raw.stops_created_at
        raw_stops = gpd.read_postgis(
            "SELECT *, geom AS geometry FROM raw.stops",  # I later added: mgeom AS geometry
            con=self.connection,
            geom_col="geometry"
        )
        stops_created_at = pd.read_sql("SELECT observer_id, created_at FROM raw.stops_created_at", con=self.connection)
        raw_stops = raw_stops.merge(stops_created_at, on="observer_id", how="left")

        if feedback:
            feedback.pushInfo("Loading raw.stops is finished.")

        # 2.1 Create bounding box
        boundary = get_boundary_from_gtfs(self.gtfs_zip_path, feedback)
        boundary_gdf = gpd.GeoDataFrame(geometry=[boundary], crs="EPSG:4326")

        if feedback:
            feedback.pushInfo("creating bounding box is finished.")

        # 2.2 Download OSM road network
        # run the function
        road_gdf = extract_osm_roads(self, boundary_gdf, feedback)

        if feedback:
            feedback.pushInfo("Loading OSM road network is finished")

        return onboard_instances, raw_stops, road_gdf
    

    def run(self,feedback=None): # we added feedback=None between brackets to know if connection ran successfully. and modified processAlgorithm
        if feedback:
            feedback.pushInfo("reading SDI data is starting.")
        # self.download_required_layers()
        logger.info("Downloading required layers...")
        # onboard_instances, raw_stops, road_gdf = self.download_required_layers()
        onboard_instances, raw_stops, road_gdf = self.download_required_layers(feedback)

        if feedback:
            feedback.pushInfo("reading SDI data is finished.")

        # Convert all to CRS 3857
        road_segments_layer = road_gdf.to_crs(3857)
        raw_stops_gdf = raw_stops.to_crs(3857)
        onboard_instances_gdf = onboard_instances.to_crs(3857)

        if feedback:
            feedback.pushInfo("converting all layers to CRS 3857 is finished.")


        logger.info("Extracting vehicle appearances")
        (
            veh_app,
            transport_model,
            candidates_base,
            candidate_trip_segments,
            matched_trip_segments,
            trips_segments_gdf,
            stops,
            trip_instances,
            temp2
        ) = utils.get_vehicle_appearances(

            self.gtfs_zip_path,
            road_segments_gdf=road_segments_layer,
            config=self.analysis_config,
        )

        all_vehicle_occurrences = gpd.GeoDataFrame(  #debug this function
            pd.concat(
                [
                    veh_app.assign(gtfs="value"),
                ]
            )
        )

        logger.info("Generating time intervals")
        intervals_df = (
            pd.DataFrame(
                self.interval_boundaries,
                columns=["start", "end"],
            )
            .assign(
                interval_start_secs=lambda df: (
                    df.start.astype("datetime64[ns]")
                    - df.start.astype("datetime64[ns]").dt.floor("d")
                ).dt.total_seconds()
            )
            .assign(
                interval_end_secs=lambda df: (
                    df.end.astype("datetime64[ns]")
                    - df.end.astype("datetime64[ns]").dt.floor("d")
                ).dt.total_seconds()
            )
            .assign(
                interval_name=[
                    "morning_peak",
                    "afternoon"
                ]
            )
        )

        logger.info("Classifying vehicle occurrences into intervals")
        # Remove geometry for compatibility
        veh_occ_df = all_vehicle_occurrences.drop(columns="geometry")

        # Cartesian join and filter manually
        veh_occ_df["key"] = 1
        intervals_df["key"] = 1

        merged_df = pd.merge(veh_occ_df, intervals_df, on="key").drop(columns="key")
        all_vehicle_occurrences_intervals = merged_df[
            (merged_df["timeofday_secs"] > merged_df["interval_start_secs"]) &
            (merged_df["timeofday_secs"] < merged_df["interval_end_secs"])
        ].copy()



        logger.info("Generating vehicle flow layers")
        veh_flow = (
            road_segments_layer[["gid", "geometry"]]
            .merge(
                all_vehicle_occurrences_intervals.groupby(
                    ["gid", "interval_name", "gtfs"], as_index=False
                )["timeofday_secs"]
                .count()
                .pivot_table(
                    values="timeofday_secs", columns="gtfs", index=["gid", "interval_name"]
                )
                .reset_index()
                .reset_index(drop=True),
                on="gid",
                how="left",
            )
            .fillna(0)
        )


        logger.info("Estimating average vehicle occupancy")
        trips_segments_gdf, stops = utils.get_trips_segments_from_gtfs(self.gtfs_zip_path)

        (
            avg_occupancy_per_trip_segment_per_interval,
            onboard_segments_with_occupancy,
            matched_stops,
            filtered_stops,
        ) = utils.get_avg_occupancy_per_segment_v2(
            trips_segments_gdf,
            intervals_df=intervals_df,
            stops_gdf=stops,
            raw_stops_gdf=raw_stops_gdf,
            onboard_instances_gdf=onboard_instances_gdf,
        )

        trips_gdf = utils.read_trips_gdf_from_gtfs(self.gtfs_zip_path)

        logger.info("Interpolating for trips without onboard data")
        trips_without_onboard_survey = trips_gdf[
            ~trips_gdf["trip_id"].isin(
                avg_occupancy_per_trip_segment_per_interval["trip_id"].unique()
            )
        ]

        segments_of_trips_without_onboard_surveys = trips_segments_gdf[
            ~trips_segments_gdf["trip_id"].isin(
                avg_occupancy_per_trip_segment_per_interval["trip_id"].unique()
            )
        ].assign(geometry=lambda gdf: gdf.geometry.simplify(100))

        segments_of_trips_with_onboard_surveys = trips_segments_gdf[
            trips_segments_gdf["trip_id"].isin(
                avg_occupancy_per_trip_segment_per_interval["trip_id"].unique()
            )
        ].assign(geometry=lambda gdf: gdf.geometry.simplify(100))

        # --- helper: reduce any geometry to a single LineString (longest component) ---
        def _to_single_linestring(geom):
            if geom is None or geom.is_empty:
                return None
            if isinstance(geom, LineString):
                return geom
            if isinstance(geom, MultiLineString):
                if not geom.geoms:
                    return None
                return max(geom.geoms, key=lambda g: g.length)
            if isinstance(geom, GeometryCollection):
                # prefer LineStrings
                lines = [g for g in geom.geoms if isinstance(g, LineString)]
                if lines:
                    return max(lines, key=lambda g: g.length)
                # fallback: longest MultiLineString sub-geom
                mls = [g for g in geom.geoms if isinstance(g, MultiLineString)]
                if mls:
                    return max(mls, key=lambda g: g.length)
                return None
            # other geom types (Point/Polygon) aren't usable here
            return None

        # --- build the candidates (same spatial logic, safer geometry handling) ---
        joined = (
            gpd.overlay(
                segments_of_trips_without_onboard_surveys[["trip_id", "segment_order", "geometry"]]
                .assign(geom_buffered=lambda gdf: gdf.geometry.buffer(100))
                .set_geometry("geom_buffered"),
                segments_of_trips_with_onboard_surveys[["trip_id", "segment_order", "geometry"]],
                keep_geom_type=False,
            )
            .merge(
                trips_segments_gdf[["trip_id", "segment_order", "geometry"]],
                left_on=["trip_id_2", "segment_order_2"],
                right_on=["trip_id", "segment_order"],
                suffixes=["_1", "_2"],
            )
        )

        # compute a safe Frechet distance column (as a 1D Series) + lengths
        def _safe_frechet(row):
            g1 = _to_single_linestring(row.get("geometry_1"))
            g2 = _to_single_linestring(row.get("geometry_2"))
            if g1 is None or g2 is None:
                return np.nan
            try:
                P = list(g1.coords)
                Q = list(g2.coords)
                if not P or not Q:
                    return np.nan
                # try both directions to be orientation-agnostic
                return min(
                    discrete_frechet(P, Q),
                    discrete_frechet(list(reversed(P)), Q)
                )
            except Exception:
                return np.nan

        # create columns one by one (avoid assign/apply DataFrame-shape pitfalls)
        vals = [ _safe_frechet(row) for _, row in joined.iterrows() ]
        joined["frechet_dist"] = pd.to_numeric(pd.Series(vals, index=joined.index), errors="coerce")


        # length of the 'from' geometry used in the ratio
        joined["segment_length_meters"] = joined["geometry_1"].apply(
            lambda g: (_to_single_linestring(g).length if _to_single_linestring(g) is not None else np.nan)
        )

        # acceptance rule (guard against division by zero/NaN)
        joined["accepted"] = (joined["frechet_dist"] / joined["segment_length_meters"]) < 0.7




        # continue with your original flow
        joined_top_matches = (
            joined.groupby(["trip_id_1", "trip_id_2"], as_index=False)["segment_order_2"]
            .count()
            .rename(columns={"segment_order_2": "n_matched_segments"})
            .sort_values(["trip_id_1", "n_matched_segments"], ascending=False)
            .groupby("trip_id_1")
            .head(4)
        )

        joined_filtered = joined.merge(joined_top_matches, on=["trip_id_1", "trip_id_2"])

        # all_vehicle_occurrences_intervals.interval_start_secs.unique()

        avg_occupancy_per_trip_segment_per_interval_interpolated = (
            joined_filtered[
                ["trip_id_1", "segment_order_1", "trip_id_2", "segment_order_2"]
            ]
            .merge(
                avg_occupancy_per_trip_segment_per_interval,
                left_on=["trip_id_2", "segment_order_2"],
                right_on=["trip_id", "segment_order"],
            )
            .groupby(
                [
                    "trip_id_1",
                    "segment_order_1",
                    "interval_name",
                    "interval_start_secs",
                    "interval_end_secs",
                ]
            )
            .agg({"vehicle_occupancy": "median"})
            .reset_index()
            .rename(columns={"trip_id_1": "trip_id", "segment_order_1": "segment_order"})
        )

        avg_occupancy_per_trip_segment_per_interval_all = pd.concat(
            [
                avg_occupancy_per_trip_segment_per_interval.assign(interpolated=False),
                avg_occupancy_per_trip_segment_per_interval_interpolated.assign(
                    interpolated=True
                ),
            ]
        )

        logger.info("Joining occupancy data to vehicle appearances")
        vehicle_appearances_with_avg_occupancy = all_vehicle_occurrences_intervals.merge(
            avg_occupancy_per_trip_segment_per_interval_all.assign(
                vehicle_occupancy=lambda df: df.vehicle_occupancy.clip(lower=0)
            ),
            on=[
                "trip_id",
                "segment_order",
                "interval_name",
                "interval_start_secs",
                "interval_end_secs",
            ],
            how="left",
        )

        vehicle_appearances_with_avg_occupancy.loc[
            vehicle_appearances_with_avg_occupancy['gtfs']=='value', 'vehicle_occupancy'
        ] = vehicle_appearances_with_avg_occupancy.query("gtfs=='value'").assign(
            vehicle_occupancy=lambda df: df.vehicle_occupancy.fillna(df.vehicle_occupancy.median())
        )['vehicle_occupancy']

        vehicle_appearances_with_avg_occupancy = vehicle_appearances_with_avg_occupancy.assign(
            interpolated= lambda df: df.interpolated.fillna(True)
        )

        logger.info("Generating passenger flow layers")
        passenger_flow = (
            road_segments_layer[["gid", "geometry"]]
            .merge(
                vehicle_appearances_with_avg_occupancy.groupby(
                    ["gid", "interval_name", "gtfs"], as_index=False
                )["vehicle_occupancy"]
                .sum()
                .assign(vehicle_occupancy= lambda df: df.vehicle_occupancy.apply(np.ceil))
                .pivot_table(
                    values="vehicle_occupancy",
                    columns="gtfs",
                    index=["gid", "interval_name"]
                )
                .reset_index()
                .reset_index(drop=True),
                on="gid",
                how="left",
            )
            .fillna(0)
        )

        veh_gpkg = os.path.join(self.output_folder, "veh_flow.gpkg")
        save_gdf_with_qgis_writer(
            veh_flow.query("interval_name == 'morning_peak'").set_crs(3857, inplace=False),
            veh_gpkg, 
            "morning_peak", 
            feedback
            )
        save_gdf_with_qgis_writer(
            veh_flow.query("interval_name == 'afternoon'").set_crs(3857, inplace=False),  
            veh_gpkg, 
            "afternoon",   
            feedback
            )

        pass_gpkg = os.path.join(self.output_folder, "passenger_flow.gpkg")
        save_gdf_with_qgis_writer(
            passenger_flow.query("interval_name == 'morning_peak'").set_crs(3857, inplace=False),
            pass_gpkg, 
            "morning_peak", 
            feedback
            )
        save_gdf_with_qgis_writer(
            passenger_flow.query("interval_name == 'afternoon'").set_crs(3857, inplace=False),   
            pass_gpkg, 
            "afternoon",   
            feedback
            )


        for name, df in [
            ("veh_flow total", veh_flow),
            ("veh morning", veh_flow.query("interval_name == 'morning_peak'")),
            ("veh afternoon", veh_flow.query("interval_name == 'afternoon'")),
        ]:
            if feedback:
                geom_ok = 0 if df.empty else (~df.geometry.is_empty).sum()
                feedback.pushInfo(f"[DEBUG] {name}: rows={len(df)}, non-empty geoms={geom_ok}, crs={df.crs}")


        # logger.success("Vehicle and passenger flow estimation complete.") # WE REPLACED logger.success() WITH THE NEXT LINE
        logger.info("✅ SUCCESS: Vehicle and passenger flow estimation complete.")