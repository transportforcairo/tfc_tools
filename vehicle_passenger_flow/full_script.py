# %%
import geopandas as gpd
import numpy as np
import pandas as pd
import os
from .utils import utils
from shapely.geometry import LineString, box, Point, MultiLineString, GeometryCollection
import json
from shapely.geometry import mapping

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

    # --- FIELD DEFINITIONS: infer type from pandas dtypes ---
    fields = []
    for col in gdf.columns:
        if col == gdf.geometry.name:
            continue
        qtype = _qvariant_type_for_series(gdf[col])

        fields.append(QgsField(str(col), qtype))
    pr.addAttributes(fields)

    vl.updateFields()

    # --- FEATURES / ATTRIBUTES: keep native types (no str() cast) ---
    for _, row in gdf.iterrows():
        feat = QgsFeature()
        feat.setGeometry(QgsGeometry.fromWkt(row.geometry.wkt))

        attrs = []
        for col in gdf.columns:
            if col == gdf.geometry.name:
                continue
            val = row[col]
            # Let QGIS handle None as NULL
            if pd.isna(val):
                attrs.append(None)
            else:
                attrs.append(val)

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


class FlowEstimator:
    def __init__(self, connection, output_folder, sdi_mode='postgres', gpkg_path=None, analysis_config=None):
        self.connection = connection
        self.output_folder = output_folder
        self.sdi_mode = sdi_mode
        self.gpkg_path = gpkg_path

        # Defaults (previously hard-coded)
        self.analysis_config = {
            "trip_segments_segmentization_threshold_meters": 300
        }

        # Allow algorithm UI to override.
        if isinstance(analysis_config, dict):
            self.analysis_config.update({k: analysis_config[k] for k in analysis_config if k in self.analysis_config})

    def download_required_layers(self, feedback=None):
        if feedback:
            feedback.pushInfo("Loading required SDI layers")

        if self.sdi_mode == 'postgres':
            # Postgres: use engine-based readers where available
            onboard_instances_df = pd.read_sql(
                "SELECT *, ST_AsText(geometry) AS geometry_wkt FROM raw.onboard_instances",
                con=self.connection
            )
            onboard_instances_df["geometry"] = gpd.GeoSeries.from_wkt(onboard_instances_df["geometry_wkt"])
            onboard_instances = gpd.GeoDataFrame(
                onboard_instances_df.drop(columns=["geometry_wkt"]),
                geometry="geometry",
                crs="EPSG:4326"
            )

            raw_stops = gpd.read_postgis(
                "SELECT * FROM raw.stops",
                con=self.connection,
                geom_col="geom"
            )
            raw_stops = raw_stops.rename(columns={raw_stops.geometry.name: "geometry"}).set_geometry("geometry")
            raw_stops = raw_stops.set_crs("EPSG:4326", allow_override=True)

            intervals = pd.read_sql("SELECT * FROM transit.intervals", con=self.connection)
            trips_intervals = pd.read_sql("SELECT * FROM transit.trips_intervals", con=self.connection)
            trip_stops_sequence = pd.read_sql("SELECT * FROM transit.trip_stops_sequence", con=self.connection)
            agencies = pd.read_sql("SELECT * FROM transit.agencies", con=self.connection)
            vehicles = pd.read_sql("SELECT * FROM transit.vehicles", con=self.connection)
            od_stats = pd.read_sql("SELECT * FROM transit.od_stats", con=self.connection)

            trips_view = gpd.read_postgis(
                "SELECT * FROM transit.trips_view",
                con=self.connection,
                geom_col="geom"
            )
            trips_view = trips_view.rename(columns={trips_view.geometry.name: "geometry"}).set_geometry("geometry")
            trips_view = trips_view.set_crs("EPSG:4326", allow_override=True)
            # After loading transit_trips_df (or trips_df)
            trip_xwalk = trips_view[["gid", "observer_id"]].drop_duplicates()
            trip_xwalk = trip_xwalk.rename(columns={"gid": "t_id", "observer_id": "trip_id"})
            trip_xwalk["trip_id"] = trip_xwalk["trip_id"].astype(str).str.strip()

        else:
            from ..tfc_tools_common.sdi_io import SDISource, read_df, read_gdf
            source = SDISource(mode='gpkg', gpkg_path=self.gpkg_path)

            onboard_instances = read_gdf(source, "raw.onboard_instances", geom_col="geometry")
            raw_stops = read_gdf(source, "raw.stops", geom_col="geom")
            raw_stops = raw_stops.rename(columns={raw_stops.geometry.name: "geometry"}).set_geometry("geometry")
            onboard_instances = onboard_instances.rename(columns={onboard_instances.geometry.name: "geometry"}).set_geometry("geometry")

            intervals = read_df(source, "transit.intervals")
            trips_intervals = read_df(source, "transit.trips_intervals")
            trip_stops_sequence = read_df(source, "transit.trip_stops_sequence")
            agencies = read_df(source, "transit.agencies")
            vehicles = read_df(source, "transit.vehicles")
            od_stats = read_df(source, "transit.od_stats")
            trips_view = read_gdf(source, "transit.trips_view", geom_col="geom")
            trips_view = trips_view.rename(columns={trips_view.geometry.name: "geometry"}).set_geometry("geometry")

        if feedback:
            feedback.pushInfo("Finished loading SDI layers.")
        
        # ============================================================
        # SDI NORMALIZATION: canonical trip_id = observer trip id (str)
        # Place this right after loading SDI layers in full_script.py
        # ============================================================

        def _norm_str(s: pd.Series) -> pd.Series:
            return s.astype(str).str.strip()

        # ---- 1) Normalize trips_view -> ensure trip_id exists and is canonical ----
        # Expect trips_view has observer_id (string id). Make it trip_id.
        if "trip_id" not in trips_view.columns:
            if "observer_id" in trips_view.columns:
                trips_view = trips_view.rename(columns={"observer_id": "trip_id"})
            else:
                raise ValueError("trips_view_gdf must have observer_id or trip_id")

        trips_view["trip_id"] = _norm_str(trips_view["trip_id"])

        # ---- 2) Normalize trip_stops_sequence -> ensure trip_id matches trips_view.trip_id ----
        # Prefer observer_trip_id if present (this is the matching universe).
        if "trip_id" not in trip_stops_sequence.columns:
            if "observer_trip_id" in trip_stops_sequence.columns:
                trip_stops_sequence = trip_stops_sequence.rename(columns={"observer_trip_id": "trip_id"})
            elif "trip_id" in trip_stops_sequence.columns:
                pass
            else:
                raise ValueError("trip_stops_sequence_df must have observer_trip_id or trip_id")

        trip_stops_sequence["trip_id"] = _norm_str(trip_stops_sequence["trip_id"])

        # Keep internal id as t_id if it exists, but never treat it as trip_id
        if "t_id" not in trip_stops_sequence.columns and "trip_id" in trip_stops_sequence.columns:
            # do nothing; just avoid renaming internal ids incorrectly
            pass

        # ---- 3) Trip crosswalk (if t_id exists) for any tables still using internal ids ----
        if "t_id" in trip_stops_sequence.columns:
            trip_xwalk = trip_stops_sequence[["t_id", "trip_id"]].drop_duplicates()

        # ---- 4) Normalize other tables that may carry trip identifiers ----
        # onboard_instances usually has trip_id already. Make it string canonical.
        if "trip_id" in onboard_instances.columns:
            onboard_instances["trip_id"] = _norm_str(onboard_instances["trip_id"])

        # ---- 5) Hard fail early if universes don't overlap ----
        common = set(trips_view["trip_id"]).intersection(set(trip_stops_sequence["trip_id"]))
        if len(common) == 0:
            raise ValueError(
                "No overlap between trips_view.trip_id and trip_stops_sequence.trip_id. "
                "You likely have mixed trip id universes (observer_id vs internal t_id)."
            )

        if feedback:
            feedback.pushInfo(f"Normalized trip_id. trips_view unique trips: {trips_view['trip_id'].nunique()}")
            feedback.pushInfo(f"Normalized trip_id. tss unique trips: {trip_stops_sequence['trip_id'].nunique()}")
            feedback.pushInfo(f"Trip_id overlap count: {len(common)}")

        return {
            "onboard_instances": onboard_instances,
            "raw_stops": raw_stops,
            "intervals": intervals,
            "trips_intervals": trips_intervals,
            "trip_stops_sequence": trip_stops_sequence,
            "agencies": agencies,
            "vehicles": vehicles,
            "od_stats": od_stats,
            "trips_view": trips_view,
        }
    

    def run(self,feedback=None): # we added feedback=None between brackets to know if connection ran successfully. and modified processAlgorithm
        if feedback:
            feedback.pushInfo("reading SDI data is starting.")
        # self.download_required_layers()
        logger.info("Downloading required layers...")
        logger.info("Loading SDI layers…")
        layers = self.download_required_layers(feedback)

        # CRS 3857 everywhere for metric work
        raw_stops_gdf = layers["raw_stops"].to_crs(3857)
        onboard_instances_gdf = layers["onboard_instances"].to_crs(3857)
        trips_view_gdf = layers["trips_view"].to_crs(3857)

        intervals_raw = layers["intervals"].copy()
        if "name" in intervals_raw.columns and "interval_name" not in intervals_raw.columns:
            intervals_raw = intervals_raw.rename(columns={"name": "interval_name"})
        if "active" in intervals_raw.columns:
            intervals_raw = intervals_raw[intervals_raw["active"].isin([True, 1, "1", "t", "T", "true", "TRUE"])].copy()

        # build interval seconds
        intervals_df = (
            intervals_raw.assign(
                interval_start_secs=lambda df: pd.to_timedelta(df["start_time"].astype(str)).dt.total_seconds(),
                interval_end_secs=lambda df: pd.to_timedelta(df["end_time"].astype(str)).dt.total_seconds(),
            )
            .sort_values("interval_start_secs")
            .reset_index(drop=True)
        )

        # ---- SDI-only: build trip stop-pair segments (replaces GTFS stop_times + shapes) ----
        trip_segments = utils.build_trip_stop_pair_segments_from_sdi(
            trips_view_gdf=trips_view_gdf,
            trip_stops_sequence_df=layers["trip_stops_sequence"],
            agencies_df=layers["agencies"],
            vehicles_df=layers["vehicles"],
        )

        # deduplicated stop-pair segments backbone for final outputs
        flow_stop_pair_segments = utils.build_flow_stop_pair_segments(trip_segments)

        # ---- Vehicle flow (supply) from SDI ----
        trips_intervals = layers["trips_intervals"].copy()
        if "headway" in trips_intervals.columns and "headway_secs" not in trips_intervals.columns:
            trips_intervals = trips_intervals.rename(columns={"headway": "headway_secs"})

        # resolve interval_id join
        intervals_join = intervals_raw.copy()
        if "gid" in intervals_join.columns and "interval_id" not in intervals_join.columns:
            intervals_join = intervals_join.rename(columns={"gid": "interval_id"})
        if "id" in intervals_join.columns and "interval_id" not in intervals_join.columns:
            intervals_join = intervals_join.rename(columns={"id": "interval_id"})

        intervals_join = intervals_join.merge(
            intervals_df[["interval_name", "interval_start_secs", "interval_end_secs"]],
            on="interval_name",
            how="left",
        )
        intervals_join["interval_duration_secs"] = intervals_join["interval_end_secs"] - intervals_join["interval_start_secs"]

        veh_supply = trips_intervals.merge(
            intervals_join[["interval_id", "interval_name", "interval_duration_secs"]],
            on="interval_id",
            how="inner",
        )

        if "trip_id" in veh_supply.columns:
            veh_supply["trip_id"] = veh_supply["trip_id"].astype(str).str.strip()
        else:
            raise ValueError("veh_supply must have trip_id (canonical) or t_id (mappable).")

        veh_supply["vehicle_trips_in_interval"] = np.floor(
            (veh_supply["interval_duration_secs"] / veh_supply["headway_secs"]).replace([np.inf, -np.inf], np.nan).fillna(0)
        )
        

        # ---- Occupancy using SDI segments + time proxy from transit.od_stats.duration ----
        od_stats_df = layers["od_stats"].copy()

        avg_occ, onboard_segments_with_occupancy, matched_stops, filtered_stops = utils.get_avg_occupancy_per_segment_v3_sdi_timeproxy(
            trip_segments_gdf=trip_segments,
            intervals_df=intervals_df,
            raw_stops_gdf=raw_stops_gdf,
            onboard_instances_gdf=onboard_instances_gdf,
            od_stats_df=od_stats_df,
            feedback=feedback
        )

        
        trip_xwalk = trips_view_gdf[["gid", "trip_id"]].drop_duplicates()
        trip_xwalk = trip_xwalk.rename(columns={"gid": "t_id"})
        trip_xwalk["t_id"] = trip_xwalk["t_id"].astype(str).str.strip()        # Interpret current trip_segments.trip_id as internal id (gid)

        # 3) Replace with observer trip_id
        veh_supply = veh_supply.rename(columns={"trip_id": "t_id"})
        veh_supply["t_id"] = veh_supply["t_id"].astype(str).str.strip()
        veh_supply = (
            veh_supply
            .merge(trip_xwalk[["trip_id","t_id"]], on="t_id", how="left")
        )

        # 4) Fail fast if mapping didn’t work
        # if trip_segments["trip_id"].isna().any():
        #     raise ValueError("trip_segments: gid->observer_id mapping produced NaNs. Check that trip_segments.trip_id truly holds gid.")--
        
        disagg = (
            trip_segments[["trip_id", "segment_order", "from_id", "to_id", "vehicle_name", "geometry"]]
            .merge(veh_supply[["trip_id", "interval_name", "vehicle_trips_in_interval"]], on="trip_id", how="left")
            .merge(
                avg_occ.rename(columns={"vehicle_occupancy": "occ_median"}),
                on=["trip_id", "segment_order", "interval_name"],
                how="left",
            )
        )

        disagg["vehicle_flow"] = disagg["vehicle_trips_in_interval"].fillna(0)
        # conservative occupancy fallback (same spirit as your existing behavior)
        disagg["occ_median"] = disagg["occ_median"].fillna(disagg["occ_median"].median())
        disagg["occ_median"] = disagg["occ_median"].clip(lower=0)
        disagg["passenger_flow"] = np.ceil(disagg["vehicle_flow"] * disagg["occ_median"]).fillna(0)

        disagg_gdf = gpd.GeoDataFrame(disagg, geometry="geometry", crs="EPSG:3857")

        # ---- Aggregate to stop-pair × interval (TALL) ----
        tall = (
            disagg_gdf.groupby(["from_id", "to_id", "interval_name"], as_index=False)[["vehicle_flow", "passenger_flow"]]
            .sum()
        )

        flow_stop_pair_segments_overall = (
            flow_stop_pair_segments
            .drop(columns=["vehicle_name", "flow_seg_id"], errors="ignore")
            .drop_duplicates(subset=["from_id", "to_id"])
            .copy()
        )

        flow_stop_pair_segments_overall["flow_seg_id"] = (
        flow_stop_pair_segments_overall["from_id"].astype(str)
        + "__"
        + flow_stop_pair_segments_overall["to_id"].astype(str)
        )

        tall_gdf = flow_stop_pair_segments_overall.merge(
            tall,
            on=["from_id", "to_id"],
            how="left"
        ).fillna(0)
        tall_gdf = gpd.GeoDataFrame(tall_gdf, geometry="geometry", crs="EPSG:3857")

        # ---- Pivot to WIDE (veh__<interval>, pax__<interval>) ----
        def _safe(s):
            return (
                str(s).strip().lower()
                .replace(" ", "_")
                .replace("-", "_")
                .replace("/", "_")
            )

        wide = flow_stop_pair_segments_overall.copy()
        for iname in intervals_df["interval_name"].unique():
            key = _safe(iname)
            sub = tall_gdf[tall_gdf["interval_name"] == iname][["from_id", "to_id", "vehicle_flow", "passenger_flow"]]
            wide = wide.merge(
                sub.rename(columns={
                    "vehicle_flow": f"veh__{key}",
                    "passenger_flow": f"pax__{key}",
                }),
                on=["from_id", "to_id"],
                how="left",
            )
        wide = wide.fillna(0)
        wide_gdf = gpd.GeoDataFrame(wide, geometry="geometry", crs="EPSG:3857")

        # ---- Write ONE output GeoPackage with ALL layers ----
        out_gpkg = os.path.join(self.output_folder, "vehicle_passenger_flow.gpkg")

        save_gdf_with_qgis_writer(flow_stop_pair_segments_overall, out_gpkg, "analysis.flow_stop_pair_segments", feedback)
        save_gdf_with_qgis_writer(disagg_gdf, out_gpkg, "analysis.flow_disagg_trip_interval", feedback)
        save_gdf_with_qgis_writer(tall_gdf, out_gpkg, "analysis.flow_stop_pair_interval_tall", feedback)
        save_gdf_with_qgis_writer(wide_gdf, out_gpkg, "analysis.flow_stop_pair_interval_wide", feedback)



        # logger.success("Vehicle and passenger flow estimation complete.") # WE REPLACED logger.success() WITH THE NEXT LINE
        logger.info("✅ SUCCESS: Vehicle and passenger flow estimation complete.")