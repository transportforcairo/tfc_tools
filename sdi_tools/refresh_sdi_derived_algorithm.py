# -*- coding: utf-8 -*-

from __future__ import annotations

import os
import sqlite3
import math
from collections import Counter
from dataclasses import dataclass
from typing import Iterable, Optional

import pandas as pd
import geopandas as gpd
from shapely.geometry import LineString, Point
from shapely.ops import substring

from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import (
    QgsProcessingAlgorithm,
    QgsProcessingParameterEnum,
    QgsProcessingParameterProviderConnection,
    QgsProcessingParameterFile,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterDefinition,
)

from ..tfc_tools_common import ensure_paths, ensure_deps
from ..tfc_tools_common.sdi_io import SDISource, postgres_engine_from_qgis_connection

# for icon
from qgis.PyQt.QtGui import QIcon

def _icon_path(*parts):
    # __file__ is inside tfc_tools/rl2sdi/
    here = os.path.dirname(__file__)
    root = os.path.abspath(os.path.join(here, ".."))     # tfc_tools/
    return os.path.join(root, *parts)

ensure_paths()


def _sqlite_type_from_series(s: pd.Series) -> str:
    if pd.api.types.is_integer_dtype(s):
        return "INTEGER"
    if pd.api.types.is_float_dtype(s):
        return "REAL"
    if pd.api.types.is_bool_dtype(s):
        return "INTEGER"
    return "TEXT"


def _write_df(conn: sqlite3.Connection, name: str, df: pd.DataFrame):
    df = df.copy()
    if "gid" in df.columns:
        df["gid"] = pd.to_numeric(df["gid"], errors="coerce")
        if df["gid"].isna().any():
            raise ValueError(f"Table '{name}': existing 'gid' contains null/non-numeric values.")
        if df["gid"].duplicated().any():
            raise ValueError(f"Table '{name}': existing 'gid' is not unique.")
        df["gid"] = df["gid"].astype("int64")
    else:
        df.insert(0, "gid", range(1, len(df) + 1))

    cols_sql = []
    for col in df.columns:
        if col == "gid":
            cols_sql.append('"gid" INTEGER PRIMARY KEY')
        else:
            cols_sql.append(f'"{col}" {_sqlite_type_from_series(df[col])}')

    conn.execute(f'DROP TABLE IF EXISTS "{name}"')
    conn.execute(f'CREATE TABLE "{name}" ({", ".join(cols_sql)})')
    df.to_sql(name=name, con=conn, if_exists="append", index=False)


def _write_gdf(gpkg_path: str, layer: str, gdf: gpd.GeoDataFrame):
    gdf = gdf.copy()
    if gdf.crs is None:
        gdf = gdf.set_crs("EPSG:4326")
    else:
        gdf = gdf.to_crs("EPSG:4326")
    gdf.to_file(gpkg_path, layer=layer, driver="GPKG", index=False)


@dataclass
class LayerRequirement:
    postgres_name: str
    gpkg_name: str
    spatial: bool = False
    geom_col: str = "geom"
    required: bool = True


class RefreshSDIDerivedLayersAlgorithm(QgsProcessingAlgorithm):
    SDI_SOURCE = "SDI_SOURCE"
    SDI_CONNECTION = "SDI_CONNECTION"
    SDI_GPKG = "SDI_GPKG"
    INCLUDE_QA = "INCLUDE_QA"

    def tr(self, s: str) -> str:
        return QCoreApplication.translate(self.__class__.__name__, s)

    def name(self):
        return "refresh_sdi_derived_layers"

    def displayName(self):
        return self.tr("Refresh SDI derived layers")

    def group(self):
        return self.tr("02 GIS Tools")

    def groupId(self):
        return "gis_tools"

    def createInstance(self):
        return RefreshSDIDerivedLayersAlgorithm()
    
    def shortHelpString(self):
        return """
        <b>Purpose of the Plugin</b>
        This tool refreshes SDI-derived layers in either:
        • a Postgres SDI database connection, or
        • an SDI GeoPackage exported by TfC Tools.
        It rebuilds the derived layers that are within scope of the refresher.<br>

        <b>Core refreshed layers:</b>
        • <code>transit_stop_clusters</code>
        • <code>transit_stops_auto</code>
        • <code>transit_trips_view</code>
        • <code>transit_trip_stops_sequence</code>
        • <code>transit_od_stats</code><br>

        <b>How to Use the Plugin</b>
        The plugin requires the following inputs:
        1. Input type / source selector
        2. Postgres connection OR GeoPackage path
        3. Include QA layers: if checked, the tool also attempts to refresh QA outputs<br>

        <b>Outputs</b>
        Updated tables inside the SDI or GeoPackage<br>

        <b>Notes</b>
        • In GeoPackage mode, layers/tables are overwritten in place.
        • If optional QA layers are requested but missing, the tool warns and continues.
        • This tool only refreshes layers that have an equivalent supported in the GeoPackage workflow.<br>

        <b>Documentation</b>
        For more information, refer to the User Guide.
        <a href="https://github.com/transportforcairo/tfc_tools/blob/main/tfc_tools_user_guide.pdf">TfC Tools User Guide</a>
        """

    def icon(self):
        return QIcon(_icon_path("icons", "GIS-icon.svg"))
    # Optional: some QGIS builds prefer this for SVG
    def svgIconPath(self):
        return _icon_path("icons", "GIS-icon.svg")
    
    def initAlgorithm(self, config=None):
        src_param = QgsProcessingParameterEnum(
            self.SDI_SOURCE,
            self.tr("SDI data source"),
            options=[
                "PostGIS (QGIS connection)",
                "GeoPackage (RouteLab / SDI export)",
            ],
            defaultValue=0,
        )
        src_param.setMetadata(
            {
                "widget_wrapper": {
                    "class": "tfc_tools.tfc_tools_common.sdi_source_toggle_wrapper.SDISourceToggleWidgetWrapper",
                    "enable_map": {
                        "0": [self.SDI_CONNECTION],
                        "1": [self.SDI_GPKG],
                    },
                }
            }
        )
        self.addParameter(src_param)

        self.addParameter(
            QgsProcessingParameterProviderConnection(
                self.SDI_CONNECTION,
                self.tr("SDI PostGIS connection (required when SDI data source = PostGIS)"),
                provider="postgres",
            )
        )
        self.parameterDefinition(self.SDI_CONNECTION).setFlags(
            self.parameterDefinition(self.SDI_CONNECTION).flags() | QgsProcessingParameterDefinition.FlagOptional
        )

        self.addParameter(
            QgsProcessingParameterFile(
                self.SDI_GPKG,
                self.tr("SDI GeoPackage (required when SDI data source = GeoPackage)"),
                behavior=QgsProcessingParameterFile.File,
                fileFilter="GeoPackage (*.gpkg)",
                optional=True,
            )
        )

        self.addParameter(
            QgsProcessingParameterBoolean(
                self.INCLUDE_QA,
                self.tr("Include QA layers (_stop_clusters, stops_auto)"),
                defaultValue=False,
            )
        )

    def checkParameterValues(self, parameters, context):
        source_mode = self.parameterAsEnum(parameters, self.SDI_SOURCE, context)
        connection_name = self.parameterAsString(parameters, self.SDI_CONNECTION, context)
        gpkg_path = self.parameterAsFile(parameters, self.SDI_GPKG, context)

        if source_mode == 0:
            if not connection_name:
                return (False, self.tr("Please choose an SDI PostGIS connection."))
        else:
            if not gpkg_path:
                return (False, self.tr("Please select an SDI GeoPackage file."))
            if not os.path.exists(gpkg_path):
                return (False, self.tr("The selected GeoPackage path does not exist."))
            if not gpkg_path.lower().endswith(".gpkg"):
                return (False, self.tr("The selected file is not a .gpkg GeoPackage."))

        return super().checkParameterValues(parameters, context)

    def processAlgorithm(self, parameters, context, feedback):
        source_mode = self.parameterAsEnum(parameters, self.SDI_SOURCE, context)
        include_qa = self.parameterAsBool(parameters, self.INCLUDE_QA, context)

        if source_mode == 0:
            ensure_deps(show_ui=True)
            conn_name = self.parameterAsConnectionName(parameters, self.SDI_CONNECTION, context)
            self._refresh_postgres(conn_name, include_qa, feedback)
            return {"SDI_CONNECTION": conn_name}

        gpkg_path = self.parameterAsFile(parameters, self.SDI_GPKG, context)
        self._refresh_gpkg(gpkg_path, include_qa, feedback)
        return {"SDI_GPKG": gpkg_path}

    # ---------------- Postgres path ----------------

    def _refresh_postgres(self, conn_name: str, include_qa: bool, feedback):
        engine = postgres_engine_from_qgis_connection(conn_name)
        statements = []
        if include_qa:
            statements.extend([
                "REFRESH MATERIALIZED VIEW transit._stop_clusters",
                "REFRESH MATERIALIZED VIEW transit.stops_auto",
            ])
        statements.extend([
            "REFRESH MATERIALIZED VIEW transit.trips_view",
            "REFRESH MATERIALIZED VIEW transit.trip_stops_sequence",
            "REFRESH MATERIALIZED VIEW transit.od_stats",
        ])

        with engine.begin() as conn:
            for sql in statements:
                label = sql.replace("REFRESH MATERIALIZED VIEW ", "")
                try:
                    feedback.pushInfo(f"Refreshing {label} …")
                    conn.exec_driver_sql(sql)
                except Exception as e:
                    if include_qa and label in ("transit._stop_clusters", "transit.stops_auto"):
                        feedback.reportError(f"Warning: could not refresh {label}: {e}")
                        continue
                    raise
        engine.dispose()
        feedback.pushInfo("Refresh completed.")

    # ---------------- GeoPackage path ----------------

    def _refresh_gpkg(self, gpkg_path: str, include_qa: bool, feedback):
        source = SDISource(mode="gpkg", gpkg_path=gpkg_path)
        layers = self._gpkg_layer_names(gpkg_path)

        if include_qa:
            stop_clusters = self._try_build_stop_clusters(source, layers, feedback)
            stops_auto = self._try_build_stops_auto(source, layers, stop_clusters, feedback)
            if stop_clusters is not None:
                if stop_clusters.geometry.name != "centroid":
                    stop_clusters = stop_clusters.set_geometry("centroid")
                _write_gdf(gpkg_path, "transit__stop_clusters", stop_clusters)
                feedback.pushInfo("Overwrote transit__stop_clusters")
            if stops_auto is not None:
                _write_gdf(gpkg_path, "transit_stops_auto", stops_auto)
                feedback.pushInfo("Overwrote transit_stops_auto")

        trips_view = self._try_build_trips_view(source, layers, feedback)
        if trips_view is not None:
            _write_gdf(gpkg_path, "transit_trips_view", trips_view)
            feedback.pushInfo("Overwrote transit_trips_view")
        else:
            trips_view = self._load_existing_gdf(gpkg_path, "transit_trips_view")

        trip_stops_sequence = self._build_trip_stops_sequence(source, layers, feedback)
        with sqlite3.connect(gpkg_path) as conn:
            _write_df(conn, "transit_trip_stops_sequence", trip_stops_sequence)
        feedback.pushInfo("Overwrote transit_trip_stops_sequence")

        od_stats = self._build_od_stats(source, layers, trips_view, trip_stops_sequence, feedback)
        _write_gdf(gpkg_path, "transit_od_stats", od_stats)
        feedback.pushInfo("Overwrote transit_od_stats")
        feedback.pushInfo("GeoPackage refresh completed.")

    def _gpkg_layer_names(self, gpkg_path: str) -> set[str]:
        with sqlite3.connect(gpkg_path) as conn:
            df = pd.read_sql_query(
                """
                SELECT name AS table_name
                FROM sqlite_master
                WHERE type IN ('table', 'view')
                AND name NOT LIKE 'gpkg_%'
                AND name NOT LIKE 'rtree_%'
                AND name NOT LIKE 'sqlite_%'
                """,
                conn,
            )
        return set(df["table_name"].astype(str))

    def _load_existing_gdf(self, gpkg_path: str, layer: str) -> gpd.GeoDataFrame:
        gdf = gpd.read_file(gpkg_path, layer=layer)
        if gdf.crs is None:
            gdf = gdf.set_crs("EPSG:4326")
        return gdf

    def _load_required_gdf(self, source: SDISource, layers: set[str], layer: str) -> gpd.GeoDataFrame:
        if layer not in layers:
            raise FileNotFoundError(layer)
        gdf = gpd.read_file(source.gpkg_path, layer=layer)
        if gdf.crs is None:
            gdf = gdf.set_crs("EPSG:4326")
        return gdf

    def _load_required_df(self, source: SDISource, layers: set[str], layer: str) -> pd.DataFrame:
        with sqlite3.connect(source.gpkg_path) as conn:
            exists = pd.read_sql_query(
                """
                SELECT COUNT(*) AS n
                FROM sqlite_master
                WHERE type IN ('table', 'view') AND name = ?
                """,
                conn,
                params=[layer],
            )["n"].iloc[0]
            if not exists:
                raise FileNotFoundError(layer)
        return pd.read_sql_query(f'SELECT * FROM "{layer}"', conn)

    def _warn_skip(self, feedback, label: str, reason: str):
        feedback.reportError(f"Warning: {label} was not refreshed: {reason}")

    def _try_build_stop_clusters(self, source: SDISource, layers: set[str], feedback) -> Optional[gpd.GeoDataFrame]:
        if "raw_stops" not in layers:
            self._warn_skip(feedback, "transit__stop_clusters", "raw_stops layer not found")
            return None

        raw_stops = self._load_required_gdf(source, layers, "raw_stops")
        name_col = "name" if "name" in raw_stops.columns else ("raw_name" if "raw_name" in raw_stops.columns else None)
        if name_col is None:
            raw_stops["raw_name"] = None
            name_col = "raw_name"

        pts = raw_stops[[name_col, raw_stops.geometry.name]].copy()
        pts = pts.rename(columns={name_col: "raw_name"}).set_geometry(raw_stops.geometry.name)
        pts = pts[pts.geometry.notna()].to_crs(3857).reset_index(drop=True)
        if pts.empty:
            self._warn_skip(feedback, "transit__stop_clusters", "raw_stops has no valid geometry")
            return None

        pts["cluster_id"] = self._dbscan_points(pts.geometry, eps=150.0, minpoints=3)
        valid = pts[pts["cluster_id"].notna()].copy()
        if valid.empty:
            self._warn_skip(feedback, "transit__stop_clusters", "no clusters met the DBSCAN threshold")
            return None

        rows = []
        for cid, grp in valid.groupby("cluster_id"):
            names = [str(v).strip() for v in grp["raw_name"].tolist() if pd.notna(v) and str(v).strip()]
            mode_name = Counter(names).most_common(1)[0][0] if names else "Unnamed"
            centroid = grp.geometry.unary_union.centroid
            rows.append({
                "cluster_id": int(cid),
                "n_points": int(len(grp)),
                "mode_name": mode_name,
                "centroid": centroid,
            })
        out = gpd.GeoDataFrame(rows, geometry="centroid", crs="EPSG:3857").to_crs(4326)
        return out[["cluster_id", "n_points", "mode_name", "centroid"]]

    def _try_build_stops_auto(self, source: SDISource, layers: set[str], stop_clusters: Optional[gpd.GeoDataFrame], feedback) -> Optional[gpd.GeoDataFrame]:
        if stop_clusters is None:
            self._warn_skip(feedback, "transit_stops_auto", "transit__stop_clusters could not be rebuilt")
            return None
        if "transit_trips" not in layers:
            self._warn_skip(feedback, "transit_stops_auto", "transit_trips layer not found")
            return None

        trips = self._load_required_gdf(source, layers, "transit_trips")
        trips = trips[trips.geometry.notna()].to_crs(3857).reset_index(drop=True)
        if trips.empty:
            self._warn_skip(feedback, "transit_stops_auto", "transit_trips has no valid geometry")
            return None

        sc = stop_clusters.to_crs(3857).copy()
        cell_m = 120.0
        sc["gx"] = (sc.geometry.x / cell_m).apply(math.floor)
        sc["gy"] = (sc.geometry.y / cell_m).apply(math.floor)
        sc = sc.sort_values(["gx", "gy", "n_points", "cluster_id"], ascending=[True, True, False, True])
        spaced = sc.groupby(["gx", "gy"], as_index=False).first()
        spaced = gpd.GeoDataFrame(spaced, geometry=stop_clusters.geometry.name, crs=sc.crs)

        results = []
        for _, row in spaced.iterrows():
            pt = row[stop_clusters.geometry.name]
            dists = trips.geometry.distance(pt)
            idx = dists.idxmin()
            dist_m = float(dists.loc[idx])
            if dist_m > 30.0:
                continue
            trip_geom = trips.loc[idx].geometry
            snap = trip_geom.interpolate(trip_geom.project(pt))
            start_pt = Point(list(trip_geom.coords)[0])
            end_pt = Point(list(trip_geom.coords)[-1])
            stop_type = "Terminal" if (snap.distance(start_pt) <= 75.0 or snap.distance(end_pt) <= 75.0) else "Informal"
            frac = trip_geom.project(snap, normalized=True)
            a = trip_geom.interpolate(max(0.0, frac - 0.0005), normalized=True)
            b = trip_geom.interpolate(min(1.0, frac + 0.0005), normalized=True)
            bearing = self._bearing_degrees(a, b)
            dir_bin = 1 if 45 <= int((bearing + 360) % 360) <= 225 else 0
            results.append({
                "stop_name": row["mode_name"],
                "location_type": 2 if stop_type == "Terminal" else 0,
                "stop_desc": f"{row['mode_name']} {'(Dir 0)' if dir_bin == 0 else '(Dir 1)'}",
                "stop_type": stop_type,
                "cluster_id": int(row["cluster_id"]),
                "double": int(dir_bin),
                "stop_lon": float(snap.x),
                "stop_lat": float(snap.y),
                "geometry": snap,
            })
        if not results:
            self._warn_skip(feedback, "transit_stops_auto", "no stop clusters could be snapped to transit_trips")
            return None
        out = gpd.GeoDataFrame(results, geometry="geometry", crs="EPSG:3857").to_crs(4326)
        out["stop_lon"] = out.geometry.x
        out["stop_lat"] = out.geometry.y
        return out[["stop_name", "location_type", "stop_desc", "stop_type", "cluster_id", "double", "stop_lon", "stop_lat", "geometry"]]

    def _try_build_trips_view(self, source: SDISource, layers: set[str], feedback) -> Optional[gpd.GeoDataFrame]:
        if "transit_trips" not in layers:
            self._warn_skip(feedback, "transit_trips_view", "transit_trips layer not found; keeping existing view")
            return None

        trips = self._load_required_gdf(source, layers, "transit_trips")
        agencies = self._load_required_df(source, layers, "transit_agencies") if "transit_agencies" in layers else pd.DataFrame()
        vehicles = self._load_required_df(source, layers, "transit_vehicles") if "transit_vehicles" in layers else pd.DataFrame()
        terminals = self._load_required_gdf(source, layers, "transit_terminals") if "transit_terminals" in layers else gpd.GeoDataFrame()

        df = trips.copy()
        if not agencies.empty:
            ag = agencies.copy()
            keep = [c for c in ["gid", "agency_id", "common_name", "vehicle_id"] if c in ag.columns]
            ag = ag[keep].rename(columns={"gid": "agency_gid"})
            df = df.merge(ag, left_on="agency_id", right_on="agency_gid", how="left")
        if not vehicles.empty and "vehicle_id" in df.columns:
            veh = vehicles.copy()
            keep = [c for c in ["gid", "name", "passenger_capacity"] if c in veh.columns]
            veh = veh[keep].rename(columns={"gid": "vehicle_gid", "name": "vehicle_name"})
            df = df.merge(veh, left_on="vehicle_id", right_on="vehicle_gid", how="left")
        if not terminals.empty:
            tt = terminals[[c for c in ["gid", "name"] if c in terminals.columns]].copy()
            df = df.merge(tt.rename(columns={"gid": "o_id", "name": "origin"}), on="o_id", how="left")
            df = df.merge(tt.rename(columns={"gid": "d_id", "name": "destination"}), on="d_id", how="left")

        agency_common = df["common_name"] if "common_name" in df.columns else pd.Series([None] * len(df), index=df.index)
        agency_serial = df["agency_serial"] if "agency_serial" in df.columns else pd.Series([None] * len(df), index=df.index)
        df["route_short"] = (agency_common.fillna("").astype(str).str.strip() + " " + agency_serial.fillna("").astype(str).str.strip()).str.strip()
        df["route_long"] = df.apply(
            lambda r: f"{r['destination']} - {r['origin']}" if str(r.get('origin', '')).lower() > str(r.get('destination', '')).lower() else f"{r.get('origin', '')} - {r.get('destination', '')}",
            axis=1,
        )
        df["trip_short"] = df.apply(lambda r: f"{r.get('origin', '')}-{r.get('destination', '')}", axis=1)
        df["len_km"] = df.to_crs(3857).geometry.length / 1000.0
        keep = [
            "gid", "observer_route_id", "route_type", "service_id", "route_short", "route_long",
            "observer_id", "direction_id", "o_id", "d_id", "origin", "destination", "agency_id",
            "vehicle_name", "passenger_capacity", "trip_short", "fare", "len_km", df.geometry.name,
        ]
        keep = [c for c in keep if c in df.columns]
        out = df[keep].rename(columns={"observer_route_id": "route_id", df.geometry.name: "geometry"})
        return gpd.GeoDataFrame(out, geometry="geometry", crs=trips.crs)

    def _build_trip_stops_sequence(self, source: SDISource, layers: set[str], feedback) -> pd.DataFrame:
        trips_layer = "transit_trips" if "transit_trips" in layers else "transit_trips_view"
        trips = self._load_required_gdf(source, layers, trips_layer)
        stops = self._load_required_gdf(source, layers, "transit_stops")
        agencies = self._load_required_df(source, layers, "transit_agencies") if "transit_agencies" in layers else pd.DataFrame()
        vehicles = self._load_required_df(source, layers, "transit_vehicles") if "transit_vehicles" in layers else pd.DataFrame()

        trips = trips[trips.geometry.notna()].to_crs(3857).copy()
        stops = stops[stops.geometry.notna()].to_crs(3857).copy()
        if trips.empty or stops.empty:
            raise ValueError("transit_trips/transit_trips_view and transit_stops must both contain geometry.")

        trip_cols = [c for c in ["gid", "observer_id", "agency_id"] if c in trips.columns]
        trips_meta = trips[trip_cols].copy()
        if not agencies.empty and "agency_id" in trips_meta.columns:
            ag = agencies[[c for c in ["gid", "vehicle_id"] if c in agencies.columns]].rename(columns={"gid": "agency_gid"})
            trips_meta = trips_meta.merge(ag, left_on="agency_id", right_on="agency_gid", how="left")
        if not vehicles.empty and "vehicle_id" in trips_meta.columns:
            veh = vehicles[[c for c in ["gid", "name"] if c in vehicles.columns]].rename(columns={"gid": "vehicle_gid", "name": "vehicle_name"})
            trips_meta = trips_meta.merge(veh, left_on="vehicle_id", right_on="vehicle_gid", how="left")
        meta_by_trip = trips_meta.set_index("gid").to_dict("index") if "gid" in trips_meta.columns else {}

        rows = []
        gid_counter = 1
        stop_geom_name = stops.geometry.name
        stop_name_col = "stop_name" if "stop_name" in stops.columns else ("name" if "name" in stops.columns else None)
        stop_id_col = "stop_id" if "stop_id" in stops.columns else ("gid" if "gid" in stops.columns else None)
        for _, trip in trips.iterrows():
            trip_geom = trip.geometry
            t_gid = trip.get("gid")
            observer_trip_id = trip.get("observer_id")
            vehicle_name = meta_by_trip.get(t_gid, {}).get("vehicle_name") if t_gid in meta_by_trip else trip.get("vehicle_name")
            stop_hits = []
            for _, stop in stops.iterrows():
                if trip_geom.distance(stop.geometry) > 1.0:
                    continue
                frac = trip_geom.project(stop.geometry, normalized=True)
                dist = trip_geom.project(stop.geometry, normalized=False)
                stop_hits.append({
                    "t_id": t_gid,
                    "observer_trip_id": observer_trip_id,
                    "stop_id": str(stop.get(stop_id_col)),
                    "stop_name": stop.get(stop_name_col),
                    "distance": float(dist),
                    "distance_frac": float(frac),
                    "vehicle_name": vehicle_name,
                })
            stop_hits.sort(key=lambda r: (r["distance_frac"], r["stop_id"]))
            prev_dist = None
            stop_seq = 1
            for hit in stop_hits:
                distance_from_prev = None if prev_dist is None else hit["distance"] - prev_dist
                prev_dist = hit["distance"]
                if distance_from_prev is not None and distance_from_prev < 100.0:
                    continue
                hit["gid"] = gid_counter
                hit["distance_from_prev"] = distance_from_prev
                hit["stop_sequence"] = stop_seq
                gid_counter += 1
                stop_seq += 1
                rows.append(hit)

        if not rows:
            raise ValueError("No stop-to-trip matches were found while rebuilding transit_trip_stops_sequence.")
        cols = ["gid", "t_id", "observer_trip_id", "stop_id", "stop_name", "distance", "distance_frac", "vehicle_name", "distance_from_prev", "stop_sequence"]
        return pd.DataFrame(rows)[cols]

    def _build_od_stats(self, source: SDISource, layers: set[str], trips_view: gpd.GeoDataFrame, trip_stops_sequence: pd.DataFrame, feedback) -> gpd.GeoDataFrame:
        raw_trackpoints = self._load_required_gdf(source, layers, "raw_trackpoints")
        raw_onboard = self._load_required_gdf(source, layers, "raw_onboard_instances")
        stops = self._load_required_gdf(source, layers, "transit_stops")
        intervals = self._load_required_df(source, layers, "transit_intervals")

        # ---- rebuild od_segments (same structure as export SQL) ----
        tss = trip_stops_sequence.copy().sort_values(["t_id", "stop_sequence"])
        tv = trips_view.to_crs(3857).copy()
        trip_geom_by_gid = tv.set_index("gid").geometry.to_dict() if "gid" in tv.columns else {}

        segments = []
        for t_id, grp in tss.groupby("t_id"):
            grp = grp.sort_values("stop_sequence").reset_index(drop=True)
            trip_geom = trip_geom_by_gid.get(t_id)
            if trip_geom is None:
                continue

            for i in range(len(grp) - 1):
                a = grp.iloc[i]
                b = grp.iloc[i + 1]

                if pd.isna(a["distance_frac"]) or pd.isna(b["distance_frac"]):
                    continue
                if float(a["distance_frac"]) >= 1:
                    continue

                geom = substring(
                    trip_geom,
                    float(a["distance_frac"]),
                    float(b["distance_frac"]),
                    normalized=True,
                )
                if geom.is_empty:
                    continue

                segments.append(
                    {
                        "o_id": str(a["stop_id"]),
                        "d_id": str(b["stop_id"]),
                        "vehicle_name": a.get("vehicle_name"),
                        "dist": float(b["distance"] - a["distance"]),
                        "geometry": geom,
                        "geom_key": geom.wkb_hex,
                    }
                )

        if not segments:
            raise ValueError("No OD segments could be rebuilt.")

        seg_df = pd.DataFrame(segments)
        od_segments = seg_df.groupby(
            ["o_id", "d_id", "vehicle_name", "geom_key"],
            dropna=False,
            as_index=False,
        ).agg({"dist": "mean", "geometry": "first"})
        od_segments = gpd.GeoDataFrame(od_segments, geometry="geometry", crs="EPSG:3857")

        # ---- valid/finished onboard instances ----
        raw_onboard = raw_onboard.copy()
        if "valid" in raw_onboard.columns:
            raw_onboard = raw_onboard[
                raw_onboard["valid"].astype(str).str.lower().isin(["true", "1", "t", "yes"])
            ]
        if "status" in raw_onboard.columns:
            raw_onboard = raw_onboard[
                raw_onboard["status"].astype(str).str.lower().eq("finished")
            ]

        if raw_onboard.empty:
            raise ValueError("raw_onboard_instances has no finished/valid rows for rebuilding transit_od_stats.")

        valid_ids = set(raw_onboard["id"].astype(str))

        # ---- trackpoints near stops ----
        stops3857 = (
            stops.to_crs(3857)[["gid", stops.geometry.name]]
            .copy()
            .rename(columns={stops.geometry.name: "geometry"})
            .set_geometry("geometry")
        )

        track3857 = raw_trackpoints.to_crs(3857).copy()
        if "timestamp" not in track3857.columns or "onboard_instance_id" not in track3857.columns:
            raise ValueError("raw_trackpoints must contain timestamp and onboard_instance_id columns.")

        track3857["instance_id"] = track3857["onboard_instance_id"].astype(str)
        track3857 = track3857[track3857["instance_id"].isin(valid_ids)].copy()
        track3857["time"] = pd.to_datetime(track3857["timestamp"], errors="coerce")
        track3857 = track3857[track3857["time"].notna()].copy()

        if track3857.empty:
            raise ValueError("raw_trackpoints has no valid timestamps for finished/valid onboard instances.")

        tp_near = []
        for _, tp in track3857.iterrows():
            dists = stops3857.geometry.distance(tp.geometry)
            near_idx = dists[dists <= 30.0].index.tolist()
            for idx in near_idx:
                stop_row = stops3857.loc[idx]
                tp_near.append(
                    {
                        "stop_id": str(stop_row["gid"]),
                        "gtfs_id": str(stop_row["gid"]),
                        "instance_id": tp["instance_id"],
                        "time": tp["time"],
                    }
                )

        if not tp_near:
            raise ValueError("No raw trackpoints were found within 30m of transit_stops.")

        tp_near_df = pd.DataFrame(tp_near)

        # ---- mirror export SQL: join stop hits to od_segments, carry vehicle_name from segments ----
        o_lookup = tp_near_df.merge(
            od_segments[["o_id", "vehicle_name"]].drop_duplicates(),
            left_on="stop_id",
            right_on="o_id",
            how="inner",
        )
        d_lookup = tp_near_df.merge(
            od_segments[["d_id", "vehicle_name"]].drop_duplicates(),
            left_on="stop_id",
            right_on="d_id",
            how="inner",
        )

        o_ts = (
            o_lookup.groupby(["o_id", "instance_id", "vehicle_name"], as_index=False)
            .agg({"gtfs_id": "min", "time": "max"})
            .rename(columns={"gtfs_id": "from_id", "time": "o_time"})
        )

        d_ts = (
            d_lookup.groupby(["d_id", "instance_id", "vehicle_name"], as_index=False)
            .agg({"gtfs_id": "min", "time": "max"})
            .rename(columns={"gtfs_id": "to_id", "time": "d_time"})
        )

        od_ts = o_ts.merge(d_ts, on=["instance_id", "vehicle_name"], how="inner")
        od_ts["time_diff"] = (od_ts["d_time"] - od_ts["o_time"]).dt.total_seconds()
        od_ts = od_ts[od_ts["time_diff"] > 0].copy()

        if od_ts.empty:
            raise ValueError("No positive OD travel times could be computed from raw_onboard_instances/raw_trackpoints.")

        # ---- tolerant interval parsing ----
        intervals = intervals.copy()
        if "active" in intervals.columns:
            intervals = intervals[
                intervals["active"].astype(str).str.lower().isin(["true", "1", "t", "yes"])
            ].copy()

        def _to_time(v):
            if pd.isna(v):
                return None
            s = str(v).strip()
            ts = pd.to_datetime(s, errors="coerce")
            if pd.notna(ts):
                return ts.time()
            for fmt in ("%H:%M", "%H:%M:%S"):
                try:
                    return pd.to_datetime(s, format=fmt, errors="raise").time()
                except Exception:
                    pass
            return None

        intervals["start_time"] = intervals["start_time"].apply(_to_time)
        intervals["end_time"] = intervals["end_time"].apply(_to_time)
        intervals = intervals[
            intervals["start_time"].notna() & intervals["end_time"].notna()
        ].copy()

        if feedback:
            feedback.pushInfo(f"Active intervals parsed: {len(intervals)}")
            if len(intervals):
                sample_cols = [c for c in ["gid", "start_time", "end_time", "active"] if c in intervals.columns]
                feedback.pushInfo(f"Interval sample: {intervals[sample_cols].head(5).to_dict('records')}")

        if intervals.empty:
            raise ValueError("No active intervals with parseable start/end times were found.")

        avg_rows = []
        for _, row in od_ts.iterrows():
            o_time = row["o_time"].time()
            for _, iv in intervals.iterrows():
                if iv["start_time"] <= o_time <= iv["end_time"]:
                    avg_rows.append(
                        {
                            "o_id": row["o_id"],
                            "d_id": row["d_id"],
                            "interval_id": iv["gid"],
                            "interval_start": iv["start_time"].strftime("%H:%M:%S"),
                            "vehicle_name": row["vehicle_name"],
                            "from_id": row["from_id"],
                            "to_id": row["to_id"],
                            "time_diff": row["time_diff"],
                        }
                    )

        avg_df = pd.DataFrame(avg_rows)
        if avg_df.empty:
            raise ValueError("No OD times matched active intervals.")

        avg_df = avg_df.groupby(
            ["o_id", "d_id", "interval_id", "interval_start", "vehicle_name"],
            dropna=False,
            as_index=False,
        ).agg({"from_id": "min", "to_id": "min", "time_diff": "mean"})

        avg_df["duration"] = avg_df["time_diff"].round().astype(int)

        merged = od_segments.merge(
            avg_df,
            on=["o_id", "d_id", "vehicle_name"],
            how="inner",
        )
        if merged.empty:
            raise ValueError("No OD segment durations matched rebuilt OD geometries.")

        merged = gpd.GeoDataFrame(merged, geometry="geometry", crs="EPSG:3857").to_crs(4326)
        merged.insert(0, "gid", range(1, len(merged) + 1))
        merged["speed"] = (merged["dist"] / merged["duration"]) * 3.6

        keep = [
            "gid",
            "o_id",
            "d_id",
            "interval_id",
            "interval_start",
            "vehicle_name",
            "dist",
            "duration",
            "speed",
            "geometry",
        ]
        return merged[keep]
    # ---------------- small geometry helpers ----------------

    def _dbscan_points(self, geoms: Iterable[Point], eps: float, minpoints: int) -> list[Optional[int]]:
        geoms = list(geoms)
        n = len(geoms)
        neighbors = []
        for i, g in enumerate(geoms):
            nbrs = []
            for j, h in enumerate(geoms):
                if g.distance(h) <= eps:
                    nbrs.append(j)
            neighbors.append(nbrs)

        UNVISITED = -99
        NOISE = None
        labels: list[Optional[int] | int] = [UNVISITED] * n
        cluster_id = 0
        for i in range(n):
            if labels[i] != UNVISITED:
                continue
            nbrs = neighbors[i]
            if len(nbrs) < minpoints:
                labels[i] = NOISE
                continue
            labels[i] = cluster_id
            seeds = [j for j in nbrs if j != i]
            k = 0
            while k < len(seeds):
                j = seeds[k]
                if labels[j] is NOISE:
                    labels[j] = cluster_id
                if labels[j] != UNVISITED:
                    k += 1
                    continue
                labels[j] = cluster_id
                jnbrs = neighbors[j]
                if len(jnbrs) >= minpoints:
                    for q in jnbrs:
                        if q not in seeds:
                            seeds.append(q)
                k += 1
            cluster_id += 1
        return [None if v == UNVISITED else v for v in labels]

    def _bearing_degrees(self, a: Point, b: Point) -> float:
        dx = b.x - a.x
        dy = b.y - a.y
        return math.degrees(math.atan2(dx, dy))
