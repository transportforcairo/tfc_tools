# -*- coding: utf-8 -*-

from __future__ import annotations

import os
import sqlite3
import math
from collections import Counter
from dataclasses import dataclass
from typing import Iterable, Optional

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely import STRtree
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
from ..tfc_tools_common.stop_params import (
    StopParams,
    add_stop_params_to_algorithm,
    read_stop_params_from_algorithm,
)

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

        # Stop-extraction parameters (all advanced, all with sensible defaults).
        # Must match the defaults used by rl2sdi so that re-running refresh
        # after editing stops/trips in GIS reproduces the same pipeline.
        add_stop_params_to_algorithm(self)

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

        # Read stop-extraction parameters (defaults preserved if untouched).
        stop_params = read_stop_params_from_algorithm(self, parameters, context)
        feedback.pushInfo(f"Stop-extraction parameters: {stop_params.describe()}")

        if source_mode == 0:
            ensure_deps(show_ui=True)
            conn_name = self.parameterAsConnectionName(parameters, self.SDI_CONNECTION, context)
            self._refresh_postgres(conn_name, include_qa, stop_params, feedback)
            return {"SDI_CONNECTION": conn_name}

        gpkg_path = self.parameterAsFile(parameters, self.SDI_GPKG, context)
        self._refresh_gpkg(gpkg_path, include_qa, stop_params, feedback)
        return {"SDI_GPKG": gpkg_path}

    # ---------------- Postgres path ----------------

    def _refresh_postgres(self, conn_name: str, include_qa: bool, stop_params: StopParams, feedback):
        engine = postgres_engine_from_qgis_connection(conn_name)

        if stop_params.is_default():
            # Fast path: parameters match the values that were baked into the
            # materialized views at creation time, so a cheap REFRESH suffices.
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
            return

        # Slow path: user changed one or more stop-extraction parameters.
        # A plain REFRESH would re-run the frozen SQL, so we must drop and
        # re-create the materialized views from the shared rl2sdi SQL files
        # with placeholders substituted. Views re-created in dependency order.
        feedback.pushInfo("Custom stop-extraction parameters detected — rebuilding materialized views with new values.")

        # The 3 SQL files used by rl2sdi live next to its plugin.
        rl2sdi_dir = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "..", "rl2sdi")
        )

        # Same placeholder translation used by rl2sdi/script4plugin.py.
        sql_subs = {k: v for k, v in stop_params.as_dict().items()
                    if k != "distinguish_speeds_by_vehicle"}
        if stop_params.distinguish_speeds_by_vehicle:
            sql_subs["vehicle_group_expr"] = "obs.vehicle_name"
            sql_subs["vehicle_join_condition"] = "AND od_segments.vehicle_name = avg_durations.vehicle_name"
        else:
            sql_subs["vehicle_group_expr"] = "'_pooled_'::text"
            sql_subs["vehicle_join_condition"] = ""

        def _load_sql(filename):
            with open(os.path.join(rl2sdi_dir, filename), "r") as f:
                return f.read().format(**sql_subs)

        # The materialized view refresh above also covered transit.trips_view,
        # which lives in updated_trips_view.sql. Re-run it too so the full
        # downstream chain stays consistent.
        sql_files = [
            ("transit._stop_clusters / transit.stops / transit.stops_auto",
             "create_processed_stops.sql"),
            ("transit.trips_view",            "updated_trips_view.sql"),
            ("transit.trip_stops_sequence",   "trips_stops_sequence.sql"),
            ("transit.od_stats",              "od_stats.sql"),
        ]

        with engine.begin() as conn:
            for label, fname in sql_files:
                feedback.pushInfo(f"Rebuilding {label} from {fname} …")
                sql = _load_sql(fname)
                try:
                    conn.exec_driver_sql(sql)
                except Exception as e:
                    if not include_qa and "stops_auto" in str(e).lower():
                        # QA-only failure; acceptable
                        feedback.reportError(f"Warning while running {fname}: {e}")
                        continue
                    raise

        engine.dispose()
        feedback.pushInfo("Rebuild completed.")

    # ---------------- GeoPackage path ----------------

    def _refresh_gpkg(self, gpkg_path: str, include_qa: bool, stop_params: StopParams, feedback):
        source = SDISource(mode="gpkg", gpkg_path=gpkg_path)
        layers = self._gpkg_layer_names(gpkg_path)

        if include_qa:
            stop_clusters = self._try_build_stop_clusters(source, layers, stop_params, feedback)
            stops_auto = self._try_build_stops_auto(source, layers, stop_clusters, stop_params, feedback)
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

        trip_stops_sequence = self._build_trip_stops_sequence(source, layers, stop_params, feedback)
        with sqlite3.connect(gpkg_path) as conn:
            _write_df(conn, "transit_trip_stops_sequence", trip_stops_sequence)
        feedback.pushInfo("Overwrote transit_trip_stops_sequence")

        od_stats = self._build_od_stats(source, layers, trips_view, trip_stops_sequence, stop_params, feedback)
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

    @staticmethod
    def _safe_identifier(name: str) -> str:
        """Validate a SQLite identifier for safe inline use in a query.

        SQLite does not accept bound parameters for table or column names, so
        when we select from a user-provided layer we must interpolate the name
        into the SQL text. To prevent any form of SQL injection we reject any
        name that is not a simple identifier (letters, digits, underscores).
        Callers should only pass hard-coded names from elsewhere in this
        codebase; this guard is defensive.
        """
        import re
        if not isinstance(name, str) or not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", name):
            # Intentionally not an f-string: Bandit's B608 pattern matcher
            # flags f-strings containing 'SELECT' / 'FROM' / etc. even when
            # they're inside a ValueError message, because it cannot tell the
            # difference. Keeping this as a plain constant keeps the scanner
            # quiet without concatenating the untrusted identifier into the
            # message at all.
            raise ValueError("Unsafe SQL identifier rejected")
        return name

    def _load_required_df(self, source: SDISource, layers: set[str], layer: str) -> pd.DataFrame:
        safe_layer = self._safe_identifier(layer)
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
        # safe_layer has been validated to match [A-Za-z_][A-Za-z0-9_]* above,
        # so the only possible values are plain SQL identifiers.
        return pd.read_sql_query(f'SELECT * FROM "{safe_layer}"', conn)  # nosec B608 — identifier validated above

    def _warn_skip(self, feedback, label: str, reason: str):
        feedback.reportError(f"Warning: {label} was not refreshed: {reason}")

    def _try_build_stop_clusters(self, source: SDISource, layers: set[str], stop_params: StopParams, feedback) -> Optional[gpd.GeoDataFrame]:
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

        pts["cluster_id"] = self._dbscan_points(
            pts.geometry,
            eps=stop_params.dbscan_eps_m,
            minpoints=stop_params.dbscan_minpoints,
        )
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

    def _try_build_stops_auto(self, source: SDISource, layers: set[str], stop_clusters: Optional[gpd.GeoDataFrame], stop_params: StopParams, feedback) -> Optional[gpd.GeoDataFrame]:
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
        cell_m = float(stop_params.cell_m)
        sc["gx"] = (sc.geometry.x / cell_m).apply(math.floor)
        sc["gy"] = (sc.geometry.y / cell_m).apply(math.floor)
        sc = sc.sort_values(["gx", "gy", "n_points", "cluster_id"], ascending=[True, True, False, True])
        spaced = sc.groupby(["gx", "gy"], as_index=False).first()
        spaced = gpd.GeoDataFrame(spaced, geometry=stop_clusters.geometry.name, crs=sc.crs)

        snap_max_m = float(stop_params.snap_max_m)
        terminal_m = float(stop_params.terminal_m)

        # Build an STRtree over the trip geometries and ask for each cluster
        # centroid's NEAREST trip in a single vectorized query. This replaces
        # the former `for cluster: trips.geometry.distance(pt)` loop — roughly
        # 500 × 1400 distance calls done one-at-a-time in Python — with one
        # C-accelerated call.
        trip_geoms = np.asarray(trips.geometry.values)
        tree = STRtree(trip_geoms)
        cluster_geoms = np.asarray(spaced.geometry.values)

        # query_nearest with max_distance trims non-matches upfront. On ties
        # (cluster equidistant from two trips) shapely returns all nearest
        # matches; we pick the smallest trip index to match the old
        # pandas.Series.idxmin tie-breaking semantics.
        pairs, dists = tree.query_nearest(
            cluster_geoms, max_distance=snap_max_m, return_distance=True, all_matches=True
        )
        # pairs shape: [[cluster_idx, ...], [trip_idx, ...]]
        # For each cluster index, gather all candidate (trip_idx, dist) and pick
        # the smallest trip_idx on ties.
        from collections import defaultdict
        cand = defaultdict(list)
        for k in range(len(pairs[0])):
            cidx = int(pairs[0][k])
            tidx = int(pairs[1][k])
            d = float(dists[k])
            cand[cidx].append((d, tidx))
        # Now build per-cluster best assignment: min distance, tie-break on smallest trip index.
        nearest_by_cluster = {c: min(v)[1] for c, v in cand.items()}
        nearest_dist_by_cluster = {c: min(v)[0] for c, v in cand.items()}

        # Positional arrays for spaced cluster attributes (cheap column access).
        spaced_mode_names = spaced["mode_name"].values
        spaced_cluster_ids = spaced["cluster_id"].values

        results = []
        for cidx in range(len(spaced)):
            if cidx not in nearest_by_cluster:
                continue  # no trip within snap_max_m
            tidx = nearest_by_cluster[cidx]
            dist_m = nearest_dist_by_cluster[cidx]
            if dist_m > snap_max_m:
                continue
            pt = cluster_geoms[cidx]
            trip_geom = trip_geoms[tidx]
            snap = trip_geom.interpolate(trip_geom.project(pt))
            coords = list(trip_geom.coords)
            start_pt = Point(coords[0])
            end_pt = Point(coords[-1])
            stop_type = "Terminal" if (snap.distance(start_pt) <= terminal_m or snap.distance(end_pt) <= terminal_m) else "Informal"
            frac = trip_geom.project(snap, normalized=True)
            a = trip_geom.interpolate(max(0.0, frac - 0.0005), normalized=True)
            b = trip_geom.interpolate(min(1.0, frac + 0.0005), normalized=True)
            bearing = self._bearing_degrees(a, b)
            dir_bin = 1 if 45 <= int((bearing + 360) % 360) <= 225 else 0
            results.append({
                "stop_name": spaced_mode_names[cidx],
                "location_type": 2 if stop_type == "Terminal" else 0,
                "stop_desc": f"{spaced_mode_names[cidx]} {'(Dir 0)' if dir_bin == 0 else '(Dir 1)'}",
                "stop_type": stop_type,
                "cluster_id": int(spaced_cluster_ids[cidx]),
                "double": int(dir_bin),
                # stop_lon / stop_lat deliberately NOT stored: the geometry is
                # the single source of truth. Downstream consumers (e.g. the
                # gis2gtfs stops builder) must derive coordinates from geom.
                "geometry": snap,
            })
        if not results:
            self._warn_skip(feedback, "transit_stops_auto", "no stop clusters could be snapped to transit_trips")
            return None
        out = gpd.GeoDataFrame(results, geometry="geometry", crs="EPSG:3857").to_crs(4326)
        return out[["stop_name", "location_type", "stop_desc", "stop_type", "cluster_id", "double", "geometry"]]

    def _try_build_trips_view(self, source: SDISource, layers: set[str], feedback) -> Optional[gpd.GeoDataFrame]:
        if "transit_trips" not in layers:
            # GeoPackage exports from rl2sdi only ship the view, not the base table.
            # In that case there's nothing to rebuild from; the existing
            # transit_trips_view is authoritative. Downstream steps handle this.
            if "transit_trips_view" in layers:
                feedback.pushInfo(
                    "transit_trips base layer not present (GeoPackage export). "
                    "Keeping the existing transit_trips_view as-is."
                )
            else:
                self._warn_skip(feedback, "transit_trips_view",
                                "neither transit_trips nor transit_trips_view found")
            return None

        trips = self._load_required_gdf(source, layers, "transit_trips")
        agencies = self._load_required_df(source, layers, "transit_agencies") if "transit_agencies" in layers else pd.DataFrame()
        vehicles = self._load_required_df(source, layers, "transit_vehicles") if "transit_vehicles" in layers else pd.DataFrame()
        terminals = self._load_required_gdf(source, layers, "transit_terminals") if "transit_terminals" in layers else gpd.GeoDataFrame()

        df = trips.copy()

        # If an existing transit_trips_view is in the GeoPackage (from rl2gpkg),
        # its vehicle_name / passenger_capacity / origin / destination columns
        # are the authoritative values (they were joined at export time using
        # proper PostgreSQL joins). Pull them into our working df keyed on gid
        # so the downstream agencies→vehicles chain merge only fills gaps,
        # never overwrites correct values with NULLs from broken FK links.
        if "transit_trips_view" in layers and "gid" in df.columns:
            try:
                existing_view = self._load_existing_gdf(source.gpkg_path, "transit_trips_view")
                enrich_cols = [c for c in ["vehicle_name", "passenger_capacity", "origin", "destination"]
                               if c in existing_view.columns and c not in df.columns]
                if enrich_cols and "gid" in existing_view.columns:
                    ev = existing_view[["gid"] + enrich_cols].copy()
                    # Coerce gid on both sides to str to tolerate type mismatches.
                    ev["gid"] = ev["gid"].astype(str)
                    df["gid"] = df["gid"].astype(str)
                    df = df.merge(ev, on="gid", how="left")
                    if feedback:
                        feedback.pushInfo(
                            f"Preserved columns from existing transit_trips_view: {enrich_cols}"
                        )
            except Exception as e:
                if feedback:
                    feedback.pushInfo(f"Could not read existing transit_trips_view for enrichment: {e}")

        if not agencies.empty:
            ag = agencies.copy()
            # Prefer joining on agency_id (text code) when both sides have it;
            # fall back to gid otherwise. Coerce to str either way to tolerate
            # int/text mismatches between exports.
            if "agency_id" in df.columns and "agency_id" in ag.columns:
                keep = [c for c in ["agency_id", "common_name", "vehicle_id"] if c in ag.columns]
                ag = ag[keep].copy()
                ag["agency_id"] = ag["agency_id"].astype(str)
                df["agency_id"] = df["agency_id"].astype(str)
                df = df.merge(ag, on="agency_id", how="left")
            else:
                keep = [c for c in ["gid", "agency_id", "common_name", "vehicle_id"] if c in ag.columns]
                ag = ag[keep].rename(columns={"gid": "agency_gid"})
                if "agency_gid" in ag.columns:
                    ag["agency_gid"] = ag["agency_gid"].astype(str)
                if "agency_id" in df.columns:
                    df["agency_id"] = df["agency_id"].astype(str)
                df = df.merge(ag, left_on="agency_id", right_on="agency_gid", how="left")
        if not vehicles.empty and "vehicle_id" in df.columns:
            veh = vehicles.copy()
            keep = [c for c in ["gid", "name", "passenger_capacity"] if c in veh.columns]
            veh = veh[keep].rename(columns={"gid": "vehicle_gid", "name": "vehicle_name"})
            if "vehicle_gid" in veh.columns:
                veh["vehicle_gid"] = veh["vehicle_gid"].astype(str)
            df["vehicle_id"] = df["vehicle_id"].astype(str)
            # Preserve any vehicle_name / passenger_capacity already present on
            # trips (rl2gpkg populates them directly) so the chained merge below
            # can't silently overwrite authoritative values with NULLs when
            # agency_id / vehicle_id mappings are incomplete. After the merge we
            # combine: original wins where populated, merged fills the gaps.
            _orig_vn = df["vehicle_name"] if "vehicle_name" in df.columns else None
            _orig_pc = df["passenger_capacity"] if "passenger_capacity" in df.columns else None
            if _orig_vn is not None:
                df = df.drop(columns=["vehicle_name"])
            if _orig_pc is not None:
                df = df.drop(columns=["passenger_capacity"])
            df = df.merge(veh, left_on="vehicle_id", right_on="vehicle_gid", how="left")
            if _orig_vn is not None:
                df["vehicle_name"] = _orig_vn.where(_orig_vn.notna(), df.get("vehicle_name"))
            if _orig_pc is not None:
                df["passenger_capacity"] = _orig_pc.where(_orig_pc.notna(), df.get("passenger_capacity"))
        if not terminals.empty:
            tt = terminals[[c for c in ["gid", "name"] if c in terminals.columns]].copy()
            # Same preservation pattern as vehicle_name above: if trips already
            # has origin / destination from the enriched transit_trips_view, keep
            # the originals and let the terminal merge only fill gaps.
            _orig_origin = df["origin"] if "origin" in df.columns else None
            _orig_destination = df["destination"] if "destination" in df.columns else None
            if _orig_origin is not None:
                df = df.drop(columns=["origin"])
            if _orig_destination is not None:
                df = df.drop(columns=["destination"])
            df = df.merge(tt.rename(columns={"gid": "o_id", "name": "origin"}), on="o_id", how="left")
            df = df.merge(tt.rename(columns={"gid": "d_id", "name": "destination"}), on="d_id", how="left")
            if _orig_origin is not None:
                df["origin"] = _orig_origin.where(_orig_origin.notna(), df.get("origin"))
            if _orig_destination is not None:
                df["destination"] = _orig_destination.where(_orig_destination.notna(), df.get("destination"))

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

    def _build_trip_stops_sequence(self, source: SDISource, layers: set[str], stop_params: StopParams, feedback) -> pd.DataFrame:
        # Prefer transit_trips_view: it carries vehicle_name, passenger_capacity,
        # origin/destination etc. directly from the export join, so we don't have
        # to reconstruct them from the (sometimes incomplete) agencies→vehicles
        # chain. transit_trips is a bare base table without those columns, so
        # falling back to it loses data. _refresh_gpkg runs _try_build_trips_view
        # immediately before this, so the view is guaranteed fresh.
        trips_layer = "transit_trips_view" if "transit_trips_view" in layers else "transit_trips"
        trips = self._load_required_gdf(source, layers, trips_layer)
        stops = self._load_required_gdf(source, layers, "transit_stops")
        agencies = self._load_required_df(source, layers, "transit_agencies") if "transit_agencies" in layers else pd.DataFrame()
        vehicles = self._load_required_df(source, layers, "transit_vehicles") if "transit_vehicles" in layers else pd.DataFrame()

        trips = trips[trips.geometry.notna()].to_crs(3857).copy()
        stops = stops[stops.geometry.notna()].to_crs(3857).copy()
        if trips.empty or stops.empty:
            raise ValueError("transit_trips/transit_trips_view and transit_stops must both contain geometry.")

        stop_trip_buffer_m = float(stop_params.stop_trip_buffer_m)
        min_stop_spacing_m = float(stop_params.min_stop_spacing_m)

        trip_cols = [c for c in ["gid", "observer_id", "agency_id"] if c in trips.columns]
        trips_meta = trips[trip_cols].copy()
        if not agencies.empty and "agency_id" in trips_meta.columns:
            # trips_view.agency_id carries the agency CODE (e.g. 'AMUGA'), not the
            # surrogate gid. The preferred join is agency_id (code) on both sides.
            # Fall back to gid join if the code column is missing in either table.
            if "agency_id" in agencies.columns:
                ag_cols = [c for c in ["agency_id", "vehicle_id"] if c in agencies.columns]
                ag = agencies[ag_cols].copy()
                # Coerce both sides to str to sidestep int/text merge errors.
                ag["agency_id"] = ag["agency_id"].astype(str)
                trips_meta["agency_id"] = trips_meta["agency_id"].astype(str)
                trips_meta = trips_meta.merge(ag, on="agency_id", how="left")
            else:
                ag = agencies[[c for c in ["gid", "vehicle_id"] if c in agencies.columns]].rename(
                    columns={"gid": "agency_gid"}
                )
                # Coerce join keys to a common string type to tolerate type mismatches.
                ag["agency_gid"] = ag["agency_gid"].astype(str)
                trips_meta["agency_id"] = trips_meta["agency_id"].astype(str)
                trips_meta = trips_meta.merge(ag, left_on="agency_id", right_on="agency_gid", how="left")
        if not vehicles.empty and "vehicle_id" in trips_meta.columns:
            veh = vehicles[[c for c in ["gid", "name"] if c in vehicles.columns]].rename(columns={"gid": "vehicle_gid", "name": "vehicle_name"})
            # Coerce join keys to a common string type to tolerate int/text mismatches.
            veh["vehicle_gid"] = veh["vehicle_gid"].astype(str)
            trips_meta["vehicle_id"] = trips_meta["vehicle_id"].astype(str)
            trips_meta = trips_meta.merge(veh, left_on="vehicle_id", right_on="vehicle_gid", how="left")
        meta_by_trip = trips_meta.set_index("gid").to_dict("index") if "gid" in trips_meta.columns else {}

        rows = []
        gid_counter = 1
        stop_name_col = "stop_name" if "stop_name" in stops.columns else ("name" if "name" in stops.columns else None)
        stop_id_col = "stop_id" if "stop_id" in stops.columns else ("gid" if "gid" in stops.columns else None)

        # Build an STRtree over the stop geometries once. Each bulk query returns
        # candidate (trip_index, stop_index) pairs whose stop lies within
        # stop_trip_buffer_m of the trip line — replacing the former
        # O(trips × stops) Python distance loop with a single C-accelerated call.
        stop_geoms = np.asarray(stops.geometry.values)
        tree = STRtree(stop_geoms)

        # Pre-extract stop attribute arrays for O(1) lookup in the hot loop.
        stop_ids_arr = stops[stop_id_col].values if stop_id_col else None
        stop_names_arr = stops[stop_name_col].values if stop_name_col else None

        # Iterate trips in the same order as the original implementation so that
        # gid_counter, row ordering, and observer_trip_id tiebreaks are preserved
        # bit-for-bit vs the previous pure-Python loop.
        trip_geoms_arr = np.asarray(trips.geometry.values)
        trip_gids_arr = trips["gid"].values if "gid" in trips.columns else np.full(len(trips), None)
        trip_obs_arr = trips["observer_id"].values if "observer_id" in trips.columns else np.full(len(trips), None)

        # Single bulk query: returns [[trip_idx...], [stop_idx...]] of all
        # stops within stop_trip_buffer_m of any trip.
        pairs = tree.query(trip_geoms_arr, predicate="dwithin", distance=stop_trip_buffer_m)
        # Group stop indices by trip index.
        from collections import defaultdict
        cand_by_trip: dict[int, list[int]] = defaultdict(list)
        for tidx, sidx in zip(pairs[0].tolist(), pairs[1].tolist()):
            cand_by_trip[tidx].append(sidx)

        for trip_pos in range(len(trips)):
            trip_geom = trip_geoms_arr[trip_pos]
            t_gid = trip_gids_arr[trip_pos] if trip_gids_arr is not None else None
            # Convert numpy scalars back to native Python for downstream compatibility.
            if hasattr(t_gid, "item"):
                t_gid = t_gid.item()
            observer_trip_id = trip_obs_arr[trip_pos]
            if hasattr(observer_trip_id, "item"):
                observer_trip_id = observer_trip_id.item()
            # vehicle_name resolution order:
            #   1. If trips already has a populated vehicle_name column (rl2gpkg
            #      writes one directly on transit_trips_view), trust that. It is
            #      authoritative because the export SQL joined it from vehicles
            #      at the source.
            #   2. Otherwise, use the chained agencies→vehicles lookup via
            #      meta_by_trip. This is only needed for legacy GPKGs that don't
            #      carry the direct column.
            # Previously this preferred (2) over (1), which caused vehicle_name
            # to silently become NULL whenever the agency→vehicle chain broke
            # (e.g. unmapped agency_id, missing vehicle_id), even though the
            # correct value was sitting right there in trips.vehicle_name.
            direct_vn = None
            if "vehicle_name" in trips.columns:
                direct_vn = trips.iloc[trip_pos].get("vehicle_name")
                if pd.isna(direct_vn):
                    direct_vn = None
            if direct_vn is not None:
                vehicle_name = direct_vn
            else:
                vehicle_name = meta_by_trip.get(t_gid, {}).get("vehicle_name") if t_gid in meta_by_trip else None

            cand_idxs = cand_by_trip.get(trip_pos, [])
            if not cand_idxs:
                continue

            stop_hits = []
            for sidx in cand_idxs:
                stop_geom = stop_geoms[sidx]
                frac = trip_geom.project(stop_geom, normalized=True)
                dist = trip_geom.project(stop_geom, normalized=False)
                sid_val = stop_ids_arr[sidx] if stop_ids_arr is not None else None
                sname_val = stop_names_arr[sidx] if stop_names_arr is not None else None
                stop_hits.append({
                    "t_id": t_gid,
                    "observer_trip_id": observer_trip_id,
                    "stop_id": str(sid_val),
                    "stop_name": sname_val,
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
                if distance_from_prev is not None and distance_from_prev < min_stop_spacing_m:
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

    def _build_od_stats(self, source: SDISource, layers: set[str], trips_view: gpd.GeoDataFrame, trip_stops_sequence: pd.DataFrame, stop_params: StopParams, feedback) -> gpd.GeoDataFrame:
        raw_trackpoints = self._load_required_gdf(source, layers, "raw_trackpoints")
        raw_onboard = self._load_required_gdf(source, layers, "raw_onboard_instances")
        stops = self._load_required_gdf(source, layers, "transit_stops")
        intervals = self._load_required_df(source, layers, "transit_intervals")

        trackpoint_buffer_m = float(stop_params.trackpoint_buffer_m)

        # ---- rebuild od_segments (same structure as export SQL) ----
        tss = trip_stops_sequence.copy().sort_values(["t_id", "stop_sequence"])
        tv = trips_view.to_crs(3857).copy()
        trip_geom_by_gid = tv.set_index("gid").geometry.to_dict() if "gid" in tv.columns else {}

        segments = []
        # Per-trip segment construction. The previous version used .iloc[i] /
        # .iloc[i+1] per pair which triggers Series construction on every access
        # (Pandas overhead — dominated profiling time). Here we pull every column
        # into flat numpy arrays once per trip-group and index positionally.
        for t_id, grp in tss.groupby("t_id", sort=False):
            grp = grp.sort_values("stop_sequence")
            trip_geom = trip_geom_by_gid.get(t_id)
            if trip_geom is None:
                continue
            n = len(grp)
            if n < 2:
                continue

            stop_ids = grp["stop_id"].values
            dist_fracs = grp["distance_frac"].values
            distances = grp["distance"].values
            vehicle_names = grp["vehicle_name"].values if "vehicle_name" in grp.columns else np.full(n, None)

            # Pair (i, i+1) adjacent rows
            for i in range(n - 1):
                af = dist_fracs[i]
                bf = dist_fracs[i + 1]
                # pd.isna is slow per-scalar; use numpy-style NaN test
                if af != af or bf != bf:  # NaN check
                    continue
                af = float(af)
                if af >= 1.0:
                    continue
                bf = float(bf)

                geom = substring(trip_geom, af, bf, normalized=True)
                if geom.is_empty:
                    continue

                segments.append(
                    {
                        "o_id": str(stop_ids[i]),
                        "d_id": str(stop_ids[i + 1]),
                        "vehicle_name": vehicle_names[i],
                        "dist": float(distances[i + 1] - distances[i]),
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

        # Build an STRtree over stop geometries and do one bulk dwithin query
        # for all trackpoints at once. Original code did
        # `for tp: stops.distance(tp)` which is O(T×S) in Python. For typical
        # surveys (tens of thousands of trackpoints × ~1k stops) this is the
        # single most expensive step in refresh. Replacing it with a single
        # C-accelerated STRtree query reduces it by ~100×.
        stop_geoms_3857 = np.asarray(stops3857.geometry.values)
        stop_gids_arr = stops3857["gid"].values
        tree = STRtree(stop_geoms_3857)

        tp_geoms = np.asarray(track3857.geometry.values)
        # Iterate in track3857's row order so the resulting tp_near rows
        # match the order produced by the previous nested loop (same input
        # iteration order + same per-point result ordering from STRtree).
        pairs = tree.query(tp_geoms, predicate="dwithin", distance=trackpoint_buffer_m)
        # pairs[0] = trackpoint indices, pairs[1] = stop indices
        tp_instance_ids = track3857["instance_id"].values
        tp_times = track3857["time"].values

        tp_near = []
        for tp_idx, sidx in zip(pairs[0].tolist(), pairs[1].tolist()):
            gid_val = stop_gids_arr[sidx]
            tp_near.append({
                "stop_id": str(gid_val),
                "gtfs_id": str(gid_val),
                "instance_id": tp_instance_ids[tp_idx],
                "time": tp_times[tp_idx],
            })

        if not tp_near:
            raise ValueError(f"No raw trackpoints were found within {trackpoint_buffer_m}m of transit_stops.")

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
            # The pipeline through this point was:
            #   raw_trackpoints → track3857 (filtered to valid instances)
            #                   → tp_near (trackpoints within N m of any stop)
            #                   → o_lookup / d_lookup (joined to od_segments' o_id / d_id)
            #                   → o_ts / d_ts (aggregated per instance)
            #                   → od_ts (inner join on instance_id + vehicle_name)
            #                   → filtered to time_diff > 0
            # If od_ts ends up empty, one of those stages zeroed out. Point to
            # the right one so the real cause is obvious in the feedback log.
            diag = [
                "No positive OD travel times could be computed from raw_onboard_instances/raw_trackpoints.",
                "Diagnostics:",
                f"  valid finished onboard instances (raw_onboard['id']): {len(valid_ids)}",
                f"  trackpoints after instance+timestamp filter:         {len(track3857)}",
                f"  trackpoint→stop matches (tp_near):                   {len(tp_near_df)}",
                f"  o_lookup rows (tp_near ⨝ od_segments.o_id):          {len(o_lookup)}",
                f"  d_lookup rows (tp_near ⨝ od_segments.d_id):          {len(d_lookup)}",
                f"  o_ts rows (per-instance origin aggregates):          {len(o_ts)}",
                f"  d_ts rows (per-instance destination aggregates):     {len(d_ts)}",
                f"  od_ts rows before time_diff filter:                  "
                f"{len(o_ts.merge(d_ts, on=['instance_id', 'vehicle_name'], how='inner'))}",
                "Likely causes depending on where the counts drop to zero:",
                "  - 0 after instance filter → raw_onboard['id'] doesn't match raw_trackpoints['onboard_instance_id']. "
                "Check column names in your raw_onboard_instances layer.",
                "  - 0 in tp_near → trackpoint_buffer_m is too small, or stops/trackpoints are in different CRS.",
                "  - 0 in o_lookup or d_lookup → stop_ids in od_segments don't match gids in transit_stops.",
                "  - 0 after time_diff filter → every instance visited stops in reverse time order (data issue).",
            ]
            if feedback:
                for line in diag:
                    feedback.reportError(line)
            raise ValueError(diag[0])

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

        # Previously this was a nested Python iterrows (od_ts × intervals),
        # ~782k × N_intervals iterations. Replace with direct numpy-style
        # column access — same logic, no Series construction overhead.
        ot_o_times   = od_ts["o_time"].values           # datetime64[ns]
        ot_o_ids     = od_ts["o_id"].values
        ot_d_ids     = od_ts["d_id"].values
        ot_vehicles  = od_ts["vehicle_name"].values
        ot_from_ids  = od_ts["from_id"].values
        ot_to_ids    = od_ts["to_id"].values
        ot_time_diffs = od_ts["time_diff"].values

        iv_starts   = intervals["start_time"].values    # datetime.time objects
        iv_ends     = intervals["end_time"].values
        iv_gids     = intervals["gid"].values
        # Pre-format interval_start strings once
        iv_start_strs = [iv.strftime("%H:%M:%S") for iv in iv_starts]

        avg_rows = []
        # Iterate od_ts row-by-row at the array level (no Series instantiation).
        for k in range(len(od_ts)):
            # pandas Timestamp stored as numpy datetime64; convert to datetime.time
            ts = ot_o_times[k]
            if hasattr(ts, "time"):
                o_time = ts.time()
            else:
                # numpy.datetime64 -> datetime -> time
                o_time = pd.Timestamp(ts).time()
            for j in range(len(intervals)):
                if iv_starts[j] <= o_time <= iv_ends[j]:
                    avg_rows.append(
                        {
                            "o_id": ot_o_ids[k],
                            "d_id": ot_d_ids[k],
                            "interval_id": iv_gids[j],
                            "interval_start": iv_start_strs[j],
                            "vehicle_name": ot_vehicles[k],
                            "from_id": ot_from_ids[k],
                            "to_id": ot_to_ids[k],
                            "time_diff": ot_time_diffs[k],
                        }
                    )

        avg_df = pd.DataFrame(avg_rows)
        if avg_df.empty:
            raise ValueError("No OD times matched active intervals.")

        # Vehicle-pooling setting from StopParams. When disabled (pool mode), we
        # average durations across all vehicle types serving the same OD+interval,
        # giving a larger sample for cases where per-vehicle trackpoints are sparse.
        if stop_params.distinguish_speeds_by_vehicle:
            duration_group_keys = ["o_id", "d_id", "interval_id", "interval_start", "vehicle_name"]
            merge_keys = ["o_id", "d_id", "vehicle_name"]
        else:
            duration_group_keys = ["o_id", "d_id", "interval_id", "interval_start"]
            merge_keys = ["o_id", "d_id"]

        # Tier 1 — CALCULATED. Group raw trackpoint-derived durations and count
        # the number of observations behind each row (for n_samples).
        # Use named aggregations to get flat column names directly.
        tier1 = avg_df.groupby(
            duration_group_keys,
            dropna=False,
            as_index=False,
        ).agg(
            from_id=("from_id", "min"),
            to_id=("to_id", "min"),
            duration=("time_diff", "mean"),
            n_samples=("time_diff", "count"),
        )
        tier1["duration"] = tier1["duration"].round().astype(int)
        tier1["calc_method"] = "calculated"

        # Build the full candidate set: every (o_segment × interval). In pooled
        # mode the vehicle_name column is kept from od_segments (one row per
        # vehicle x geometry), and durations from tier1 are merged without a
        # vehicle_name join key so that every vehicle category gets the pooled
        # per-(o_id,d_id,interval) duration.
        intervals_for_merge = intervals[
            [c for c in ["gid", "name", "start_time"] if c in intervals.columns]
        ].copy()
        intervals_for_merge = intervals_for_merge.rename(
            columns={"gid": "interval_id", "name": "interval_name"}
        )
        intervals_for_merge["interval_start"] = intervals_for_merge["start_time"].apply(
            lambda t: t.strftime("%H:%M:%S") if t is not None else None
        )
        # Keep interval_name alongside id/start so the downstream cartesian product
        # carries it through to the final od_stats columns. This matches the schema
        # emitted by the Postgres materialized view and by rl2gpkg, so consumers
        # (veh-pas-flow, gis2gtfs) see the same columns regardless of producer.
        keep_intervals = [c for c in ["interval_id", "interval_name", "interval_start"]
                          if c in intervals_for_merge.columns]
        intervals_for_merge = intervals_for_merge[keep_intervals]

        # Cartesian product of od_segments × intervals — each row represents
        # an (o_id, d_id, vehicle_name, geometry, interval_id) combination
        # that COULD exist. We'll then annotate duration/calc_method per row.
        candidates = od_segments.assign(_key=1).merge(
            intervals_for_merge.assign(_key=1), on="_key"
        ).drop(columns="_key")

        # Merge Tier 1 durations onto candidates.
        merged = candidates.merge(
            tier1[duration_group_keys + ["from_id", "to_id", "duration", "n_samples", "calc_method"]],
            on=merge_keys + [k for k in duration_group_keys if k in ("interval_id", "interval_start")],
            how="left",
        )

        # ---- Tier 2 — INTERPOLATED FROM TRIP-SEQUENCE NEIGHBORS ----
        # For each trip (t_id) in trip_stops_sequence, walk its stop_sequence
        # and for any (o_id→d_id) segment in this trip that has a missing
        # duration for a given interval, estimate its speed by looking at the
        # nearest calculated neighbor segments (before and/or after) on the
        # same trip in the same interval. Duration = segment_distance / speed.
        # If only one side has a calculated neighbor, use it unilaterally.
        merged = self._fill_tier2_trip_neighbors(merged, trip_stops_sequence, stop_params, feedback)

        # ---- Tier 3 — SAME SEGMENT, OTHER INTERVALS ----
        # For any (o_id, d_id, vehicle_name) still missing a duration for a
        # given interval, average its calculated durations across other
        # intervals. Speed = dist/duration, interpolated into the missing slot.
        merged = self._fill_tier3_other_intervals(merged, merge_keys)

        # ---- Tier 4 — SAME TRIP ROUTE AVERAGE SPEED ----
        # For a still-missing row, take the average speed across all
        # calculated-or-estimated segments on the same trip (same vehicle,
        # same interval). Multiply by this segment's dist to get duration.
        merged = self._fill_tier4_trip_avg_speed(merged, trip_stops_sequence, merge_keys)

        # ---- Tier 5 — SAME VEHICLE TYPE + INTERVAL CITYWIDE ----
        merged = self._fill_tier5_vehicle_citywide(merged)

        # ---- Tier 6 — CITYWIDE AVERAGE SPEED (LAST RESORT) ----
        merged = self._fill_tier6_global(merged)

        # Any row still missing duration after all tiers — either the whole
        # refresh has no trackpoint data at all, or a pathological case. Drop.
        before = len(merged)
        merged = merged[merged["duration"].notna()].copy()
        dropped = before - len(merged)
        if dropped and feedback:
            feedback.pushInfo(f"Dropped {dropped} OD segment rows that could not be filled by any tier.")

        if merged.empty:
            raise ValueError("No OD segment durations could be calculated or estimated.")

        # Fill n_samples=0 for any row that went through an estimation tier.
        merged["n_samples"] = merged["n_samples"].fillna(0).astype(int)
        merged["duration"] = merged["duration"].round().astype(int)

        merged = gpd.GeoDataFrame(merged, geometry="geometry", crs="EPSG:3857").to_crs(4326)
        merged.insert(0, "gid", range(1, len(merged) + 1))
        merged["speed"] = (merged["dist"] / merged["duration"]) * 3.6

        if feedback:
            tier_counts = merged["calc_method"].value_counts().to_dict()
            feedback.pushInfo(f"OD stats by calc_method: {tier_counts}")

        keep = [
            "gid",
            "o_id",
            "d_id",
            "interval_id",
            "interval_name",
            "interval_start",
            "vehicle_name",
            "dist",
            "duration",
            "speed",
            "calc_method",
            "n_samples",
            "geometry",
        ]
        # Defensive: if intervals didn't have a `name` column (unexpected but
        # possible with upstream schema drift), fall back to interval_id as the
        # human-readable label so downstream consumers still find the column.
        if "interval_name" not in merged.columns:
            merged["interval_name"] = merged["interval_id"].astype(str)
        return merged[keep]

    # ---- tiered estimation helpers ----

    def _fill_tier2_trip_neighbors(
        self,
        merged: pd.DataFrame,
        trip_stops_sequence: pd.DataFrame,
        stop_params: StopParams,
        feedback,
    ) -> pd.DataFrame:
        """Tier 2: fill missing durations by interpolating between the nearest
        calculated neighbor segments on the same trip (t_id) and interval.

        For each trip, walk its stops in stop_sequence order. Between each pair
        of consecutive stops i and i+1 forms a segment. If that segment is
        missing a duration for a given interval, look forward and backward along
        the sequence for the nearest segment that HAS a calculated speed for
        that same interval + vehicle_name; interpolate speed, convert to
        duration via segment dist. Uses single-sided neighbors if only one side
        has a match.
        """
        # Build a lookup: (t_id, stop_sequence_i, interval_id, vehicle_name) → speed
        # using only Tier-1 calculated rows. The segment between stop_i and
        # stop_{i+1} is identified by the stop_sequence of the ORIGIN stop.
        if trip_stops_sequence.empty:
            return merged

        tss = trip_stops_sequence.sort_values(["t_id", "stop_sequence"]).copy()
        tss["o_id"] = tss["stop_id"].astype(str)
        tss["next_stop_id"] = tss.groupby("t_id")["stop_id"].shift(-1).astype("object")
        tss = tss[tss["next_stop_id"].notna()].copy()
        tss["d_id"] = tss["next_stop_id"].astype(str)

        # For each interval, prepare speeds of calculated segments.
        tier1_rows = merged[merged["calc_method"] == "calculated"].copy()
        if tier1_rows.empty:
            return merged

        tier1_rows["speed_mps"] = tier1_rows["dist"] / tier1_rows["duration"]

        # Attach trip sequence info to merged rows so we know each row's (t_id, stop_sequence).
        # A single (o_id, d_id, vehicle_name) may appear on multiple trips — we need to process
        # each trip independently, propagating estimates to all candidate rows that share the
        # same (o_id, d_id, vehicle_name, interval).
        tier1_speed_by_key = tier1_rows.set_index(
            ["o_id", "d_id", "vehicle_name", "interval_id"]
        )["speed_mps"].to_dict()

        # Walk each trip, identifying missing segments and the surrounding calculated neighbors.
        estimates: dict[tuple, float] = {}  # (o_id, d_id, vehicle_name, interval_id) → speed_mps

        interval_ids = sorted(merged["interval_id"].dropna().unique().tolist())

        for t_id, grp in tss.groupby("t_id", sort=False):
            grp = grp.sort_values("stop_sequence").reset_index(drop=True)
            # For this trip, each row is an (o_id, d_id, vehicle_name) segment.
            o_ids = grp["o_id"].values
            d_ids = grp["d_id"].values
            vehicles = grp["vehicle_name"].values if "vehicle_name" in grp.columns else np.array([None] * len(grp))
            n = len(grp)
            if n == 0:
                continue

            # Per interval, find which segments are calculated and which are missing, and fill.
            for iv in interval_ids:
                # Build the speed array aligned to grp order.
                speeds = np.full(n, np.nan, dtype=float)
                for i in range(n):
                    k = (o_ids[i], d_ids[i], vehicles[i], iv)
                    s = tier1_speed_by_key.get(k)
                    if s is not None:
                        speeds[i] = s
                # If nothing calculated on this trip+interval, skip — tier 4 handles it.
                if not np.any(np.isfinite(speeds)):
                    continue

                # For each missing position, find nearest finite speeds on both sides.
                finite_positions = np.where(np.isfinite(speeds))[0]
                if len(finite_positions) == 0:
                    continue

                for i in range(n):
                    if np.isfinite(speeds[i]):
                        continue
                    # Nearest left and right finite positions
                    left_positions = finite_positions[finite_positions < i]
                    right_positions = finite_positions[finite_positions > i]
                    left = left_positions[-1] if len(left_positions) > 0 else None
                    right = right_positions[0] if len(right_positions) > 0 else None

                    if left is not None and right is not None:
                        # Distance-weighted average: closer neighbor has more weight.
                        # Weighting by inverse sequence distance approximates "near" on a route.
                        w_l = 1.0 / (i - left)
                        w_r = 1.0 / (right - i)
                        est_speed = (w_l * speeds[left] + w_r * speeds[right]) / (w_l + w_r)
                    elif left is not None:
                        est_speed = float(speeds[left])
                    elif right is not None:
                        est_speed = float(speeds[right])
                    else:
                        continue

                    # Record (o_id, d_id, vehicle_name, interval_id) → speed
                    # Use setdefault — first-writer-wins, so same segment appearing on
                    # multiple trips keeps the first estimate (deterministic).
                    key = (o_ids[i], d_ids[i], vehicles[i], iv)
                    estimates.setdefault(key, est_speed)

        # Apply estimates to merged rows that are still missing duration
        if not estimates:
            return merged

        mask_missing = merged["duration"].isna()
        for idx in merged.index[mask_missing]:
            key = (
                merged.at[idx, "o_id"],
                merged.at[idx, "d_id"],
                merged.at[idx, "vehicle_name"],
                merged.at[idx, "interval_id"],
            )
            speed = estimates.get(key)
            if speed is None or speed <= 0:
                continue
            dist = merged.at[idx, "dist"]
            if pd.isna(dist) or dist <= 0:
                continue
            merged.at[idx, "duration"] = float(dist) / float(speed)
            merged.at[idx, "calc_method"] = "interpolated_segment_neighbors"
            iv_match = merged.at[idx, "interval_id"]
            # Attach interval_start for newly filled rows if null
            if pd.isna(merged.at[idx, "interval_start"]):
                iv_row = merged[merged["interval_id"] == iv_match].dropna(subset=["interval_start"]).head(1)
                if len(iv_row):
                    merged.at[idx, "interval_start"] = iv_row["interval_start"].iloc[0]

        return merged

    def _fill_tier3_other_intervals(self, merged: pd.DataFrame, merge_keys: list) -> pd.DataFrame:
        """Tier 3: fill missing rows using the same segment's CALCULATED
        durations from other intervals.

        Draws from Tier-1 (calculated) rows ONLY — no cascade from prior
        estimation tiers. Same segment = same (o_id, d_id, vehicle_name).
        Computes speed = dist/duration on calculated rows, averages across
        intervals for each segment, applies that speed to the missing row's
        distance.
        """
        mask_missing = merged["duration"].isna()
        if not mask_missing.any():
            return merged

        filled_rows = merged[merged["calc_method"] == "calculated"].copy()
        if filled_rows.empty:
            return merged
        filled_rows["speed_mps"] = filled_rows["dist"] / filled_rows["duration"]
        # Average speed for the same segment across other intervals
        segment_speed = filled_rows.groupby(merge_keys, dropna=False)["speed_mps"].mean().to_dict()

        for idx in merged.index[mask_missing]:
            key = tuple(merged.at[idx, k] for k in merge_keys)
            speed = segment_speed.get(key)
            if speed is None or speed <= 0 or pd.isna(speed):
                continue
            dist = merged.at[idx, "dist"]
            if pd.isna(dist) or dist <= 0:
                continue
            merged.at[idx, "duration"] = float(dist) / float(speed)
            merged.at[idx, "calc_method"] = "interpolated_same_segment_other_interval"

        return merged

    def _fill_tier4_trip_avg_speed(
        self,
        merged: pd.DataFrame,
        trip_stops_sequence: pd.DataFrame,
        merge_keys: list,
    ) -> pd.DataFrame:
        """Tier 4: fill missing rows using the average CALCULATED speed across
        segments on the SAME trip (t_id), within the same (vehicle_name,
        interval_id).

        Draws from Tier-1 (calculated) rows ONLY — no cascade. Requires mapping
        from (o_id, d_id, vehicle_name) to the trip(s) on which they appear.
        """
        mask_missing = merged["duration"].isna()
        if not mask_missing.any() or trip_stops_sequence.empty:
            return merged

        # Build segment→trip map: (o_id, d_id, vehicle_name) → set of t_ids
        tss = trip_stops_sequence.copy()
        tss["o_id"] = tss["stop_id"].astype(str)
        tss["next_stop_id"] = tss.groupby("t_id")["stop_id"].shift(-1).astype("object")
        tss = tss[tss["next_stop_id"].notna()].copy()
        tss["d_id"] = tss["next_stop_id"].astype(str)
        seg_to_trips = tss.groupby(["o_id", "d_id", "vehicle_name"], dropna=False)["t_id"].apply(list).to_dict()

        # For each (trip, vehicle, interval), compute avg speed of available rows.
        filled_rows = merged[merged["calc_method"] == "calculated"].copy()
        if filled_rows.empty:
            return merged
        filled_rows["speed_mps"] = filled_rows["dist"] / filled_rows["duration"]

        # Attach trip ids to filled_rows via a segment→trips explode
        # Use merge_keys without interval: match on (o_id, d_id, vehicle_name)
        filled_rows["_trips"] = filled_rows.apply(
            lambda r: seg_to_trips.get((r["o_id"], r["d_id"], r["vehicle_name"]), []),
            axis=1,
        )
        exploded = filled_rows.explode("_trips").dropna(subset=["_trips"])
        exploded = exploded.rename(columns={"_trips": "t_id"})
        if exploded.empty:
            return merged
        trip_speed = exploded.groupby(["t_id", "vehicle_name", "interval_id"], dropna=False)["speed_mps"].mean().to_dict()

        for idx in merged.index[mask_missing]:
            key_seg = (merged.at[idx, "o_id"], merged.at[idx, "d_id"], merged.at[idx, "vehicle_name"])
            iv = merged.at[idx, "interval_id"]
            vehicle = merged.at[idx, "vehicle_name"]
            t_ids = seg_to_trips.get(key_seg, [])
            # Try each trip this segment belongs to, average their speeds
            speeds = []
            for t_id in t_ids:
                s = trip_speed.get((t_id, vehicle, iv))
                if s is not None and s > 0 and not pd.isna(s):
                    speeds.append(s)
            if not speeds:
                continue
            avg = sum(speeds) / len(speeds)
            dist = merged.at[idx, "dist"]
            if pd.isna(dist) or dist <= 0:
                continue
            merged.at[idx, "duration"] = float(dist) / float(avg)
            merged.at[idx, "calc_method"] = "estimated_route_avg"

        return merged

    def _fill_tier5_vehicle_citywide(self, merged: pd.DataFrame) -> pd.DataFrame:
        """Tier 5: fill missing rows with the citywide average speed for the
        same (vehicle_name, interval_id).
        """
        mask_missing = merged["duration"].isna()
        if not mask_missing.any():
            return merged

        filled_rows = merged[merged["calc_method"] == "calculated"].copy()
        if filled_rows.empty:
            return merged
        filled_rows["speed_mps"] = filled_rows["dist"] / filled_rows["duration"]
        by_vehicle = filled_rows.groupby(["vehicle_name", "interval_id"], dropna=False)["speed_mps"].mean().to_dict()

        for idx in merged.index[mask_missing]:
            key = (merged.at[idx, "vehicle_name"], merged.at[idx, "interval_id"])
            speed = by_vehicle.get(key)
            if speed is None or speed <= 0 or pd.isna(speed):
                continue
            dist = merged.at[idx, "dist"]
            if pd.isna(dist) or dist <= 0:
                continue
            merged.at[idx, "duration"] = float(dist) / float(speed)
            merged.at[idx, "calc_method"] = "estimated_vehicle_avg"

        return merged

    def _fill_tier6_global(self, merged: pd.DataFrame) -> pd.DataFrame:
        """Tier 6: fill remaining rows with a single citywide average speed
        (across all vehicles, across all intervals). Last resort.
        """
        mask_missing = merged["duration"].isna()
        if not mask_missing.any():
            return merged

        filled_rows = merged[merged["calc_method"] == "calculated"].copy()
        if filled_rows.empty:
            return merged
        filled_rows["speed_mps"] = filled_rows["dist"] / filled_rows["duration"]
        # Per-interval global average so time-of-day patterns are preserved
        by_interval = filled_rows.groupby("interval_id", dropna=False)["speed_mps"].mean().to_dict()
        global_avg = float(filled_rows["speed_mps"].mean()) if len(filled_rows) else None

        for idx in merged.index[mask_missing]:
            iv = merged.at[idx, "interval_id"]
            speed = by_interval.get(iv, global_avg)
            if speed is None or speed <= 0 or pd.isna(speed):
                continue
            dist = merged.at[idx, "dist"]
            if pd.isna(dist) or dist <= 0:
                continue
            merged.at[idx, "duration"] = float(dist) / float(speed)
            merged.at[idx, "calc_method"] = "estimated_citywide"

        return merged

    # ---------------- small geometry helpers ----------------

    def _dbscan_points(self, geoms: Iterable[Point], eps: float, minpoints: int) -> list[Optional[int]]:
        """DBSCAN over 2D Points.

        Identical result to a naive O(N²) Python implementation, but neighbor
        discovery is delegated to scipy's cKDTree (C-accelerated, O(N log N)).
        For N=1000 points this is tens of thousands of times faster. Falls
        back to the naive algorithm only if scipy is unavailable.
        """
        geoms = list(geoms)
        n = len(geoms)
        if n == 0:
            return []

        # Extract coordinates once into a (n, 2) numpy array.
        coords = np.array([(g.x, g.y) for g in geoms], dtype=float)

        # Build the neighbor index (cKDTree) and ask for all-pairs within eps.
        # query_ball_tree returns a list of lists of neighbor indices — same
        # shape/contents as the naive version but computed in C.
        try:
            from scipy.spatial import cKDTree
            tree = cKDTree(coords)
            neighbors = tree.query_ball_tree(tree, r=eps)
        except ImportError:
            # Scipy should be present in QGIS 3.40; if somehow not, keep the
            # slow-but-correct path as a safety net.
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
