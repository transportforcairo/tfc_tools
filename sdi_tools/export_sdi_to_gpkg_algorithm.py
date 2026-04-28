# -*- coding: utf-8 -*-

"""Export SDI (PostGIS) to GeoPackage.

This is a bridge tool: it reads the already-migrated SDI schema (raw.*, transit.*)
from a PostGIS database and writes a single GeoPackage which follows the TfC
RouteLab SDI GeoPackage naming convention (transit_*, raw_*, tfc_meta).

It is intentionally read-only to the SDI database.
"""

from __future__ import annotations

import os
import sqlite3
from datetime import datetime

import pandas as pd
import geopandas as gpd

from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import (
    QgsProcessingAlgorithm,
    QgsProcessingParameterProviderConnection,
    QgsProcessingParameterFileDestination,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterNumber,
)

from ..tfc_tools_common.deps import ensure_deps
from ..tfc_tools_common.sdi_io import SDISource, read_df, read_gdf, gpkg_table_name

# for icon
from qgis.PyQt.QtGui import QIcon

def _icon_path(*parts):
    # __file__ is inside tfc_tools/rl2sdi/
    here = os.path.dirname(__file__)
    root = os.path.abspath(os.path.join(here, ".."))     # tfc_tools/
    return os.path.join(root, *parts)


def _icon_path(*parts):
    # __file__ is inside tfc_tools/rl2sdi/
    here = os.path.dirname(__file__)
    root = os.path.abspath(os.path.join(here, ".."))     # tfc_tools/
    return os.path.join(root, *parts)

def _sqlite_type_from_series(s: pd.Series) -> str:
    if pd.api.types.is_integer_dtype(s):
        return "INTEGER"
    elif pd.api.types.is_float_dtype(s):
        return "REAL"
    elif pd.api.types.is_bool_dtype(s):
        return "INTEGER"
    else:
        return "TEXT"


def _write_df(conn: sqlite3.Connection, name: str, df: pd.DataFrame):
    df = df.copy()

    # Use existing gid if present; otherwise create one
    if "gid" in df.columns:
        df["gid"] = pd.to_numeric(df["gid"], errors="coerce")
        if df["gid"].isna().any():
            raise ValueError(f"Table '{name}': existing 'gid' contains null/non-numeric values.")
        if df["gid"].duplicated().any():
            raise ValueError(f"Table '{name}': existing 'gid' is not unique.")
        df["gid"] = df["gid"].astype("int64")
    else:
        df.insert(0, "gid", range(1, len(df) + 1))

    columns_sql = []
    for col in df.columns:
        if col == "gid":
            columns_sql.append('"gid" INTEGER PRIMARY KEY')
        else:
            sql_type = _sqlite_type_from_series(df[col])
            columns_sql.append(f'"{col}" {sql_type}')

    conn.execute(f'DROP TABLE IF EXISTS "{name}"')
    conn.execute(f'CREATE TABLE "{name}" ({", ".join(columns_sql)})')

    df.to_sql(name=name, con=conn, if_exists="append", index=False)


def _write_gdf(gpkg_path: str, layer: str, gdf: gpd.GeoDataFrame):
    if gdf is None:
        return
    if gdf.crs is None:
        gdf = gdf.set_crs("EPSG:4326")
    else:
        gdf = gdf.to_crs("EPSG:4326")
    # Always write with replace semantics
    gdf.to_file(gpkg_path, layer=layer, driver="GPKG", index=False)


class ExportSDIToGeoPackageAlgorithm(QgsProcessingAlgorithm):
    SDI_CONN = "sdi_connection"
    OUT_GPKG = "out_gpkg"
    OVERWRITE = "overwrite"
    INCLUDE_RAW = "include_raw"
    INCLUDE_QA = "include_qa"
    FALLBACK_HEADWAY = "fallback_headway"

    def tr(self, s: str) -> str:
        return QCoreApplication.translate(self.__class__.__name__, s)

    def name(self):
        return "export_sdi_to_gpkg"

    def displayName(self):
        return self.tr("Export SDI (PostGIS) to GeoPackage")

    def group(self):
        return self.tr("01 RouteLab Tools")

    def groupId(self):
        return "rl_tools"
    
    def shortHelpString(self):
        return self.tr("""
    <b>Purpose of the Plugin</b>
    This tool exports a PostgreSQL/PostGIS SDI database into a GeoPackage format.
    The exported file preserves the schema structure and can be used independently of the database.<br>

    <b>How to Use the Plugin</b>
    The plugin requires the following inputs:
    1. PostgreSQL SDI connection
    2. Output folder<br>
                       
    <b>Advanced Parameters</b>
    • Headway fallback value (optional, for trips with missing headway data).<br>

    <b>Outputs</b>
    • GeoPackage replicating the SDI schema<br>

    <b>Use Cases</b>
    • Data sharing
    • Offline workflows
    • Backup of SDI datasets<br>

    <b>Documentation</b>
    For more information, refer to the User Guide.
    <a href="https://github.com/transportforcairo/tfc_tools/blob/main/tfc_tools_user_guide.pdf">TfC Tools User Guide</a>
    """)

    def createInstance(self):
        return ExportSDIToGeoPackageAlgorithm()
    
    def icon(self):
        return QIcon(_icon_path("icons", "RL-icon.svg"))
    # Optional: some QGIS builds prefer this for SVG
    def svgIconPath(self):
        return _icon_path("icons", "RL-icon.svg")

    def initAlgorithm(self, config=None):
        self.addParameter(
            QgsProcessingParameterProviderConnection(
                self.SDI_CONN,
                self.tr("SDI PostGIS connection"),
                provider="postgres",
            )
        )
        self.addParameter(
            QgsProcessingParameterFileDestination(
                self.OUT_GPKG,
                self.tr("Output GeoPackage"),
                fileFilter="GeoPackage (*.gpkg)",
            )
        )
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.OVERWRITE,
                self.tr("Overwrite output if exists"),
                defaultValue=True,
            )
        )
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.INCLUDE_RAW,
                self.tr("Include raw_* layers"),
                defaultValue=True,
            )
        )
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.INCLUDE_QA,
                self.tr("Include QA layers (clusters, stops_auto)"),
                defaultValue=True,
            )
        )

        p = QgsProcessingParameterNumber(
            self.FALLBACK_HEADWAY,
            self.tr("Headway (seconds) for empty values (optional)"),
            type=QgsProcessingParameterNumber.Integer,
            defaultValue=None,
            optional=True,
        )
        try:
            p.setFlags(p.flags() | QgsProcessingParameterNumber.FlagAdvanced)
        except Exception:
            pass
        self.addParameter(p)

    def processAlgorithm(self, parameters, context, feedback):
        ensure_deps(show_ui=True)

        conn_name = self.parameterAsConnectionName(parameters, self.SDI_CONN, context)
        out_gpkg = self.parameterAsFileOutput(parameters, self.OUT_GPKG, context)
        overwrite = self.parameterAsBool(parameters, self.OVERWRITE, context)
        include_raw = self.parameterAsBool(parameters, self.INCLUDE_RAW, context)
        include_qa = self.parameterAsBool(parameters, self.INCLUDE_QA, context)
        fallback = self.parameterAsInt(parameters, self.FALLBACK_HEADWAY, context)
        if fallback is not None and fallback <= 0:
            fallback = None

        if os.path.exists(out_gpkg):
            if overwrite:
                os.remove(out_gpkg)
            else:
                raise Exception(f"Output already exists: {out_gpkg}")

        source = SDISource(mode="postgres", conn_name=conn_name)

        feedback.pushInfo("Reading SDI tables from PostGIS…")

        # Non-spatial tables
        vehicles = read_df(source, "transit.vehicles")
        agencies = read_df(source, "transit.agencies")
        intervals = read_df(source, "transit.intervals")
        trips_intervals = read_df(source, "transit.trips_intervals")
        trip_stops_sequence = read_df(source, "transit.trip_stops_sequence")

        # Apply fallback headway if provided
        if fallback is not None and "headway_secs" in trips_intervals.columns:
            missing = trips_intervals["headway_secs"].isna()
            if missing.any():
                trips_intervals.loc[missing, "headway_secs"] = int(fallback)
                if "headway_estimation_method" in trips_intervals.columns:
                    trips_intervals.loc[missing, "headway_estimation_method"] = "from_user_default"

        # Spatial layers
        terminals = read_gdf(source, "transit.terminals", geom_col="geom")
        stops = read_gdf(source, "transit.stops", geom_col="geom")
        trips_view = read_gdf(source, "transit.trips_view", geom_col="geom")
        od_stats = read_gdf(source, "transit.od_stats", geom_col="geom")

        # Optional QA layers
        stop_clusters = None
        stops_auto = None
        if include_qa:
            try:
                stop_clusters = read_gdf(source, "transit._stop_clusters", geom_col="centroid")
            except Exception:
                stop_clusters = None
            try:
                stops_auto = read_gdf(source, "transit.stops_auto", geom_col="geom")
            except Exception:
                stops_auto = None

        # Optional raw layers
        raw_onboard = None
        raw_stops = None
        raw_trackpoints = None
        raw_identification = None
        raw_frequency = None
        if include_raw:
            try:
                raw_onboard = read_gdf(source, "raw.onboard_instances", geom_col="geometry")
            except Exception:
                # Fallback to reading as table
                raw_onboard = None
            try:
                raw_stops = read_gdf(source, "raw.stops", geom_col="geom")
            except Exception:
                raw_stops = None
            try:
                raw_trackpoints = read_gdf(source, "raw.trackpoints", geom_col="geom")
            except Exception:
                raw_trackpoints = None
            try:
                raw_identification = read_gdf(source, "raw.identification_instances", geom_col="geometry")
            except Exception:
                raw_identification = None
            try:
                raw_frequency = read_gdf(source, "raw.frequency_instances", geom_col="geometry")
            except Exception:
                raw_frequency = None

        feedback.pushInfo("Writing GeoPackage…")

        _write_gdf(out_gpkg, "transit_terminals", terminals)

        with sqlite3.connect(out_gpkg) as conn:
            meta = pd.DataFrame(
                [
                    {
                        "schema_name": "tfc_rl_sdi_gpkg",
                        "schema_version": "0.2",
                        "exported_at": datetime.now().isoformat(timespec="seconds"),
                        "generator": "TfC Tools: Export SDI (PostGIS) to GeoPackage",
                        "notes": None,
                    }
                ]
            )
            _write_df(conn, "tfc_meta", meta)
            _write_df(conn, "transit_vehicles", vehicles)
            _write_df(conn, "transit_agencies", agencies)
            _write_df(conn, "transit_intervals", intervals)
            _write_df(conn, "transit_trips_intervals", trips_intervals)
            _write_df(conn, "transit_trip_stops_sequence", trip_stops_sequence)

        # Spatial layers written after the sqlite connection closes
        
        _write_gdf(out_gpkg, "transit_stops", stops)
        _write_gdf(out_gpkg, "transit_trips_view", trips_view)
        _write_gdf(out_gpkg, "transit_od_stats", od_stats)

        if include_qa:
            if stop_clusters is not None:
                # Rename to match schema contract
                if "centroid" in stop_clusters.columns and stop_clusters.geometry.name != "centroid":
                    stop_clusters = stop_clusters.set_geometry("centroid")
                _write_gdf(out_gpkg, "transit__stop_clusters", stop_clusters)
            if stops_auto is not None:
                _write_gdf(out_gpkg, "transit_stops_auto", stops_auto)

        if include_raw:
            if raw_stops is not None:
                _write_gdf(out_gpkg, "raw_stops", raw_stops)
            if raw_trackpoints is not None:
                _write_gdf(out_gpkg, "raw_trackpoints", raw_trackpoints)
            if raw_onboard is not None:
                # raw_onboard may have many columns; keep geometry as active
                _write_gdf(out_gpkg, "raw_onboard_instances", raw_onboard)
            if raw_identification is not None:
                _write_gdf(out_gpkg, "raw_identification_instances", raw_identification)
            if raw_frequency is not None:
                _write_gdf(out_gpkg, "raw_frequency_instances", raw_frequency)

        feedback.pushInfo(f"GeoPackage written: {out_gpkg}")
        return {"OUT_GPKG": out_gpkg}
