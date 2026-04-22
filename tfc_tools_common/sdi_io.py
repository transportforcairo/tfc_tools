"""tfc_tools_common.sdi_io

Shared helpers for reading TfC SDI data either from:
1) A QGIS-managed Postgres connection (SDI PostGIS), or
2) A local GeoPackage that follows the TfC RouteLab SDI GeoPackage schema.

This module is designed to be imported from QGIS' Python environment.
"""

from __future__ import annotations

import os
import sqlite3
import urllib.parse
from dataclasses import dataclass
from typing import Optional

import pandas as pd
import geopandas as gpd

from qgis.core import QgsProviderRegistry, QgsDataSourceUri

# NOTE: We intentionally import SQLAlchemy lazily so the plugin can run in
# GeoPackage mode without any Postgres dependencies installed. When Postgres
# mode is used, we install the pinned deps into the plugin-local ./libs folder
# (see tfc_tools_common/deps.py) and then import SQLAlchemy.


def _get_sqlalchemy_create_engine():
    """Return sqlalchemy.create_engine, installing plugin-local deps if needed."""
    try:
        from sqlalchemy import create_engine  # type: ignore

        return create_engine
    except Exception:
        # Install pinned deps into plugin-local ./libs and retry.
        from .deps import ensure_deps, ensure_paths

        ensure_paths()
        ensure_deps(show_ui=True)

        from sqlalchemy import create_engine  # type: ignore

        return create_engine


@dataclass
class SDISource:
    """Represents the SDI source configuration."""

    mode: str  # 'postgres' | 'gpkg'
    conn_name: Optional[str] = None
    gpkg_path: Optional[str] = None

    def validate(self):
        if self.mode not in ("postgres", "gpkg"):
            raise ValueError(f"Unsupported SDI mode: {self.mode}")
        if self.mode == "postgres" and not self.conn_name:
            raise ValueError("Postgres mode requires conn_name")
        if self.mode == "gpkg":
            if not self.gpkg_path:
                raise ValueError("GeoPackage mode requires gpkg_path")
            if not os.path.exists(self.gpkg_path):
                raise FileNotFoundError(self.gpkg_path)


def _gpkg_connect(gpkg_path: str) -> sqlite3.Connection:
    return sqlite3.connect(gpkg_path)


# Allow-list of characters for SQL identifiers that must be inlined into a
# query (neither SQLite nor Postgres accept bound parameters for table or
# column names). Any identifier that does not match this pattern is rejected
# before string interpolation, preventing SQL injection even though the
# underlying f-string form of the query is, by construction, unavoidable.
#
# Accepted: optional schema prefix, then a plain identifier
#   e.g.  "transit_stops", "raw_onboard_instances", "transit.trips_view"
# Rejected: anything with quotes, semicolons, spaces, parentheses, etc.
import re as _re

_SAFE_QUALIFIED_IDENT = _re.compile(r"^[A-Za-z_][A-Za-z0-9_]*(?:\.[A-Za-z_][A-Za-z0-9_]*)?$")
_SAFE_SIMPLE_IDENT = _re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def _safe_simple_identifier(name: str) -> str:
    """Validate a simple SQL identifier (no schema). Raises on anything unsafe."""
    if not isinstance(name, str) or not _SAFE_SIMPLE_IDENT.match(name):
        raise ValueError("Unsafe SQL identifier rejected")
    return name


def _safe_qualified_identifier(name: str) -> str:
    """Validate a [schema.]identifier. Raises on anything unsafe."""
    if not isinstance(name, str) or not _SAFE_QUALIFIED_IDENT.match(name):
        raise ValueError("Unsafe SQL identifier rejected")
    return name


def postgres_engine_from_qgis_connection(conn_name: str):
    """Build a SQLAlchemy engine from a QGIS 'postgres' connection name."""
    create_engine = _get_sqlalchemy_create_engine()

    provider = QgsProviderRegistry.instance().providerMetadata("postgres")
    connection = provider.findConnection(conn_name)
    if connection is None:
        raise ValueError(f"Could not find PostgreSQL connection named '{conn_name}' in QGIS.")

    uri = QgsDataSourceUri(connection.uri())
    db_name = uri.database()
    db_user = uri.username()
    db_password = uri.password() or ""
    db_host = uri.host()
    db_port = uri.port() or "5432"

    safe_password = urllib.parse.quote_plus(db_password)
    conn_str = f"postgresql://{db_user}:{safe_password}@{db_host}:{db_port}/{db_name}"
    return create_engine(conn_str)


def gpkg_table_name(postgres_table: str) -> str:
    """Map Postgres SDI object name to GeoPackage table name."""
    mapping = {
        "transit.agencies": "transit_agencies",
        "transit.vehicles": "transit_vehicles",
        "transit.stops": "transit_stops",
        "transit.trips_view": "transit_trips_view",
        "transit.terminals": "transit_terminals",
        "transit.intervals": "transit_intervals",
        "transit.trips_intervals": "transit_trips_intervals",
        "transit.trip_stops_sequence": "transit_trip_stops_sequence",
        "transit.od_stats": "transit_od_stats",
        "transit._stop_clusters": "transit__stop_clusters",
        "transit.stops_auto": "transit_stops_auto",
        "transit.trips": "transit_trips",
        "raw.onboard_instances": "raw_onboard_instances",
        "raw.stops": "raw_stops",
        "raw.trackpoints": "raw_trackpoints",
    }
    return mapping.get(postgres_table, postgres_table.replace(".", "_"))


def read_df(source: SDISource, postgres_table: str, schema: Optional[str] = None) -> pd.DataFrame:
    """Read a non-spatial table/view as a DataFrame."""
    source.validate()

    if source.mode == "postgres":
        engine = postgres_engine_from_qgis_connection(source.conn_name)
        if schema is None and "." in postgres_table:
            schema, table = postgres_table.split(".", 1)
        else:
            table = postgres_table
        df = pd.read_sql_table(table, con=engine, schema=schema)
        engine.dispose()
        return df

    table_name = gpkg_table_name(postgres_table)
    safe_table = _safe_simple_identifier(table_name)
    with _gpkg_connect(source.gpkg_path) as conn:
        # safe_table is validated against _SAFE_SIMPLE_IDENT above; the double
        # quotes around it defend against reserved words. There is no user
        # input in this f-string — only an allow-listed identifier.
        df = pd.read_sql_query(f'SELECT * FROM "{safe_table}"', conn)  # nosec B608
    return df


def read_sql_df(source: SDISource, sql: str) -> pd.DataFrame:
    """Read query results as DataFrame (postgres only)."""
    source.validate()
    if source.mode != "postgres":
        raise ValueError("read_sql_df is only supported for postgres mode")
    engine = postgres_engine_from_qgis_connection(source.conn_name)
    df = pd.read_sql_query(sql, con=engine)
    engine.dispose()
    return df


def _ensure_geometry_name(gdf: gpd.GeoDataFrame, target_name: str = "geometry") -> gpd.GeoDataFrame:
    """Ensure GeoDataFrame's active geometry column is named target_name."""
    if gdf is None:
        return gdf
    current = gdf.geometry.name
    if current == target_name:
        return gdf
    if hasattr(gdf, "rename_geometry"):
        return gdf.rename_geometry(target_name)
    gdf = gdf.rename(columns={current: target_name})
    return gdf.set_geometry(target_name)


def read_gdf(
    source: SDISource,
    postgres_table: str,
    where: Optional[str] = None,
    geom_col: str = "geom",
) -> gpd.GeoDataFrame:
    """Read a spatial table/view as a GeoDataFrame.

    The ``postgres_table`` argument must match ``schema.name`` where both parts
    are plain identifiers (letters, digits, underscores). The ``where``
    argument, if given, must be one of a small allow-list of known-safe
    predicates; this module never accepts arbitrary user SQL.
    """
    source.validate()

    # Allow-list of predicates used by internal callers. Extending this list
    # requires a code change, not a runtime parameter — user input can never
    # reach the WHERE clause of this query.
    _ALLOWED_WHERE = {"geom IS NOT NULL"}
    if where is not None and where not in _ALLOWED_WHERE:
        raise ValueError("Unsupported WHERE predicate")

    if source.mode == "postgres":
        safe_table = _safe_qualified_identifier(postgres_table)
        engine = postgres_engine_from_qgis_connection(source.conn_name)
        conn = engine.connect()
        # safe_table is validated against _SAFE_QUALIFIED_IDENT; where is
        # validated against _ALLOWED_WHERE. No free-form input reaches the SQL.
        if where:
            sql = f"SELECT * FROM {safe_table} WHERE {where}"  # nosec B608
        else:
            sql = f"SELECT * FROM {safe_table}"  # nosec B608
        gdf = gpd.read_postgis(sql, con=conn, geom_col=geom_col)
        conn.close()
        engine.dispose()
        return _ensure_geometry_name(gdf, "geometry")

    layer = gpkg_table_name(postgres_table)
    gdf = gpd.read_file(source.gpkg_path, layer=layer)
    if gdf.crs is None:
        gdf = gdf.set_crs("EPSG:4326")
    return _ensure_geometry_name(gdf, "geometry")
