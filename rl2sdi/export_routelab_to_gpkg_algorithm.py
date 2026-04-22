# -*- coding: utf-8 -*-

"""Export RouteLab (Observer) PostGIS to an SDI GeoPackage (read-only).

This tool connects to a RouteLab Postgres/PostGIS database (Observer schema) and
produces a GeoPackage which follows the agreed SDI GeoPackage schema
(transit_*, raw_*, tfc_meta) WITHOUT writing anything back to the RouteLab DB.

Implementation approach (Option A):
  - Use SELECT-only PostGIS queries (CTEs) to reproduce the spatial parts of RL2SDI
    (stop clustering/snapping, stop sequencing, OD stats).
  - Use pandas for headway estimation (same logic as RL2SDI), then apply an
    optional user-provided fallback headway for any remaining null values.
"""

from __future__ import annotations

import os
import contextlib
import sqlite3
from datetime import datetime

import pandas as pd
import geopandas as gpd

from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import (
    QgsProcessingAlgorithm,
    QgsProcessingParameterProviderConnection,
    QgsProcessingParameterString,
    QgsProcessingParameterFileDestination,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterNumber,
)

from ..tfc_tools_common.deps import ensure_deps
from ..tfc_tools_common.sdi_io import postgres_engine_from_qgis_connection
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

    # Prefer existing gid; create one only if missing or unusable
    use_existing_gid = False
    if "gid" in df.columns:
        gid_test = pd.to_numeric(df["gid"], errors="coerce")
        if gid_test.notna().all() and not gid_test.duplicated().any():
            df["gid"] = gid_test.astype("int64")
            use_existing_gid = True

    if not use_existing_gid:
        if "gid" in df.columns:
            df = df.drop(columns=["gid"])
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
    mode = "w" if not os.path.exists(gpkg_path) else "a"
    gdf.to_file(gpkg_path, layer=layer, driver="GPKG", index=False, mode=mode)


class ExportRouteLabToGeoPackageAlgorithm(QgsProcessingAlgorithm):
    RL_CONN = "rl_connection"
    PROJECT_ID = "project_id"
    OUT_GPKG = "out_gpkg"
    OVERWRITE = "overwrite"
    INCLUDE_RAW = "include_raw"
    INCLUDE_QA = "include_qa"
    FALLBACK_HEADWAY = "fallback_headway"

    def tr(self, s: str) -> str:
        return QCoreApplication.translate(self.__class__.__name__, s)

    def name(self):
        return "export_routelab_to_gpkg"

    def displayName(self):
        return self.tr("Export RouteLab to GeoPackage")

    def group(self):
        return self.tr("01 RouteLab Tools")

    def groupId(self):
        return "rl_tools"
    
    def shortHelpString(self):
        return self.tr("""
    <b>Purpose of the Plugin</b>
    This tool exports RouteLab field survey data into a GeoPackage format.
    The output follows TfC’s standardized schema and can be used without requiring a database connection.<br>

    <b>How to Use the Plugin</b>
    The plugin requires the following inputs:
    1. RouteLab database connection (credentials provided by TfC)
    2. Project ID — the unique identifier for your RouteLab project
    3. Headways (optional, for trips with missing headway data)
    4. Output folder<br>

    <b>Outputs</b>
    • GeoPackage containing standardized tables (raw and transit equivalents)<br>

    <b>Use Cases</b>
    • Offline analysis without PostgreSQL
    • Sharing data across teams
    • Input to GIS2GTFS, Flow, and Revenue Estimator tools<br>

    <b>Documentation</b>
    For more information, refer to the User Guide.
    <a href="https://github.com/transportforcairo/tfc_tools/blob/main/tfc_tools_user_guide.pdf">TfC Tools User Guide</a>
    """)

    def createInstance(self):
        return ExportRouteLabToGeoPackageAlgorithm()
    
    def icon(self):
        return QIcon(_icon_path("icons", "RL-icon.svg"))
    # Optional: some QGIS builds prefer this for SVG
    def svgIconPath(self):
        return _icon_path("icons", "RL-icon.svg")

    # ---- connections -------------------------------------------------
    def _dbapi_conn_from_qgis_connection(self, conn_name: str):
        """Return (engine, dbapi_conn) for a QGIS-managed Postgres connection.

        We intentionally use a DBAPI connection (engine.raw_connection()) for
        pandas/geopandas read_sql/read_postgis to avoid edge-case incompatibilities
        with SQLAlchemy engines in some QGIS/Python combos.
        """
        engine = postgres_engine_from_qgis_connection(conn_name)
        conn = engine.raw_connection()  # DBAPI connection (typically psycopg2)
        return engine, conn

    def _psycopg2_conn_from_qgis_connection(self, conn_name: str):
        """Backward-compatible helper returning a DBAPI connection.

        Some earlier iterations called this method name. We keep it to avoid
        breaking call sites, but internally we still use SQLAlchemy to build
        the connection string, then return a DBAPI connection.
        """
        engine, conn = self._dbapi_conn_from_qgis_connection(conn_name)
        # Stash engine so we can dispose it when closing the DBAPI connection.
        self._dbapi_engine = engine
        return conn

    def initAlgorithm(self, config=None):
        self.addParameter(
            QgsProcessingParameterProviderConnection(
                self.RL_CONN,
                self.tr("RouteLab PostGIS connection"),
                provider="postgres",
            )
        )
        self.addParameter(
            QgsProcessingParameterString(
                self.PROJECT_ID,
                self.tr("RouteLab Project ID"),
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
        # Older QGIS versions may not expose FlagAdvanced on this parameter type.
        with contextlib.suppress(Exception):
            p.setFlags(p.flags() | QgsProcessingParameterNumber.FlagAdvanced)
        self.addParameter(p)

        # Stop-extraction parameters (all advanced, all with sensible defaults).
        # Shared with the "Refresh SDI Derived Layers" tool so that editing
        # stops/trips in GIS and refreshing reproduces the same pipeline.
        add_stop_params_to_algorithm(self)

    # ---- SQL builders -------------------------------------------------
    # All user values are bound via pyformat placeholders (%(name)s) and
    # supplied through a params= dict at execute time. The builders below
    # expect the following binds:
    #   %(pid)s               — RouteLab project_id from the dialog
    #   %(sid)s               — resolved setting_id (UUID from settings table)
    #   %(dbscan_eps_m)s      — stop_params.dbscan_eps_m
    #   %(dbscan_minpoints)s  — stop_params.dbscan_minpoints
    #   %(snap_max_m)s        — stop_params.snap_max_m
    #   %(terminal_m)s        — stop_params.terminal_m
    #   %(cell_m)s            — stop_params.cell_m
    #   %(stop_trip_buffer_m)s — stop_params.stop_trip_buffer_m
    #   %(min_stop_spacing_m)s — stop_params.min_stop_spacing_m
    #   %(trackpoint_buffer_m)s — stop_params.trackpoint_buffer_m
    # Literal % signs in SQL (e.g. in LIKE '%test%') are escaped as %%.
    @staticmethod
    def _sql_setting_id() -> str:
        # Use ORDER BY id DESC to avoid relying on created_at.
        return (
            "SELECT id FROM settings "
            "WHERE project_id = %(pid)s AND deleted_at IS NULL "
            "ORDER BY id DESC LIMIT 1"
        )

    @staticmethod
    def _cte_vehicles() -> str:
        return """
        rl_agencies AS (
            SELECT * FROM agencies WHERE setting_id = %(sid)s AND deleted_at IS NULL
        ),
        transit_vehicles AS (
            SELECT
                ROW_NUMBER() OVER (ORDER BY name) AS gid,
                name,
                NULL::int AS passenger_capacity
            FROM (
                SELECT DISTINCT name FROM rl_agencies WHERE name IS NOT NULL
            ) x
        )
        """

    @staticmethod
    def _cte_agencies() -> str:
        return """
        rl_agencies AS (
            SELECT * FROM agencies WHERE setting_id = %(sid)s AND deleted_at IS NULL
        ),
        transit_vehicles AS (
            SELECT
                ROW_NUMBER() OVER (ORDER BY name) AS gid,
                name,
                NULL::int AS passenger_capacity
            FROM (
                SELECT DISTINCT name FROM rl_agencies WHERE name IS NOT NULL
            ) x
        ),
        transit_agencies AS (
            SELECT
                ROW_NUMBER() OVER (ORDER BY a.id) AS gid,
                a.name::text AS agency_id,
                a.name::text AS agency_name,
                NULL::text AS agency_url,
                NULL::text AS agency_timezone,
                a.name::text AS common_name,
                (a.serial IS NOT NULL) AS has_serial,
                v.gid AS vehicle_id
            FROM rl_agencies a
            LEFT JOIN transit_vehicles v ON v.name = a.name
        )
        """

    @staticmethod
    def _cte_terminals() -> str:
        return """
        transit_terminals AS (
            SELECT
                ROW_NUMBER() OVER (ORDER BY t.id) AS gid,
                t.name,
                t.id::text AS observer_id,
                t.geometry::geometry(MULTIPOLYGON,4326) AS geom
            FROM terminals t
            WHERE t.setting_id = %(sid)s
              AND t.deleted_at IS NULL
              AND t.status = 'accepted'
              AND t.geometry IS NOT NULL
              AND ST_IsValid(t.geometry)
              AND (t.name IS NULL OR lower(t.name) NOT LIKE '%%test%%')
        )
        """

    @staticmethod
    def _cte_intervals() -> str:
        return """
        transit_intervals AS (
            SELECT
                ROW_NUMBER() OVER (ORDER BY i.start, i.end, i.id) AS gid,
                to_char(i.start::time, 'HH24:MI:SS') AS start_time,
                to_char(i.end::time,   'HH24:MI:SS') AS end_time,
                i.id::text AS observer_id,
                i.name::text AS name,
                (i.deleted_at IS NULL) AS active
            FROM intervals i
            JOIN frequency_setting_intervals fsi ON i.id = fsi.interval_id
            JOIN frequency_settings fs ON fsi.frequency_setting_id = fs.id
            WHERE fs.setting_id = %(sid)s
        )
        """

    @staticmethod
    def _cte_trips() -> str:
        # trips with stable gid ordering. The f-string here only composes
        # fellow CTE builders that emit their own %%(sid)s bindings — there
        # is no user value interpolated at Python level.
        return f"""
        {ExportRouteLabToGeoPackageAlgorithm._cte_agencies()},
        {ExportRouteLabToGeoPackageAlgorithm._cte_terminals()},
        trips_raw AS (
            SELECT
                v.id::text AS observer_id,
                v.route_id::text AS observer_route_id,
                v.direction::int AS direction_id,
                v.origin_id::text AS origin_id,
                v.destination_id::text AS destination_id,
                v.agency::text AS agency,
                v.bus_number::text AS agency_serial,
                NULLIF(v.fare::text,'')::double precision AS fare,
                v.geometry::geometry(LINESTRING,4326) AS geom
            FROM v_trips_ext v
            WHERE v.setting_id = %(sid)s
              AND v.deleted_at IS NULL
              AND v.geometry IS NOT NULL
              AND ST_IsValid(v.geometry)
        ),
        transit_trips AS (
            SELECT
                ROW_NUMBER() OVER (ORDER BY tr.observer_id) AS gid,
                o.gid AS o_id,
                d.gid AS d_id,
                a.gid AS agency_id,
                tr.direction_id,
                tr.observer_id,
                tr.observer_route_id,
                COALESCE(tr.fare, (
                    SELECT CEIL(AVG(fare)) FROM trips_raw WHERE fare IS NOT NULL
                )) AS fare,
                tr.agency_serial,
                3 AS route_type,
                'Ground_Daily'::text AS service_id,
                tr.geom
            FROM trips_raw tr
            JOIN transit_terminals o ON o.observer_id = tr.origin_id
            JOIN transit_terminals d ON d.observer_id = tr.destination_id
            JOIN transit_agencies a ON a.agency_id = tr.agency
            WHERE (o.name IS NULL OR lower(o.name) NOT LIKE '%%test%%')
              AND (d.name IS NULL OR lower(d.name) NOT LIKE '%%test%%')
        )
        """  # nosec B608 — only composes sibling builders; no user values in f-string

    @staticmethod
    def _sql_trips_view() -> str:
        return f"""
        WITH
        {ExportRouteLabToGeoPackageAlgorithm._cte_trips()}
        SELECT
            t0.gid AS gid,
            t0.observer_route_id AS route_id,
            t0.route_type AS route_type,
            t0.service_id AS service_id,
            CONCAT(t1.common_name, ' ', t0.agency_serial) AS route_short,
            CASE
                WHEN LOWER(origin_terminal.name) > LOWER(dest_terminal.name)
                    THEN CONCAT(dest_terminal.name,' - ',origin_terminal.name)
                ELSE CONCAT(origin_terminal.name,' - ',dest_terminal.name)
            END AS route_long,
            t0.observer_id AS observer_id,
            t0.direction_id AS direction_id,
            t0.o_id AS o_id,
            t0.d_id AS d_id,
            t0.geom AS geom,
            origin_terminal.name AS origin,
            dest_terminal.name AS destination,
            t1.agency_id AS agency_id,
            v.name AS vehicle_name,
            v.passenger_capacity AS passenger_capacity,
            CONCAT(origin_terminal.name, '-', dest_terminal.name) AS trip_short,
            t0.fare AS fare,
            ST_Length(t0.geom::geography)/1000.0 AS len_km
        FROM transit_trips t0
        LEFT JOIN transit_agencies t1 ON t0.agency_id = t1.gid
        LEFT JOIN transit_vehicles v ON t1.vehicle_id = v.gid
        LEFT JOIN transit_terminals origin_terminal ON t0.o_id = origin_terminal.gid
        LEFT JOIN transit_terminals dest_terminal ON t0.d_id = dest_terminal.gid
        """  # nosec B608 — only composes sibling builders; no user values in f-string

    @staticmethod
    def _cte_raw_stops() -> str:
        # raw stop events from onboard_instance_stops
        return """
        r AS (
            SELECT id AS route_id FROM routes WHERE setting_id = %(sid)s AND deleted_at IS NULL
        ),
        t AS (
            SELECT id AS trip_id FROM trips JOIN r ON trips.route_id = r.route_id WHERE deleted_at IS NULL
        ),
        oi AS (
            SELECT id AS oi_id, onboard_instances.trip_id, status, valid
            FROM onboard_instances JOIN t ON onboard_instances.trip_id = t.trip_id
            WHERE deleted_at IS NULL
        ),
        raw_stops AS (
            SELECT
                ROW_NUMBER() OVER () AS gid,
                NULLIF(TRIM(tp.name), '') AS raw_name,
                tp.onboard_instance_id::text AS onboard_instance_observer_id,
                tp.created_at AS created_at,
                tp.alight As alight,
                tp.board AS board,
                NULLIF(BTRIM(tp.board_categorized  ->> 'male'),   '')::numeric AS board_male,
                NULLIF(BTRIM(tp.board_categorized  ->> 'female'), '')::numeric AS board_female,
                NULLIF(BTRIM(tp.alight_categorized ->> 'male'),   '')::numeric AS alight_male,
                NULLIF(BTRIM(tp.alight_categorized ->> 'female'), '')::numeric AS alight_female,
                tp.geometry::geometry(POINT,4326) AS geom
            FROM onboard_instance_stops tp
            JOIN oi ON tp.onboard_instance_id = oi.oi_id
            WHERE tp.deleted_at IS NULL
              AND tp.geometry IS NOT NULL
              AND ST_IsValid(tp.geometry)
        )
        """

    @staticmethod
    def _sql_stop_clusters() -> str:
        # dbscan_eps_m / dbscan_minpoints bound at execute time via params dict.
        return f"""
        WITH
        {ExportRouteLabToGeoPackageAlgorithm._cte_raw_stops()},
        pts AS (
            SELECT
                s.gid AS src_id,
                s.raw_name,
                s.geom,
                ST_Transform(s.geom, 3857) AS g3857
            FROM raw_stops s
            WHERE s.geom IS NOT NULL
        ),
        clusters AS (
            SELECT *, ST_ClusterDBSCAN(g3857, eps := %(dbscan_eps_m)s, minpoints := %(dbscan_minpoints)s) OVER () AS cluster_id
            FROM pts
        ),
        valid AS (
            SELECT * FROM clusters WHERE cluster_id IS NOT NULL
        ),
        agg AS (
            SELECT
                cluster_id,
                COUNT(*) AS n_points,
                ST_Transform(ST_Centroid(ST_Collect(g3857)), 4326) AS centroid
            FROM valid
            GROUP BY cluster_id
        ),
        name_counts AS (
            SELECT v.cluster_id, v.raw_name, COUNT(*) AS cnt
            FROM valid v
            WHERE v.raw_name IS NOT NULL
            GROUP BY v.cluster_id, v.raw_name
        ),
        name_mode AS (
            SELECT cluster_id,
                   (ARRAY_AGG(raw_name ORDER BY cnt DESC NULLS LAST, raw_name))[1] AS mode_name
            FROM name_counts
            GROUP BY cluster_id
        )
        SELECT
            a.cluster_id,
            a.n_points,
            a.centroid,
            COALESCE(n.mode_name, 'Unnamed') AS mode_name
        FROM agg a
        LEFT JOIN name_mode n USING (cluster_id)
        """  # nosec B608 — only composes sibling builders; no user values in f-string

    @staticmethod
    def _sql_stops_auto() -> str:
        # SELECT-only version of transit.stops_auto (from create_processed_stops.sql),
        # but sourced from RouteLab (raw_stops + v_trips_ext-derived transit_trips CTE).
        return f"""
        WITH
        {ExportRouteLabToGeoPackageAlgorithm._cte_trips()},
        stop_clusters AS ({ExportRouteLabToGeoPackageAlgorithm._sql_stop_clusters()}),
        const AS (
            SELECT 250::double precision AS MIN_SPACING_M,
                   %(snap_max_m)s::double precision  AS SNAP_MAX_M,
                   %(terminal_m)s::double precision  AS TERMINAL_M
        ),
        c AS (
            SELECT sc.cluster_id, sc.mode_name, sc.n_points, sc.centroid,
                   ST_Transform(sc.centroid, 3857) AS c3857
            FROM stop_clusters sc
        ),
        spaced AS (
            WITH params AS (SELECT %(cell_m)s::double precision AS cell_m),
            cells AS (
                SELECT
                    c.*,
                    (SELECT cell_m FROM params) AS cell_m,
                    FLOOR(ST_X(c.c3857) / (SELECT cell_m FROM params)) AS gx,
                    FLOOR(ST_Y(c.c3857) / (SELECT cell_m FROM params)) AS gy
                FROM c
            ),
            ranked AS (
                SELECT
                    cells.*,
                    ROW_NUMBER() OVER (PARTITION BY gx, gy ORDER BY n_points DESC, cluster_id) AS rnk
                FROM cells
            )
            SELECT cluster_id, mode_name, n_points, centroid AS geom, c3857
            FROM ranked
            WHERE rnk = 1
        ),
        nearest AS (
            SELECT s.cluster_id,
                   s.mode_name,
                   s.geom AS orig_geom,
                   t.gid AS trip_gid,
                   ST_Transform(
                     ST_ClosestPoint(ST_Transform(t.geom,3857), s.c3857), 4326
                   ) AS geom,
                   ST_Distance(s.c3857, ST_Transform(t.geom,3857)) AS dist_m,
                   t.geom AS trip_geom
            FROM spaced s
            JOIN LATERAL (
                SELECT gid, geom
                FROM transit_trips
                WHERE geom IS NOT NULL
                ORDER BY s.c3857 <-> ST_Transform(geom,3857)
                LIMIT 1
            ) t ON TRUE
        ),
        kept AS (
            SELECT * FROM nearest WHERE dist_m <= (SELECT SNAP_MAX_M FROM const)
        ),
        term AS (
            SELECT k.*,
                   CASE WHEN ST_DWithin(
                              ST_Transform(k.geom,3857),
                              ST_StartPoint(ST_Transform(k.trip_geom,3857)),
                              (SELECT TERMINAL_M FROM const)
                          )
                          OR ST_DWithin(
                              ST_Transform(k.geom,3857),
                              ST_EndPoint(ST_Transform(k.trip_geom,3857)),
                              (SELECT TERMINAL_M FROM const)
                          )
                        THEN 'Terminal' ELSE 'Informal' END AS stop_type
            FROM kept k
        ),
        bearing_calc AS (
            SELECT
                cluster_id, mode_name, stop_type, geom, trip_geom,
                degrees(
                    ST_Azimuth(
                        ST_LineInterpolatePoint(
                            ST_Transform(trip_geom,3857),
                            GREATEST(
                                ST_LineLocatePoint(ST_Transform(trip_geom,3857), ST_Transform(geom,3857)) - 0.0005,
                                0
                            )
                        ),
                        ST_LineInterpolatePoint(
                            ST_Transform(trip_geom,3857),
                            LEAST(
                                ST_LineLocatePoint(ST_Transform(trip_geom,3857), ST_Transform(geom,3857)) + 0.0005,
                                1
                            )
                        )
                    )
                ) AS bearing
            FROM term
        ),
        final AS (
            SELECT
                cluster_id,
                mode_name,
                stop_type,
                geom,
                CASE WHEN (bearing + 360)::int %% 360 BETWEEN 45 AND 225 THEN 1 ELSE 0 END AS dir_bin
            FROM bearing_calc
        )
        SELECT
            mode_name AS stop_name,
            CASE stop_type WHEN 'Terminal' THEN 2 ELSE 0 END AS location_type,
            mode_name || CASE WHEN dir_bin = 0 THEN ' (Dir 0)' ELSE ' (Dir 1)' END AS stop_desc,
            stop_type,
            cluster_id,
            dir_bin AS "double",
            -- stop_lon / stop_lat intentionally omitted: the point geometry is
            -- the single source of truth. Derive coordinates from geom at read
            -- time; storing scalar coordinates risks them desyncing from the
            -- geometry when stops are edited visually in QGIS.
            geom
        FROM final
        """  # nosec B608 — only composes sibling builders; no user values in f-string

    @staticmethod
    def _sql_stops() -> str:
        return f"""
        WITH stops_auto AS ({ExportRouteLabToGeoPackageAlgorithm._sql_stops_auto()})
        SELECT
            ROW_NUMBER() OVER (ORDER BY stop_name, cluster_id, "double") AS gid,
            (ROW_NUMBER() OVER (ORDER BY stop_name, cluster_id, "double"))::text AS stop_id,
            stop_name,
            stop_desc,
            location_type,
            "double",
            -- stop_lon / stop_lat intentionally omitted (see _sql_stops_auto).
            geom
        FROM stops_auto
        """  # nosec B608 — only composes sibling builders; no user values in f-string

    @staticmethod
    def _sql_trip_stops_sequence() -> str:
        return f"""
        WITH
        {ExportRouteLabToGeoPackageAlgorithm._cte_trips()},
        stops AS ({ExportRouteLabToGeoPackageAlgorithm._sql_stops()}),
        stop_distance_along_trip AS (
            SELECT
                ROW_NUMBER() OVER () AS gid,
                t.gid AS t_id,
                t.observer_id AS observer_trip_id,
                s.gid::text AS stop_id,
                s.stop_name AS stop_name,
                ST_LineLocatePoint(t.geom, s.geom) * ST_Length(t.geom::geography) AS distance,
                ST_LineLocatePoint(t.geom, s.geom) AS distance_frac,
                a.agency_name::text AS vehicle_name
            FROM stops s
            JOIN transit_trips t ON ST_DWithin(t.geom::geography, s.geom::geography, %(stop_trip_buffer_m)s::real)
            LEFT JOIN transit_agencies a ON t.agency_id = a.gid
        ),
        enriched_pairs AS (
            SELECT *,
                   distance - LAG(distance, 1) OVER (PARTITION BY observer_trip_id ORDER BY distance_frac) AS distance_from_prev,
                   ROW_NUMBER() OVER (PARTITION BY t_id ORDER BY distance_frac) AS stop_sequence
            FROM stop_distance_along_trip
        )
        SELECT *
        FROM enriched_pairs
        WHERE distance_from_prev >= %(min_stop_spacing_m)s OR distance_from_prev IS NULL
        """  # nosec B608 — only composes sibling builders; no user values in f-string

    @staticmethod
    def _sql_od_stats(distinguish_speeds_by_vehicle: bool) -> str:
        # SELECT-only version of od_stats.sql adapted to RouteLab base tables.
        # Mirrors the tiered estimation in
        # sdi_tools/refresh_sdi_derived_algorithm.py::_build_od_stats — every
        # (segment × interval) candidate gets a duration/calc_method by
        # cascading down 6 tiers (tiers 2-6 read Tier 1 only; no cascade
        # between estimation tiers).
        #
        # The vehicle-pooling toggle expands to several SQL snippets with
        # different alias pairs; see script4plugin.py for the mirror version
        # used by the materialized view. These fragments are SQL identifiers
        # chosen from a fixed two-branch if/else — no user input crosses a
        # trust boundary. Values (sid, trackpoint_buffer_m, etc.) ARE bound
        # via %%(name)s placeholders at execute time.
        if distinguish_speeds_by_vehicle:
            vehicle_group_expr          = "obs.vehicle_name"
            vehicle_join_condition_t1_s = "AND t1.vehicle_name = s.vehicle_name"
            vehicle_join_condition_c_t1 = "AND c.vehicle_name = t1.vehicle_name"
            tier3_group_keys            = "o_id, d_id, vehicle_name"
            tier3_join_condition_b_t3   = ("b.o_id = t3.o_id AND b.d_id = t3.d_id "
                                           "AND b.vehicle_name = t3.vehicle_name")
        else:
            vehicle_group_expr          = "'_pooled_'::text"
            vehicle_join_condition_t1_s = ""
            vehicle_join_condition_c_t1 = ""
            tier3_group_keys            = "o_id, d_id"
            tier3_join_condition_b_t3   = "b.o_id = t3.o_id AND b.d_id = t3.d_id"

        return f"""
        WITH
        {ExportRouteLabToGeoPackageAlgorithm._cte_trips()},
        {ExportRouteLabToGeoPackageAlgorithm._cte_intervals()},
        stops AS ({ExportRouteLabToGeoPackageAlgorithm._sql_stops()}),
        trip_stops_sequence AS ({ExportRouteLabToGeoPackageAlgorithm._sql_trip_stops_sequence()}),

        od_segments AS (
            SELECT DISTINCT ON (o_id, d_id, vehicle_name, geom)
                o_id,
                d_id,
                vehicle_name,
                AVG(dist) AS dist,
                geom,
                COUNT(*)
            FROM (
                SELECT * FROM (
                    SELECT
                        ts.stop_id AS o_id,
                        LEAD(ts.stop_id, 1) OVER (PARTITION BY ts.t_id ORDER BY ts.t_id, ts.stop_sequence) AS d_id,
                        ts.vehicle_name,
                        ts.distance_frac AS o_frac,
                        LEAD(ts.distance, 1) OVER (PARTITION BY ts.t_id ORDER BY ts.t_id, ts.stop_sequence) - ts.distance AS dist,
                        ST_LineSubstring(t.geom, ts.distance_frac,
                            LEAD(ts.distance_frac, 1) OVER (PARTITION BY ts.t_id ORDER BY ts.t_id, ts.stop_sequence)
                        ) AS geom
                    FROM trip_stops_sequence ts
                    JOIN transit_trips t ON ts.t_id = t.gid
                    -- distance_frac<1 must filter only the origin row (see od_stats.sql
                    -- for details); previously this was on the JOIN which also dropped
                    -- rows needed as LEAD destinations, causing last-segment-of-trip loss.
                ) x
                WHERE d_id IS NOT NULL AND o_frac < 1
            ) y
            GROUP BY o_id, d_id, vehicle_name, geom
        ),

        -- valid onboard instances for this setting
        r AS (SELECT id AS route_id FROM routes WHERE setting_id = %(sid)s AND deleted_at IS NULL),
        t AS (SELECT id AS trip_id FROM trips JOIN r ON trips.route_id = r.route_id WHERE deleted_at IS NULL),
        oi AS (
            SELECT id AS instance_id, onboard_instances.trip_id::text AS trip_id, status, valid
            FROM onboard_instances JOIN t ON onboard_instances.trip_id = t.trip_id
            WHERE deleted_at IS NULL AND valid = true AND status = 'finished'
        ),
        trackpoints AS (
            SELECT
                tp.onboard_instance_id::text AS instance_id,
                tp.timestamp AS time,
                tp.geometry::geometry(POINT,4326) AS geom
            FROM onboard_instance_track_points tp
            JOIN oi ON tp.onboard_instance_id = oi.instance_id
            WHERE tp.deleted_at IS NULL AND tp.geometry IS NOT NULL AND ST_IsValid(tp.geometry)
        ),

        trackpoints_near_stops AS (
            SELECT
                s.gid::text AS stop_id,
                s.gid::text AS gtfs_id,
                tr.instance_id,
                tr.time,
                oi.trip_id
            FROM trackpoints tr
            JOIN oi ON tr.instance_id = oi.instance_id
            JOIN stops s ON ST_DWithin(tr.geom::geography, s.geom::geography, %(trackpoint_buffer_m)s)
        ),

        o_timestamps AS (
            SELECT
                seg.o_id AS o_id,
                MIN(tns.gtfs_id) AS o_gtfs_id,
                tns.instance_id AS o_instance_id,
                MAX(tns.time) AS o_time,
                seg.vehicle_name AS vehicle_name
            FROM trackpoints_near_stops tns
            JOIN od_segments seg ON tns.stop_id::text = seg.o_id
            GROUP BY seg.o_id, tns.instance_id, seg.vehicle_name
        ),
        d_timestamps AS (
            SELECT
                seg.d_id AS d_id,
                MIN(tns.gtfs_id) AS d_gtfs_id,
                tns.instance_id AS d_instance_id,
                MAX(tns.time) AS d_time,
                seg.vehicle_name AS vehicle_name
            FROM trackpoints_near_stops tns
            JOIN od_segments seg ON tns.stop_id::text = seg.d_id
            GROUP BY seg.d_id, tns.instance_id, seg.vehicle_name
        ),
        od_timestamps AS (
            SELECT
                o_id,
                d_id,
                o_gtfs_id AS from_id,
                d_gtfs_id AS to_id,
                o_time,
                d_time,
                o_instance_id AS instance_id,
                o.vehicle_name
            FROM o_timestamps o
            JOIN d_timestamps d
              ON o.o_instance_id = d.d_instance_id AND o.vehicle_name = d.vehicle_name
        ),
        -- Per-instance OD observations with interval assignment (same in both modes).
        od_observations AS (
            SELECT
                tmp.o_id, tmp.d_id, tmp.vehicle_name, tmp.from_id, tmp.to_id,
                tmp.time_diff, i.gid AS interval_id, i.name AS interval_name,
                i.start_time AS interval_start
            FROM (
                SELECT *, EXTRACT(EPOCH FROM (d_time - o_time)) AS time_diff
                FROM od_timestamps
            ) tmp
            JOIN transit_intervals i ON (tmp.o_time::time BETWEEN i.start_time::time AND i.end_time::time) AND i.active
            WHERE tmp.time_diff > 0
        ),

        -- ===== Tier 1 — CALCULATED durations from trackpoint samples =====
        tier1_calculated AS (
            SELECT
                obs.o_id,
                obs.d_id,
                obs.interval_id,
                MIN(obs.interval_name)  AS interval_name,
                MIN(obs.vehicle_name)   AS vehicle_name,
                MIN(obs.from_id)        AS from_id,
                MIN(obs.to_id)          AS to_id,
                MIN(obs.interval_start) AS interval_start,
                AVG(obs.time_diff)::int AS duration,
                COUNT(*)                AS n_samples,
                'calculated'::text      AS calc_method
            FROM od_observations obs
            GROUP BY obs.o_id, obs.d_id, obs.interval_id, {vehicle_group_expr}
        ),

        -- ===== Phase 2: candidate set × interval and tier 2-6 estimation =====
        intervals_active AS (
            SELECT gid AS interval_id, name AS interval_name, start_time AS interval_start
            FROM transit_intervals WHERE active = true
        ),
        candidates AS (
            SELECT seg.o_id, seg.d_id, seg.vehicle_name, seg.dist, seg.geom,
                   iv.interval_id, iv.interval_name, iv.interval_start
            FROM od_segments seg CROSS JOIN intervals_active iv
        ),
        -- Pre-aggregate od_segments to one row per (o_id,d_id,vehicle_name) so the
        -- downstream tier1_with_speed join doesn't multiply rows when the same OD
        -- segment was produced with different geometric slices across trips.
        od_segments_by_key AS (
            SELECT o_id, d_id, vehicle_name, MAX(dist) AS dist
            FROM od_segments
            GROUP BY o_id, d_id, vehicle_name
        ),
        tier1_with_speed AS (
            SELECT t1.o_id, t1.d_id, t1.interval_id,
                   t1.vehicle_name AS t1_vehicle_name,
                   s.vehicle_name  AS vehicle_name,
                   s.dist          AS seg_dist,
                   t1.duration, t1.n_samples,
                   (s.dist::float / NULLIF(t1.duration, 0))::float AS speed_mps
            FROM tier1_calculated t1
            JOIN od_segments_by_key s
              ON t1.o_id = s.o_id AND t1.d_id = s.d_id
             {vehicle_join_condition_t1_s}
        ),
        base AS (
            SELECT c.o_id, c.d_id, c.vehicle_name, c.dist, c.geom,
                   c.interval_id, c.interval_name, c.interval_start,
                   t1.from_id, t1.to_id,
                   t1.duration    AS t1_duration,
                   t1.n_samples   AS t1_n_samples,
                   t1.calc_method AS t1_calc_method,
                   t1w.speed_mps  AS t1_speed_mps
            FROM candidates c
            LEFT JOIN tier1_calculated t1
              ON c.o_id = t1.o_id AND c.d_id = t1.d_id AND c.interval_id = t1.interval_id
              {vehicle_join_condition_c_t1}
            LEFT JOIN tier1_with_speed t1w
              ON c.o_id = t1w.o_id AND c.d_id = t1w.d_id AND c.interval_id = t1w.interval_id
              AND c.vehicle_name = t1w.vehicle_name
        ),
        tss_segments AS (
            SELECT t_id,
                   stop_id::text AS o_id,
                   LEAD(stop_id, 1) OVER (PARTITION BY t_id ORDER BY stop_sequence)::text AS d_id,
                   vehicle_name,
                   stop_sequence
            FROM trip_stops_sequence
        ),
        -- Expand base: one row per (candidate × trip containing its segment).
        -- Mirrors Python _fill_tier2_trip_neighbors which walks EVERY trip
        -- independently and uses setdefault() so the first trip's estimate wins.
        -- A single-trip representative (earlier version: seg_trip_one) was too
        -- restrictive — it caused candidates and their tier-1 neighbors to land
        -- in different window partitions when they were each "represented" by
        -- different trips.
        with_trip_ctx_all AS (
            SELECT b.*, st.t_id, st.stop_sequence
            FROM base b
            LEFT JOIN tss_segments st
              ON b.o_id = st.o_id AND b.d_id = st.d_id AND b.vehicle_name = st.vehicle_name
             AND st.d_id IS NOT NULL
        ),
        -- Tier 2: inverse-seq-distance weighted neighbors on same trip+interval+vehicle
        with_neighbors AS (
            SELECT w.*,
                MAX(stop_sequence) FILTER (WHERE t1_speed_mps IS NOT NULL)
                  OVER (PARTITION BY t_id, interval_id, vehicle_name
                        ORDER BY stop_sequence
                        ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW) AS prev_seq,
                (array_agg(t1_speed_mps) FILTER (WHERE t1_speed_mps IS NOT NULL)
                  OVER (PARTITION BY t_id, interval_id, vehicle_name
                        ORDER BY stop_sequence
                        ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW))[
                    array_length(
                        array_agg(t1_speed_mps) FILTER (WHERE t1_speed_mps IS NOT NULL)
                          OVER (PARTITION BY t_id, interval_id, vehicle_name
                                ORDER BY stop_sequence
                                ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW), 1)
                  ] AS prev_speed,
                MIN(stop_sequence) FILTER (WHERE t1_speed_mps IS NOT NULL)
                  OVER (PARTITION BY t_id, interval_id, vehicle_name
                        ORDER BY stop_sequence
                        ROWS BETWEEN CURRENT ROW AND UNBOUNDED FOLLOWING) AS next_seq,
                (array_agg(t1_speed_mps) FILTER (WHERE t1_speed_mps IS NOT NULL)
                  OVER (PARTITION BY t_id, interval_id, vehicle_name
                        ORDER BY stop_sequence
                        ROWS BETWEEN CURRENT ROW AND UNBOUNDED FOLLOWING))[1] AS next_speed
            FROM with_trip_ctx_all w
        ),
        t2_per_trip AS (
            SELECT o_id, d_id, vehicle_name, interval_id, t_id,
                   CASE
                       WHEN prev_seq IS NOT NULL AND next_seq IS NOT NULL THEN
                           ((1.0/(stop_sequence - prev_seq)) * prev_speed
                          + (1.0/(next_seq - stop_sequence)) * next_speed)
                         / ((1.0/(stop_sequence - prev_seq))
                          + (1.0/(next_seq - stop_sequence)))
                       WHEN prev_seq IS NOT NULL THEN prev_speed
                       WHEN next_seq IS NOT NULL THEN next_speed
                       ELSE NULL
                   END AS t2_speed
            FROM with_neighbors
            WHERE t1_duration IS NULL
              AND t_id IS NOT NULL
        ),
        t2_speed AS (
            SELECT DISTINCT ON (o_id, d_id, vehicle_name, interval_id)
                   o_id, d_id, vehicle_name, interval_id, t2_speed
            FROM t2_per_trip
            WHERE t2_speed IS NOT NULL
            ORDER BY o_id, d_id, vehicle_name, interval_id, t_id
        ),
        -- Tier 3: same segment keys, avg across other intervals
        tier3_speed AS (
            SELECT {tier3_group_keys}, AVG(speed_mps) AS t3_speed
            FROM tier1_with_speed
            GROUP BY {tier3_group_keys}
        ),
        -- Tier 4: per-(trip,vehicle,interval) avg, then avg across trips the segment uses
        t1_trip_link AS (
            SELECT t1w.o_id, t1w.d_id, t1w.vehicle_name, t1w.interval_id,
                   t1w.speed_mps, st.t_id
            FROM tier1_with_speed t1w
            JOIN tss_segments st
              ON t1w.o_id = st.o_id AND t1w.d_id = st.d_id AND t1w.vehicle_name = st.vehicle_name
            WHERE st.d_id IS NOT NULL
        ),
        tier4_trip_speed AS (
            SELECT t_id, vehicle_name, interval_id, AVG(speed_mps) AS trip_speed
            FROM t1_trip_link
            GROUP BY t_id, vehicle_name, interval_id
        ),
        seg_all_trips AS (
            SELECT DISTINCT o_id, d_id, vehicle_name, t_id
            FROM tss_segments
            WHERE d_id IS NOT NULL
        ),
        tier4_speed AS (
            SELECT sat.o_id, sat.d_id, sat.vehicle_name, tts.interval_id,
                   AVG(tts.trip_speed) AS t4_speed
            FROM seg_all_trips sat
            JOIN tier4_trip_speed tts
              ON sat.t_id = tts.t_id AND sat.vehicle_name = tts.vehicle_name
            GROUP BY sat.o_id, sat.d_id, sat.vehicle_name, tts.interval_id
        ),
        -- Tier 5: per-(vehicle, interval) citywide
        tier5_speed AS (
            SELECT vehicle_name, interval_id, AVG(speed_mps) AS t5_speed
            FROM tier1_with_speed
            GROUP BY vehicle_name, interval_id
        ),
        -- Tier 6: per-interval + overall fallback
        tier6_per_interval AS (
            SELECT interval_id, AVG(speed_mps) AS t6_speed_iv
            FROM tier1_with_speed GROUP BY interval_id
        ),
        tier6_global AS (
            SELECT AVG(speed_mps) AS t6_speed_all FROM tier1_with_speed
        ),
        assembled AS (
            SELECT
                b.o_id, b.d_id, b.vehicle_name,
                b.interval_id, b.interval_name, b.interval_start,
                b.dist, b.geom, b.from_id, b.to_id,
                b.t1_duration, b.t1_n_samples, b.t1_calc_method,
                t2.t2_speed,
                t3.t3_speed, t4.t4_speed, t5.t5_speed, t6iv.t6_speed_iv,
                (SELECT t6_speed_all FROM tier6_global) AS t6_speed_all
            FROM base b
            LEFT JOIN t2_speed t2
              ON b.o_id = t2.o_id AND b.d_id = t2.d_id
             AND b.vehicle_name = t2.vehicle_name AND b.interval_id = t2.interval_id
            LEFT JOIN tier3_speed t3
              ON {tier3_join_condition_b_t3}
            LEFT JOIN tier4_speed t4
              ON b.o_id = t4.o_id AND b.d_id = t4.d_id
             AND b.vehicle_name = t4.vehicle_name AND b.interval_id = t4.interval_id
            LEFT JOIN tier5_speed t5
              ON b.vehicle_name = t5.vehicle_name AND b.interval_id = t5.interval_id
            LEFT JOIN tier6_per_interval t6iv
              ON b.interval_id = t6iv.interval_id
        ),
        final AS (
            SELECT a.*,
                CASE
                    WHEN a.t1_duration IS NOT NULL THEN a.t1_duration
                    WHEN a.t2_speed IS NOT NULL AND a.t2_speed > 0 THEN (a.dist / a.t2_speed)::int
                    WHEN a.t3_speed IS NOT NULL AND a.t3_speed > 0 THEN (a.dist / a.t3_speed)::int
                    WHEN a.t4_speed IS NOT NULL AND a.t4_speed > 0 THEN (a.dist / a.t4_speed)::int
                    WHEN a.t5_speed IS NOT NULL AND a.t5_speed > 0 THEN (a.dist / a.t5_speed)::int
                    WHEN a.t6_speed_iv IS NOT NULL AND a.t6_speed_iv > 0 THEN (a.dist / a.t6_speed_iv)::int
                    WHEN a.t6_speed_all IS NOT NULL AND a.t6_speed_all > 0 THEN (a.dist / a.t6_speed_all)::int
                    ELSE NULL
                END AS duration_final,
                CASE
                    WHEN a.t1_duration IS NOT NULL THEN a.t1_calc_method
                    WHEN a.t2_speed IS NOT NULL AND a.t2_speed > 0 THEN 'interpolated_segment_neighbors'
                    WHEN a.t3_speed IS NOT NULL AND a.t3_speed > 0 THEN 'interpolated_same_segment_other_interval'
                    WHEN a.t4_speed IS NOT NULL AND a.t4_speed > 0 THEN 'estimated_route_avg'
                    WHEN a.t5_speed IS NOT NULL AND a.t5_speed > 0 THEN 'estimated_vehicle_avg'
                    WHEN a.t6_speed_iv IS NOT NULL AND a.t6_speed_iv > 0 THEN 'estimated_citywide'
                    WHEN a.t6_speed_all IS NOT NULL AND a.t6_speed_all > 0 THEN 'estimated_citywide'
                    ELSE NULL
                END AS calc_method_final
            FROM assembled a
        )

        SELECT
            ROW_NUMBER() OVER () AS gid,
            f.o_id, f.d_id,
            f.interval_id, f.interval_name, f.interval_start,
            f.vehicle_name,
            f.dist,
            f.duration_final AS duration,
            (f.dist / NULLIF(f.duration_final, 0)) * 3.6 AS speed,
            f.calc_method_final AS calc_method,
            COALESCE(f.t1_n_samples, 0) AS n_samples,
            f.geom::geometry(LINESTRING,4326) AS geom
        FROM final f
        WHERE f.duration_final IS NOT NULL
        """  # nosec B608 — interpolated fragments are fixed SQL constants chosen from a closed if/else; all user values use %(name)s binds

    # ---- main ---------------------------------------------------------
    def processAlgorithm(self, parameters, context, feedback):
        ensure_deps(show_ui=True)

        rl_conn = self.parameterAsConnectionName(parameters, self.RL_CONN, context)
        project_id = self.parameterAsString(parameters, self.PROJECT_ID, context)
        out_gpkg = self.parameterAsFileOutput(parameters, self.OUT_GPKG, context)
        overwrite = self.parameterAsBool(parameters, self.OVERWRITE, context)
        include_raw = self.parameterAsBool(parameters, self.INCLUDE_RAW, context)
        include_qa = self.parameterAsBool(parameters, self.INCLUDE_QA, context)
        fallback = self.parameterAsInt(parameters, self.FALLBACK_HEADWAY, context)
        if fallback is not None and fallback <= 0:
            fallback = None

        # Read stop-extraction parameters (advanced; defaults preserved if untouched).
        stop_params = read_stop_params_from_algorithm(self, parameters, context)
        feedback.pushInfo(f"Stop-extraction parameters: {stop_params.describe()}")

        if os.path.exists(out_gpkg):
            if overwrite:
                os.remove(out_gpkg)
            else:
                raise Exception(f"Output already exists: {out_gpkg}")

        # If file exists (overwrite=False), verify it's a real gpkg/sqlite
        if os.path.exists(out_gpkg):
            con = sqlite3.connect(out_gpkg)
            try:
                cur = con.cursor()
                cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='gpkg_spatial_ref_sys'")
                ok = cur.fetchone() is not None
            finally:
                con.close()
            if not ok:
                raise Exception("Output file exists but is not a valid GeoPackage. Delete it or enable overwrite.")

        pg_conn = self._psycopg2_conn_from_qgis_connection(rl_conn)
        feedback.pushInfo("Connecting to RouteLab Postgres (read-only)…")

        # Resolve setting_id — project_id bound via params, never interpolated.
        setting_df = pd.read_sql_query(
            self._sql_setting_id(),
            con=pg_conn,
            params={"pid": project_id},
        )
        if setting_df.empty:
            raise Exception(f"Could not find setting_id for project_id={project_id}")
        setting_id = str(setting_df.iloc[0]["id"])
        feedback.pushInfo(f"Using setting_id={setting_id}")

        # Single bind dict used by every subsequent query. psycopg2 ignores
        # keys a given query doesn't reference, so reusing one dict is safe.
        q = {
            "sid":                 setting_id,
            "pid":                 project_id,
            "dbscan_eps_m":        stop_params.dbscan_eps_m,
            "dbscan_minpoints":    stop_params.dbscan_minpoints,
            "snap_max_m":          stop_params.snap_max_m,
            "terminal_m":          stop_params.terminal_m,
            "cell_m":              stop_params.cell_m,
            "stop_trip_buffer_m":  stop_params.stop_trip_buffer_m,
            "min_stop_spacing_m":  stop_params.min_stop_spacing_m,
            "trackpoint_buffer_m": stop_params.trackpoint_buffer_m,
        }

        # Base tables
        vehicles_sql  = f"WITH {self._cte_vehicles()}  SELECT * FROM transit_vehicles"   # nosec B608 — composes parameterized builder
        agencies_sql  = f"WITH {self._cte_agencies()}  SELECT * FROM transit_agencies"   # nosec B608 — composes parameterized builder
        intervals_sql = f"WITH {self._cte_intervals()} SELECT * FROM transit_intervals"  # nosec B608 — composes parameterized builder
        terminals_sql = f"WITH {self._cte_terminals()} SELECT * FROM transit_terminals"  # nosec B608 — composes parameterized builder
        trips_sql     = f"WITH {self._cte_trips()}     SELECT * FROM transit_trips"      # nosec B608 — composes parameterized builder

        vehicles  = pd.read_sql_query(vehicles_sql,  con=pg_conn, params=q)
        agencies  = pd.read_sql_query(agencies_sql,  con=pg_conn, params=q)
        intervals = pd.read_sql_query(intervals_sql, con=pg_conn, params=q)
        terminals = gpd.read_postgis(terminals_sql,  con=pg_conn, geom_col="geom", params=q)
        trips     = gpd.read_postgis(trips_sql,      con=pg_conn, geom_col="geom", params=q)

        # trips_view
        feedback.pushInfo("Computing transit.trips_view …")
        trips_view = gpd.read_postgis(self._sql_trips_view(), con=pg_conn, geom_col="geom", params=q)

        # stops / QA
        feedback.pushInfo("Computing transit stops (DBSCAN → spaced → snapped) …")
        stops_auto = gpd.read_postgis(self._sql_stops_auto(), con=pg_conn, geom_col="geom", params=q)
        stops      = gpd.read_postgis(self._sql_stops(),      con=pg_conn, geom_col="geom", params=q)
        stop_clusters = None
        if include_qa:
            stop_clusters = gpd.read_postgis(self._sql_stop_clusters(), con=pg_conn, geom_col="centroid", params=q)

        # trip_stops_sequence + od_stats
        feedback.pushInfo("Computing trip stop sequence …")
        trip_stops_sequence = pd.read_sql_query(self._sql_trip_stops_sequence(), con=pg_conn, params=q)

        feedback.pushInfo("Computing OD stats …")
        od_stats = gpd.read_postgis(
            self._sql_od_stats(stop_params.distinguish_speeds_by_vehicle),
            con=pg_conn, geom_col="geom", params=q,
        )

        # Headway estimation (pandas) using frequency instances
        feedback.pushInfo("Estimating headways …")
        freq = pd.read_sql_query(
            "select * from v_frequency_instances_ext "
            "where setting_id = %(sid)s and deleted_at is null",
            con=pg_conn,
            params=q,
        )
        # Filter + aggregate (same logic as RL2SDI)
        if not freq.empty:
            freq = freq.copy()
            # Some deployments may have null origin/destination
            freq["origin"] = freq.get("origin")
            freq["destination"] = freq.get("destination")
            mask = True
            if "origin" in freq.columns:
                mask = mask & (~freq.origin.astype(str).str.lower().str.contains("test", na=False))
            if "destination" in freq.columns:
                mask = mask & (~freq.destination.astype(str).str.lower().str.contains("test", na=False))
            if "status" in freq.columns:
                mask = mask & (freq.status == "finished")
            freq2 = (
                freq.loc[mask, ["trip_id", "interval", "avg_headway_sec"]]
                .groupby(["trip_id", "interval"], as_index=False)
                .agg({"avg_headway_sec": "mean"})
            )
        else:
            freq2 = pd.DataFrame(columns=["trip_id", "interval", "avg_headway_sec"])

        intervals_map = intervals[["gid", "observer_id"]].rename(columns={"gid": "interval_gid", "observer_id": "interval"})
        trips_df = trips.drop(columns=["geom"]).copy()
        # cross join
        trips_x = trips_df.merge(intervals_map, how="cross")
        trips_x = trips_x.merge(
            freq2,
            left_on=["observer_id", "interval"],
            right_on=["trip_id", "interval"],
            how="left",
        )
        # Fill within route+interval
        trips_x["avg_headway_sec"] = trips_x.groupby(["observer_route_id", "interval"], dropna=False)["avg_headway_sec"].transform(lambda s: s.ffill().bfill())
        trips_x["avg_headway_agency_interval"] = trips_x.groupby(["agency_id", "interval"], dropna=False)["avg_headway_sec"].transform("mean")
        trips_x["ratio_headway"] = trips_x["avg_headway_sec"] / trips_x["avg_headway_agency_interval"]
        trips_x["avg_ratio_headway"] = trips_x.groupby(["observer_route_id"], dropna=False)["ratio_headway"].transform("mean")
        trips_x["headway_estimation_method"] = trips_x["avg_headway_sec"].isna().map(lambda x: "from_similar_agency_interval" if x else "from_own_freq_surveys")
        trips_x["final_headway"] = (trips_x["avg_headway_sec"].fillna(trips_x["avg_ratio_headway"] * trips_x["avg_headway_agency_interval"]).apply(pd.to_numeric, errors="coerce")).apply(lambda v: None if pd.isna(v) else float(v))
        # floor
        trips_x["final_headway"] = trips_x["final_headway"].apply(lambda v: None if pd.isna(v) else float(int(v // 1)))
        if fallback is not None:
            miss = trips_x["final_headway"].isna()
            if miss.any():
                trips_x.loc[miss, "final_headway"] = float(fallback)
                trips_x.loc[miss, "headway_estimation_method"] = "from_user_default"

        trips_intervals = trips_x[["gid", "final_headway", "interval_gid", "headway_estimation_method"]].rename(
            columns={"gid": "trip_id", "final_headway": "headway_secs", "interval_gid": "interval_id"}
        )

        # Raw layers (optional)
        raw_stops = raw_trackpoints = raw_onboard = raw_identification = raw_frequency = None
        if include_raw:
            feedback.pushInfo("Exporting raw layers …")
            raw_stops = gpd.read_postgis(
                f"""
                WITH {self._cte_raw_stops()}
                SELECT *
                FROM raw_stops
                """,  # nosec B608 — composes parameterized builder
                con=pg_conn,
                geom_col="geom",
                params=q,
            )
            raw_trackpoints = gpd.read_postgis(
                """
                WITH
                r AS (SELECT id AS route_id FROM routes WHERE setting_id = %(sid)s AND deleted_at IS NULL),
                t AS (SELECT id AS trip_id FROM trips JOIN r ON trips.route_id = r.route_id WHERE deleted_at IS NULL),
                oi AS (
                    SELECT id AS oi_id, status, valid
                    FROM onboard_instances JOIN t ON onboard_instances.trip_id = t.trip_id
                    WHERE deleted_at IS NULL
                )
                SELECT
                    ROW_NUMBER() OVER () AS gid,
                    tp.timestamp AS timestamp,
                    tp.onboard_instance_id::text AS onboard_instance_id,
                    oi.status AS onboard_instance_status,
                    oi.valid AS onboard_instance_valid,
                    tp.geometry::geometry(POINT,4326) AS geom
                FROM onboard_instance_track_points tp
                JOIN oi ON tp.onboard_instance_id = oi.oi_id
                WHERE tp.deleted_at IS NULL AND tp.geometry IS NOT NULL AND ST_IsValid(tp.geometry)
                """,
                con=pg_conn,
                geom_col="geom",
                params=q,
            )
            raw_onboard = gpd.read_postgis(
                "select * from v_onboard_instances_ext "
                "where project_id = %(pid)s and deleted_at is null "
                "and geometry is not null and ST_IsValid(geometry)",
                con=pg_conn,
                geom_col="geometry",
                params=q,
            )
            raw_identification = gpd.read_postgis(
                """
                select *
                from v_terminal_trips_ext
                where project_id = %(pid)s
                and deleted_at is null
                and geometry is not null
                and ST_IsValid(geometry)
                """,
                con=pg_conn,
                geom_col="geometry",
                params=q,
            )

            raw_frequency = gpd.read_postgis(
                """
                select *
                from v_frequency_instances_ext
                where project_id = %(pid)s
                and deleted_at is null
                and geometry is not null
                and ST_IsValid(geometry)
                """,
                con=pg_conn,
                geom_col="geometry",
                params=q,
            )

        pg_conn.close()
        # Best-effort cleanup: engine disposal can fail if already closed. Non-fatal.
        with contextlib.suppress(Exception):
            getattr(self, "_dbapi_engine", None) and self._dbapi_engine.dispose()

        # Write GeoPackage
        feedback.pushInfo("Writing GeoPackage…")

        _write_gdf(out_gpkg, "transit_terminals", terminals)

        with sqlite3.connect(out_gpkg) as conn:
            meta = pd.DataFrame(
                [
                    {
                        "schema_name": "tfc_rl_sdi_gpkg",
                        "schema_version": "0.2",
                        "exported_at": datetime.now().isoformat(timespec="seconds"),
                        "routelab_project_id": project_id,
                        "generator": "TfC Tools: Export RouteLab to GeoPackage (read-only)",
                        "notes": None,
                    }
                ]
            )
            _write_df(conn, "tfc_meta", meta)
            _write_df(conn, "transit_vehicles", vehicles)
            _write_df(conn, "transit_agencies", agencies)
            _write_df(conn, "transit_intervals", intervals.rename(columns={"name": "interval_name"}) if "name" in intervals.columns and "interval_name" not in intervals.columns else intervals)
            _write_df(conn, "transit_trips_intervals", trips_intervals)
            _write_df(conn, "transit_trip_stops_sequence", trip_stops_sequence)

        _write_gdf(out_gpkg, "transit_trips", trips)
        _write_gdf(out_gpkg, "transit_trips_view", trips_view)
        _write_gdf(out_gpkg, "transit_stops", stops)
        _write_gdf(out_gpkg, "transit_od_stats", od_stats)

        if include_qa:
            if stop_clusters is not None:
                if stop_clusters.geometry.name != "centroid" and "centroid" in stop_clusters.columns:
                    stop_clusters = stop_clusters.set_geometry("centroid")
                _write_gdf(out_gpkg, "transit__stop_clusters", stop_clusters)
            _write_gdf(out_gpkg, "transit_stops_auto", stops_auto)

        if include_raw:
            if raw_stops is not None:
                _write_gdf(out_gpkg, "raw_stops", raw_stops)
            if raw_trackpoints is not None:
                _write_gdf(out_gpkg, "raw_trackpoints", raw_trackpoints)
            if raw_onboard is not None:
                _write_gdf(out_gpkg, "raw_onboard_instances", raw_onboard)
            if raw_identification is not None:
                _write_gdf(out_gpkg, "raw_identification_instances", raw_identification)
            if raw_frequency is not None:
                _write_gdf(out_gpkg, "raw_frequency_instances", raw_frequency)

        feedback.pushInfo(f"GeoPackage written: {out_gpkg}")
        return {"OUT_GPKG": out_gpkg}
