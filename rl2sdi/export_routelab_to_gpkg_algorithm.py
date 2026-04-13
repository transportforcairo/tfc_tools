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
        try:
            p.setFlags(p.flags() | QgsProcessingParameterNumber.FlagAdvanced)
        except Exception:
            pass
        self.addParameter(p)

    # ---- SQL builders -------------------------------------------------
    @staticmethod
    def _sql_setting_id(project_id: str) -> str:
        # Use ORDER BY id DESC to avoid relying on created_at.
        pid = project_id.replace("'", "''")
        return (
            "SELECT id FROM settings "
            f"WHERE project_id='{pid}' AND deleted_at IS NULL "
            "ORDER BY id DESC LIMIT 1"
        )

    @staticmethod
    def _cte_vehicles(sid: str) -> str:
        return f"""
        rl_agencies AS (
            SELECT * FROM agencies WHERE setting_id='{sid}' AND deleted_at IS NULL
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
    def _cte_agencies(sid: str) -> str:
        return f"""
        rl_agencies AS (
            SELECT * FROM agencies WHERE setting_id='{sid}' AND deleted_at IS NULL
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
    def _cte_terminals(sid: str) -> str:
        return f"""
        transit_terminals AS (
            SELECT
                ROW_NUMBER() OVER (ORDER BY t.id) AS gid,
                t.name,
                t.id::text AS observer_id,
                t.geometry::geometry(MULTIPOLYGON,4326) AS geom
            FROM terminals t
            WHERE t.setting_id='{sid}'
              AND t.deleted_at IS NULL
              AND t.status = 'accepted'
              AND t.geometry IS NOT NULL
              AND ST_IsValid(t.geometry)
              AND (t.name IS NULL OR lower(t.name) NOT LIKE '%test%')
        )
        """

    @staticmethod
    def _cte_intervals(sid: str) -> str:
        return f"""
        transit_intervals AS (
            SELECT
                ROW_NUMBER() OVER (ORDER BY i.start, i.end, i.id) AS gid,
                i.start::time AS start_time,
                i.end::time AS end_time,
                i.id::text AS observer_id,
                i.name::text AS name,
                (i.deleted_at IS NULL) AS active
            FROM intervals i
            JOIN frequency_setting_intervals fsi ON i.id = fsi.interval_id
            JOIN frequency_settings fs ON fsi.frequency_setting_id = fs.id
            WHERE fs.setting_id='{sid}'
        )
        """

    @staticmethod
    def _cte_trips(sid: str) -> str:
        # trips with stable gid ordering
        return f"""
        {ExportRouteLabToGeoPackageAlgorithm._cte_agencies(sid)},
        {ExportRouteLabToGeoPackageAlgorithm._cte_terminals(sid)},
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
            WHERE v.setting_id='{sid}'
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
            WHERE (o.name IS NULL OR lower(o.name) NOT LIKE '%test%')
              AND (d.name IS NULL OR lower(d.name) NOT LIKE '%test%')
        )
        """

    @staticmethod
    def _sql_trips_view(sid: str) -> str:
        return f"""
        WITH
        {ExportRouteLabToGeoPackageAlgorithm._cte_trips(sid)}
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
        """

    @staticmethod
    def _cte_raw_stops(sid: str) -> str:
        # raw stop events from onboard_instance_stops
        return f"""
        r AS (
            SELECT id AS route_id FROM routes WHERE setting_id='{sid}' AND deleted_at IS NULL
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
    def _sql_stop_clusters(sid: str) -> str:
        return f"""
        WITH
        {ExportRouteLabToGeoPackageAlgorithm._cte_raw_stops(sid)},
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
            SELECT *, ST_ClusterDBSCAN(g3857, eps := 150, minpoints := 3) OVER () AS cluster_id
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
        """

    @staticmethod
    def _sql_stops_auto(sid: str) -> str:
        # SELECT-only version of transit.stops_auto (from create_processed_stops.sql),
        # but sourced from RouteLab (raw_stops + v_trips_ext-derived transit_trips CTE).
        return f"""
        WITH
        {ExportRouteLabToGeoPackageAlgorithm._cte_trips(sid)},
        stop_clusters AS ({ExportRouteLabToGeoPackageAlgorithm._sql_stop_clusters(sid)}),
        const AS (
            SELECT 250::double precision AS MIN_SPACING_M,
                   30::double precision  AS SNAP_MAX_M,
                   75::double precision  AS TERMINAL_M
        ),
        c AS (
            SELECT sc.cluster_id, sc.mode_name, sc.n_points, sc.centroid,
                   ST_Transform(sc.centroid, 3857) AS c3857
            FROM stop_clusters sc
        ),
        spaced AS (
            WITH params AS (SELECT 120.0::double precision AS cell_m),
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
                CASE WHEN (bearing + 360)::int % 360 BETWEEN 45 AND 225 THEN 1 ELSE 0 END AS dir_bin
            FROM bearing_calc
        )
        SELECT
            mode_name AS stop_name,
            CASE stop_type WHEN 'Terminal' THEN 2 ELSE 0 END AS location_type,
            mode_name || CASE WHEN dir_bin = 0 THEN ' (Dir 0)' ELSE ' (Dir 1)' END AS stop_desc,
            stop_type,
            cluster_id,
            dir_bin AS "double",
            ST_X(geom) AS stop_lon,
            ST_Y(geom) AS stop_lat,
            geom
        FROM final
        """

    @staticmethod
    def _sql_stops(sid: str) -> str:
        return f"""
        WITH stops_auto AS ({ExportRouteLabToGeoPackageAlgorithm._sql_stops_auto(sid)})
        SELECT
            ROW_NUMBER() OVER (ORDER BY stop_name, cluster_id, "double") AS gid,
            (ROW_NUMBER() OVER (ORDER BY stop_name, cluster_id, "double"))::text AS stop_id,
            stop_name,
            stop_desc,
            location_type,
            "double",
            stop_lon,
            stop_lat,
            geom
        FROM stops_auto
        """

    @staticmethod
    def _sql_trip_stops_sequence(sid: str) -> str:
        return f"""
        WITH
        {ExportRouteLabToGeoPackageAlgorithm._cte_trips(sid)},
        stops AS ({ExportRouteLabToGeoPackageAlgorithm._sql_stops(sid)}),
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
            JOIN transit_trips t ON ST_DWithin(t.geom::geography, s.geom::geography, 1::real)
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
        WHERE distance_from_prev >= 100 OR distance_from_prev IS NULL
        """

    @staticmethod
    def _sql_od_stats(sid: str) -> str:
        # SELECT-only version of od_stats.sql adapted to RouteLab base tables.
        return f"""
        WITH
        {ExportRouteLabToGeoPackageAlgorithm._cte_trips(sid)},
        {ExportRouteLabToGeoPackageAlgorithm._cte_intervals(sid)},
        stops AS ({ExportRouteLabToGeoPackageAlgorithm._sql_stops(sid)}),
        trip_stops_sequence AS ({ExportRouteLabToGeoPackageAlgorithm._sql_trip_stops_sequence(sid)}),

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
                        LEAD(ts.distance, 1) OVER (PARTITION BY ts.t_id ORDER BY ts.t_id, ts.stop_sequence) - ts.distance AS dist,
                        ST_LineSubstring(t.geom, ts.distance_frac,
                            LEAD(ts.distance_frac, 1) OVER (PARTITION BY ts.t_id ORDER BY ts.t_id, ts.stop_sequence)
                        ) AS geom
                    FROM trip_stops_sequence ts
                    JOIN transit_trips t ON ts.t_id = t.gid AND ts.distance_frac < 1
                ) x
                WHERE d_id IS NOT NULL
            ) y
            GROUP BY o_id, d_id, vehicle_name, geom
        ),

        -- valid onboard instances for this setting
        r AS (SELECT id AS route_id FROM routes WHERE setting_id='{sid}' AND deleted_at IS NULL),
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
            JOIN stops s ON ST_DWithin(tr.geom::geography, s.geom::geography, 30)
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
        avg_durations AS (
            SELECT
                ROW_NUMBER() OVER () AS gid,
                o_id,
                d_id,
                i.gid AS interval_id,
                i.name AS interval_name,
                vehicle_name,
                MIN(tmp.from_id) AS from_id,
                MIN(tmp.to_id) AS to_id,
                MIN(i.start_time) AS interval_start,
                AVG(time_diff)::int AS duration
            FROM (
                SELECT *, EXTRACT(EPOCH FROM (d_time - o_time)) AS time_diff
                FROM od_timestamps
            ) tmp
            JOIN transit_intervals i ON (tmp.o_time::time BETWEEN i.start_time AND i.end_time) AND i.active
            WHERE time_diff > 0
            GROUP BY o_id, d_id, i.gid, i.name, vehicle_name
        )

        SELECT
            ROW_NUMBER() OVER () AS gid,
            seg.o_id,
            seg.d_id,
            avg_durations.interval_id,
            avg_durations.interval_name,
            avg_durations.interval_start,
            avg_durations.vehicle_name,
            seg.dist,
            avg_durations.duration,
            (seg.dist/NULLIF(avg_durations.duration,0))*3.6 AS speed,
            seg.geom::geometry(LINESTRING,4326) AS geom
        FROM od_segments seg
        JOIN avg_durations ON seg.o_id = avg_durations.o_id AND seg.d_id = avg_durations.d_id AND seg.vehicle_name = avg_durations.vehicle_name
        """

    # ---- main ---------------------------------------------------------
    def processAlgorithm(self, parameters, context, feedback):
        ensure_deps(show_ui=True)

        rl_conn = self.parameterAsConnectionName(parameters, self.RL_CONN, context)
        project_id = self.parameterAsString(parameters, self.PROJECT_ID, context)
        safe_project_id = project_id.replace("'", "''")
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

        # Resolve setting_id
        setting_df = pd.read_sql_query(self._sql_setting_id(project_id), con=pg_conn)
        if setting_df.empty:
            raise Exception(f"Could not find setting_id for project_id={project_id}")
        setting_id = str(setting_df.iloc[0]["id"]) 
        feedback.pushInfo(f"Using setting_id={setting_id}")

        # Base tables
        vehicles_sql = f"WITH {self._cte_vehicles(setting_id)} SELECT * FROM transit_vehicles"
        agencies_sql = f"WITH {self._cte_agencies(setting_id)} SELECT * FROM transit_agencies"
        intervals_sql = f"WITH {self._cte_intervals(setting_id)} SELECT * FROM transit_intervals"
        terminals_sql = f"WITH {self._cte_terminals(setting_id)} SELECT * FROM transit_terminals"
        trips_sql = f"WITH {self._cte_trips(setting_id)} SELECT * FROM transit_trips"

        vehicles = pd.read_sql_query(vehicles_sql, con=pg_conn)
        agencies = pd.read_sql_query(agencies_sql, con=pg_conn)
        intervals = pd.read_sql_query(intervals_sql, con=pg_conn)
        terminals = gpd.read_postgis(terminals_sql, con=pg_conn, geom_col="geom")
        trips = gpd.read_postgis(trips_sql, con=pg_conn, geom_col="geom")

        # trips_view
        feedback.pushInfo("Computing transit.trips_view …")
        trips_view = gpd.read_postgis(self._sql_trips_view(setting_id), con=pg_conn, geom_col="geom")

        # stops / QA
        feedback.pushInfo("Computing transit stops (DBSCAN → spaced → snapped) …")
        stops_auto = gpd.read_postgis(self._sql_stops_auto(setting_id), con=pg_conn, geom_col="geom")
        stops = gpd.read_postgis(self._sql_stops(setting_id), con=pg_conn, geom_col="geom")
        stop_clusters = None
        if include_qa:
            stop_clusters = gpd.read_postgis(self._sql_stop_clusters(setting_id), con=pg_conn, geom_col="centroid")

        # trip_stops_sequence + od_stats
        feedback.pushInfo("Computing trip stop sequence …")
        trip_stops_sequence = pd.read_sql_query(self._sql_trip_stops_sequence(setting_id), con=pg_conn)

        feedback.pushInfo("Computing OD stats …")
        od_stats = gpd.read_postgis(self._sql_od_stats(setting_id), con=pg_conn, geom_col="geom")

        # Headway estimation (pandas) using frequency instances
        feedback.pushInfo("Estimating headways …")
        freq = pd.read_sql_query(
            f"select * from v_frequency_instances_ext where setting_id='{setting_id}' and deleted_at is null",
            con=pg_conn,
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
                WITH {self._cte_raw_stops(setting_id)}
                SELECT *
                FROM raw_stops
                """,
                con=pg_conn,
                geom_col="geom",
            )
            raw_trackpoints = gpd.read_postgis(
                f"""
                WITH
                r AS (SELECT id AS route_id FROM routes WHERE setting_id='{setting_id}' AND deleted_at IS NULL),
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
            )
            raw_onboard = gpd.read_postgis(
                f"select * from v_onboard_instances_ext where project_id='{safe_project_id}' and deleted_at is null and geometry is not null and ST_IsValid(geometry)",
                con=pg_conn,
                geom_col="geometry",
            )
            raw_identification = gpd.read_postgis(
                f"""
                select *
                from v_terminal_trips_ext
                where project_id='{safe_project_id}'
                and deleted_at is null
                and geometry is not null
                and ST_IsValid(geometry)
                """,
                con=pg_conn,
                geom_col="geometry",
            )

            raw_frequency = gpd.read_postgis(
                f"""
                select *
                from v_frequency_instances_ext
                where project_id='{safe_project_id}'
                and deleted_at is null
                and geometry is not null
                and ST_IsValid(geometry)
                """,
                con=pg_conn,
                geom_col="geometry",
            )

        pg_conn.close()
        try:
            getattr(self, "_dbapi_engine", None) and self._dbapi_engine.dispose()
        except Exception:
            pass

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
