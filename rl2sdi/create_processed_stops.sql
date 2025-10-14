/*==============================================================================
   STOP PROCESSING PIPELINE (compact; minimal persistent objects)
   Keeps only: transit._stop_clusters (QA) and transit.stops (final)
==============================================================================*/

CREATE EXTENSION IF NOT EXISTS postgis;
-- NOTE: pgcrypto no longer needed

/*----------------------------------------------------------------------------
 STEP 1: Build clusters directly from raw.stops (no extra temp views)
----------------------------------------------------------------------------*/
DROP MATERIALIZED VIEW IF EXISTS transit._stop_clusters CASCADE;
CREATE MATERIALIZED VIEW transit._stop_clusters AS
WITH pts AS (
  SELECT
    s.gid AS src_id,
    NULLIF(TRIM(s.name), '') AS raw_name,
    s.geom,
    ST_Transform(s.geom, 3857) AS g3857
  FROM raw.stops s
  WHERE s.geom IS NOT NULL
),
clusters AS (
  SELECT
    *,
    ST_ClusterDBSCAN(g3857, eps := 150, minpoints := 3) OVER () AS cluster_id
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
LEFT JOIN name_mode n USING (cluster_id);

CREATE INDEX IF NOT EXISTS _stop_clusters_centroid_gix ON transit._stop_clusters USING GIST (centroid);

/*----------------------------------------------------------------------------
 STEP 2: One-shot pipeline (CTEs only) → stops_auto (few rows; easy to inspect)
   - 250 m spacing
   - snap to trips (≤ 20 m; tune SNAP_MAX_M below)
   - terminal flag (within 75 m of trip start/end)
   - direction clone 0/1 from local bearing
----------------------------------------------------------------------------*/
DROP MATERIALIZED VIEW IF EXISTS transit.stops_auto CASCADE;
CREATE MATERIALIZED VIEW transit.stops_auto AS
WITH
const AS (
  SELECT 250::double precision AS MIN_SPACING_M,
         30::double precision  AS SNAP_MAX_M,
         75::double precision  AS TERMINAL_M
),
c AS (
  SELECT sc.cluster_id, sc.mode_name, sc.n_points, sc.centroid,
         ST_Transform(sc.centroid, 3857) AS c3857
  FROM transit._stop_clusters sc
),
-- replace your spaced AS (...) with this block
spaced AS (
  WITH params AS (SELECT 120.0::double precision AS cell_m),  -- tweak spacing here
  cells AS (
    SELECT
      c.*,
      (SELECT cell_m FROM params) AS cell_m,
      -- snap 3857 coords to a grid
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
         s.geom     AS orig_geom, 
         t.gid      AS trip_gid,
         ST_Transform(
           ST_ClosestPoint(ST_Transform(t.geom,3857), s.c3857), 4326
         ) AS geom,
         ST_Distance(s.c3857, ST_Transform(t.geom,3857)) AS dist_m,
         t.geom AS trip_geom
  FROM spaced s
  JOIN LATERAL (
    SELECT gid, geom
    FROM transit.trips
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
        ST_LineInterpolatePoint(ST_Transform(trip_geom,3857),
          GREATEST(ST_LineLocatePoint(ST_Transform(trip_geom,3857), ST_Transform(geom,3857)) - 0.0005, 0)
        ),
        ST_LineInterpolatePoint(ST_Transform(trip_geom,3857),
          LEAST(ST_LineLocatePoint(ST_Transform(trip_geom,3857), ST_Transform(geom,3857)) + 0.0005, 1)
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
  mode_name                                  AS stop_name,
  CASE stop_type WHEN 'Terminal' THEN 2 ELSE 0 END AS location_type,
  mode_name || CASE WHEN dir_bin = 0 THEN ' (Dir 0)' ELSE ' (Dir 1)' END AS stop_desc,
  stop_type,
  cluster_id                                 AS cluster_id,
  dir_bin                                    AS "double",
  ST_X(geom)                                 AS stop_lon,
  ST_Y(geom)                                 AS stop_lat,
  geom
FROM final;

CREATE INDEX IF NOT EXISTS stops_auto_name_idx ON transit.stops_auto (stop_name);
CREATE INDEX IF NOT EXISTS stops_auto_geom_gix  ON transit.stops_auto USING GIST (geom);
CREATE INDEX IF NOT EXISTS trips_geom_gix       ON transit.trips USING GIST (geom);

/*----------------------------------------------------------------------------
 STEP 3: Final table (stop_id mirrors gid) + load
----------------------------------------------------------------------------*/
DROP TABLE IF EXISTS transit.stops CASCADE;
CREATE TABLE transit.stops (
  gid           SERIAL PRIMARY KEY,
  stop_id       TEXT GENERATED ALWAYS AS (gid::text) STORED,
  stop_name     TEXT,
  stop_desc     TEXT,
  location_type INT,
  "double"      INT,
  stop_lon      DOUBLE PRECISION,
  stop_lat      DOUBLE PRECISION,
  geom          geometry(Point, 4326)
);

CREATE UNIQUE INDEX transit_stops_stop_id_uidx ON transit.stops (stop_id);
CREATE INDEX transit_stops_geom_gist ON transit.stops USING GIST (geom);

INSERT INTO transit.stops (stop_name, stop_desc, location_type, "double", stop_lon, stop_lat, geom)
SELECT stop_name, stop_desc, location_type, "double", stop_lon, stop_lat, geom
FROM transit.stops_auto;
