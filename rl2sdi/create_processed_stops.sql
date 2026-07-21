/*==============================================================================
   STOP PROCESSING PIPELINE  —  direction-aware (v2)
   ----------------------------------------------------------------------------
   Rewritten to build a directional, paired, cleanly-named stops network from
   raw observed stop events, instead of the previous "cluster → grid-thin to one
   point → compass-bin the direction" approach.

   KEY CHANGE: direction is DATA, not geometry.
     Each raw stop was observed on a specific onboard instance, which belongs to
     a specific trip, which has a known direction_id. So we attribute every raw
     stop to its surveyed direction and cluster PER DIRECTION. Where both
     directions are observed at one physical location, two stops are produced and
     linked as a pair (shared name + stop_group). This reproduces the human
     pattern (≈50% of stops paired at ~10 m, sharing a name).

   REQUIRED SOURCE COLUMNS (SDI schema — confirm names against your DB):
     raw.stops(gid, name, onboard_instance_observer_id, geom(Point,4326))
       -- 'name' is the raw operator-tagged stop name; if your column is
       --  'raw_name', change the two references below.
     raw.onboard_instances(id, trip_id)
       -- id            = onboard instance observer id (matches raw.stops link)
       -- trip_id       = observer/external trip id (matches transit.trips.observer_id)
     transit.trips(gid, observer_id, direction_id, geom(LineString,4326))

   PLACEHOLDERS (Python str.format substitution, same mechanism as before):
     {dbscan_eps_m} {dbscan_minpoints} {snap_max_m} {terminal_m} {pair_max_m}
     {rescue_gap_m}

   TWO-TIER STOP CREATION:
     Tier A — confident stops: DBSCAN(minpoints) clusters, per direction.
     Tier B — coverage rescue: isolated on-route raw stops (DBSCAN noise) kept
       as stops when no Tier-A stop is within {rescue_gap_m}, so sparsely-surveyed
       corridors still get accessibility. Straight-line gap (not along-corridor).

   Validated against the Sousse export with:
     eps=40 m, minpoints=2, snap_max=40 m, pair_max=40 m, terminal=75 m
   giving ~45% pairing, ~1.0 same-name-within-pair, ~11 m median match to the
   human layer, and zero name-whitespace defects.

   SYNTHETIC / MANUAL STOPS: rare hand-added stops (to allow boarding/alighting
   where no raw data exists) live in transit.stops with stop_source='manual' and
   are PRESERVED across rebuilds — see STEP 3.
==============================================================================*/

CREATE EXTENSION IF NOT EXISTS postgis;

/*----------------------------------------------------------------------------
 STEP 1: Direction-aware clusters  (QA / inspection matview)
   - attribute each raw stop with its surveyed direction_id
   - DBSCAN PER direction so the two carriageways are clustered independently
----------------------------------------------------------------------------*/
DROP MATERIALIZED VIEW IF EXISTS transit._stop_clusters CASCADE;
CREATE MATERIALIZED VIEW transit._stop_clusters AS
WITH attributed AS (
  SELECT
    s.gid                              AS src_id,
    NULLIF(BTRIM(s.name), '')          AS raw_name,     -- change to s.raw_name if needed
    t.direction_id                     AS direction_id,
    ST_Transform(s.geom, 3857)         AS g3857
  FROM raw.stops s
  JOIN raw.onboard_instances oi ON s.onboard_instance_observer_id = oi.id
  JOIN transit.trips         t  ON oi.trip_id = t.observer_id
  WHERE s.geom IS NOT NULL
    AND t.direction_id IS NOT NULL
),
clustered AS (
  SELECT
    a.*,
    ST_ClusterDBSCAN(g3857, eps := {dbscan_eps_m}, minpoints := {dbscan_minpoints})
      OVER (PARTITION BY direction_id) AS cluster_id
  FROM attributed a
),
valid AS (
  SELECT * FROM clustered WHERE cluster_id IS NOT NULL
),
name_counts AS (
  SELECT direction_id, cluster_id, raw_name, COUNT(*) AS cnt
  FROM valid WHERE raw_name IS NOT NULL
  GROUP BY direction_id, cluster_id, raw_name
),
name_mode AS (
  -- modal raw name per (direction, cluster); ties broken alphabetically
  SELECT direction_id, cluster_id,
         (ARRAY_AGG(raw_name ORDER BY cnt DESC, raw_name))[1] AS mode_name
  FROM name_counts
  GROUP BY direction_id, cluster_id
)
SELECT
  v.direction_id,
  v.cluster_id,
  COUNT(*)                                              AS n_points,
  ST_Transform(ST_Centroid(ST_Collect(v.g3857)), 4326) AS centroid,
  -- normalized modal name: collapse internal whitespace + trim (fixes the
  -- stray-whitespace defects seen in hand-made layers)
  COALESCE(NULLIF(regexp_replace(BTRIM(nm.mode_name), '\s+', ' ', 'g'), ''),
           'Unnamed')                                   AS mode_name
FROM valid v
LEFT JOIN name_mode nm USING (direction_id, cluster_id)
GROUP BY v.direction_id, v.cluster_id, nm.mode_name;

CREATE INDEX IF NOT EXISTS _stop_clusters_centroid_gix
  ON transit._stop_clusters USING GIST (centroid);

/*----------------------------------------------------------------------------
 STEP 2: stops_auto  —  snap each cluster to its OWN-direction trip line,
                        classify terminal, then pair the two directions.
----------------------------------------------------------------------------*/
DROP MATERIALIZED VIEW IF EXISTS transit.stops_auto CASCADE;
CREATE MATERIALIZED VIEW transit.stops_auto AS
WITH
-- Re-attribute + cluster from raw so BOTH the confident clusters (Tier A) and
-- the DBSCAN noise points (Tier B rescue candidates) are available here.
attributed AS (
  SELECT
    s.gid                      AS src_id,
    NULLIF(BTRIM(s.name), '')  AS raw_name,   -- change to s.raw_name if needed
    t.direction_id             AS direction_id,
    ST_Transform(s.geom, 3857) AS g3857
  FROM raw.stops s
  JOIN raw.onboard_instances oi ON s.onboard_instance_observer_id = oi.id
  JOIN transit.trips         t  ON oi.trip_id = t.observer_id
  WHERE s.geom IS NOT NULL AND t.direction_id IS NOT NULL
),
clustered AS (
  SELECT a.*,
         ST_ClusterDBSCAN(g3857, eps := {dbscan_eps_m}, minpoints := {dbscan_minpoints})
           OVER (PARTITION BY direction_id) AS cid
  FROM attributed a
),
-- ---- Tier A: confident clusters (cid not null), one point each ----
tierA_names AS (
  SELECT direction_id, cid,
         (ARRAY_AGG(raw_name ORDER BY cnt DESC, raw_name))[1] AS mode_name
  FROM (
    SELECT direction_id, cid, raw_name, COUNT(*) AS cnt
    FROM clustered WHERE cid IS NOT NULL AND raw_name IS NOT NULL
    GROUP BY direction_id, cid, raw_name
  ) z
  GROUP BY direction_id, cid
),
tierA AS (
  SELECT cl.direction_id, cl.cid,
         COUNT(*) AS n_points,
         ST_Centroid(ST_Collect(cl.g3857)) AS c3857,
         COALESCE(NULLIF(regexp_replace(BTRIM(nm.mode_name), '\s+', ' ', 'g'), ''),
                  'Unnamed') AS mode_name
  FROM clustered cl
  LEFT JOIN tierA_names nm USING (direction_id, cid)
  WHERE cl.cid IS NOT NULL
  GROUP BY cl.direction_id, cl.cid, nm.mode_name
),
-- Snap each Tier-A cluster onto the nearest SAME-direction trip line.
tierA_snap AS (
  SELECT a.direction_id, a.n_points, a.mode_name,
         ST_ClosestPoint(ST_Transform(t.geom, 3857), a.c3857) AS g3857,
         ST_Distance(a.c3857, ST_Transform(t.geom, 3857))     AS dist_m,
         t.geom AS trip_geom
  FROM tierA a
  JOIN LATERAL (
    SELECT geom FROM transit.trips
    WHERE geom IS NOT NULL AND direction_id = a.direction_id
    ORDER BY a.c3857 <-> ST_Transform(geom, 3857) LIMIT 1
  ) t ON TRUE
),
tierA_typed AS (
  SELECT direction_id, n_points, mode_name, g3857,
         CASE WHEN ST_DWithin(g3857, ST_StartPoint(ST_Transform(trip_geom,3857)), {terminal_m})
                OR ST_DWithin(g3857, ST_EndPoint(ST_Transform(trip_geom,3857)),   {terminal_m})
              THEN 'Terminal' ELSE 'Informal' END AS stop_type
  FROM tierA_snap
  WHERE dist_m <= {snap_max_m}
),
-- ---- Tier B: COVERAGE RESCUE ----
-- Isolated on-route raw stops (DBSCAN noise) become stops when no confident
-- stop sits within rescue_gap, so sparsely-surveyed corridors are not left
-- with zero accessibility. Off-route noise is still dropped by snap_max.
noise_snap AS (
  SELECT n.direction_id,
         COALESCE(NULLIF(regexp_replace(BTRIM(n.raw_name), '\s+', ' ', 'g'), ''),
                  'Unnamed') AS mode_name,
         ST_ClosestPoint(ST_Transform(t.geom, 3857), n.g3857) AS g3857,
         ST_Distance(n.g3857, ST_Transform(t.geom, 3857))     AS dist_m
  FROM (SELECT direction_id, raw_name, g3857 FROM clustered WHERE cid IS NULL) n
  JOIN LATERAL (
    SELECT geom FROM transit.trips
    WHERE geom IS NOT NULL AND direction_id = n.direction_id
    ORDER BY n.g3857 <-> ST_Transform(geom, 3857) LIMIT 1
  ) t ON TRUE
),
-- Dedup rescue candidates among themselves so several isolated pings in one
-- gap collapse to a single rescued stop (minpoints:=1 → every point clusters).
noise_clustered AS (
  SELECT nk.*,
         ST_ClusterDBSCAN(g3857, eps := {rescue_gap_m}, minpoints := 1)
           OVER (PARTITION BY direction_id) AS ncid
  FROM noise_snap nk
  WHERE dist_m <= {snap_max_m}
),
noise_rep AS (
  SELECT direction_id,
         COUNT(*) AS n_points,
         ST_Centroid(ST_Collect(g3857)) AS g3857,
         (ARRAY_AGG(mode_name ORDER BY (mode_name = 'Unnamed'), mode_name))[1] AS mode_name
  FROM noise_clustered
  GROUP BY direction_id, ncid
),
rescued AS (
  SELECT r.direction_id, r.n_points, r.mode_name, r.g3857,
         'Informal'::text AS stop_type   -- rescued stops are mid-corridor by definition
  FROM noise_rep r
  WHERE NOT EXISTS (
    SELECT 1 FROM tierA_typed a
    WHERE a.direction_id = r.direction_id
      AND ST_DWithin(a.g3857, r.g3857, {rescue_gap_m})
  )
),
-- ---- Union Tier A + rescued ----
auto AS (
  SELECT
    ROW_NUMBER() OVER (ORDER BY direction_id, ST_X(g3857), ST_Y(g3857)) AS auto_id,
    direction_id, n_points, mode_name, stop_type,
    ST_Transform(g3857, 4326) AS geom,
    g3857
  FROM (
    SELECT direction_id, n_points, mode_name, stop_type, g3857 FROM tierA_typed
    UNION ALL
    SELECT direction_id, n_points, mode_name, stop_type, g3857 FROM rescued
  ) u
),
-- ---- Pairing: match a dir-0 stop to its nearest dir-1 stop within pair_max ----
-- One-directional-nearest with a used-once guarantee (DISTINCT ON on the dir-1
-- side, ordered by distance) approximates the greedy 1:1 matching used in the
-- validated prototype. Directions that diverge (one-way loop routing) fall
-- outside pair_max and remain independent single stops.
d0 AS (SELECT * FROM auto WHERE direction_id = 0),
d1 AS (SELECT * FROM auto WHERE direction_id = 1),
cand AS (
  SELECT d0.auto_id AS a0, m.auto_id AS a1,
         ST_Distance(d0.g3857, m.g3857) AS d
  FROM d0
  JOIN LATERAL (
    SELECT auto_id, g3857
    FROM d1
    WHERE ST_DWithin(d1.g3857, d0.g3857, {pair_max_m})
    ORDER BY d1.g3857 <-> d0.g3857
    LIMIT 1
  ) m ON TRUE
),
-- each dir-1 stop may only be claimed once (nearest dir-0 wins)
pairs AS (
  SELECT DISTINCT ON (a1) a0, a1, d
  FROM cand
  ORDER BY a1, d
),
grouped AS (
  SELECT a.*,
         CASE WHEN p.a0  IS NOT NULL THEN 'g' || p.a0::text END  AS grp_from_d0,
         CASE WHEN p2.a0 IS NOT NULL THEN 'g' || p2.a0::text END AS grp_from_d1
  FROM auto a
  LEFT JOIN pairs p  ON a.auto_id = p.a0     -- a is a dir-0 with a partner
  LEFT JOIN pairs p2 ON a.auto_id = p2.a1    -- a is a dir-1 claimed by a dir-0
)
SELECT
  auto_id,
  mode_name                                    AS stop_name,
  direction_id,
  0                                            AS location_type,   -- GTFS platform
  stop_type,
  n_points,
  COALESCE(grp_from_d0, grp_from_d1)           AS stop_group,
  geom
FROM grouped;

CREATE INDEX IF NOT EXISTS stops_auto_name_idx ON transit.stops_auto (stop_name);
CREATE INDEX IF NOT EXISTS stops_auto_geom_gix ON transit.stops_auto USING GIST (geom);
CREATE INDEX IF NOT EXISTS trips_geom_gix       ON transit.trips USING GIST (geom);

/*----------------------------------------------------------------------------
 STEP 3: transit.stops  —  final table, PRESERVING manual/synthetic stops.
   The table is created once (if absent). On every rebuild we delete ONLY the
   auto-derived rows and re-insert them, so hand-added synthetic stops
   (stop_source='manual') survive refreshes. stop_id mirrors gid.
----------------------------------------------------------------------------*/
CREATE TABLE IF NOT EXISTS transit.stops (
  gid           SERIAL PRIMARY KEY,
  stop_id       TEXT GENERATED ALWAYS AS (gid::text) STORED,
  stop_name     TEXT,
  stop_desc     TEXT,
  location_type INT,
  direction_id  INT,
  stop_group    TEXT,
  stop_type     TEXT,
  stop_source   TEXT DEFAULT 'auto',   -- 'auto' (rebuildable) | 'manual' (preserved)
  -- stop_lon / stop_lat intentionally omitted: geom is the single source of
  -- truth; derive coordinates with ST_X/ST_Y at read time.
  geom          geometry(Point, 4326)
);

CREATE UNIQUE INDEX IF NOT EXISTS transit_stops_stop_id_uidx ON transit.stops (stop_id);
CREATE INDEX IF NOT EXISTS transit_stops_geom_gist ON transit.stops USING GIST (geom);

-- Rebuild only the auto rows; manual/synthetic rows are left untouched.
DELETE FROM transit.stops WHERE stop_source = 'auto';

INSERT INTO transit.stops
  (stop_name, stop_desc, location_type, direction_id, stop_group, stop_type, stop_source, geom)
SELECT
  sa.stop_name,
  sa.stop_name AS stop_desc,
  sa.location_type,
  sa.direction_id,
  sa.stop_group,
  sa.stop_type,
  'auto',
  sa.geom
FROM transit.stops_auto sa;

-- Cross-direction name reconciliation for paired stops: give both members of a
-- stop_group the modal name of whichever side has more observations (and prefer
-- any real name over 'Unnamed').
WITH grp_pick AS (
  SELECT stop_group,
         (ARRAY_AGG(stop_name ORDER BY (stop_name = 'Unnamed'), n_pts DESC, stop_name))[1] AS best_name
  FROM (
    SELECT s.stop_group, s.stop_name, s.direction_id,
           COALESCE(sa.n_points, 0) AS n_pts
    FROM transit.stops s
    LEFT JOIN transit.stops_auto sa
      ON sa.stop_group = s.stop_group AND sa.direction_id = s.direction_id
    WHERE s.stop_group IS NOT NULL AND s.stop_source = 'auto'
  ) q
  GROUP BY stop_group
)
UPDATE transit.stops s
SET stop_name = gp.best_name,
    stop_desc = gp.best_name
FROM grp_pick gp
WHERE s.stop_group = gp.stop_group
  AND s.stop_source = 'auto';
