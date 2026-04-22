DROP MATERIALIZED VIEW IF EXISTS transit.od_stats CASCADE;
CREATE MATERIALIZED VIEW transit.od_stats AS (
    WITH RECURSIVE
    -- ========================================================================
    -- PHASE 1: Build od_segments + compute Tier-1 (calculated) durations.
    -- Identical to the previous (non-tiered) version up through od_observations.
    -- ========================================================================
    od_segments AS (
        SELECT DISTINCT ON (o_id, d_id, vehicle_name, geom)
            o_id,
            d_id,
            vehicle_name,
            avg(dist) as dist,
            geom,
            count(*)
        FROM (
            SELECT * FROM (
                SELECT
                    ts.stop_id as o_id,
                    LEAD(ts.stop_id, 1) OVER (PARTITION BY ts.t_id ORDER BY ts.t_id, ts.stop_sequence) as d_id,
                    ts.vehicle_name,
                    ts.distance_frac as o_frac,
                    LEAD(ts.distance, 1) OVER (PARTITION BY ts.t_id ORDER BY ts.t_id, ts.stop_sequence) - ts.distance as dist,
                    st_linesubstring(trip.geom, ts.distance_frac,
                        LEAD(ts.distance_frac, 1) OVER (PARTITION BY ts.t_id ORDER BY ts.t_id, ts.stop_sequence)) as geom
                FROM transit.trip_stops_sequence ts
                JOIN transit.trips_view trip ON ts.t_id = trip.gid
                -- Note: distance_frac<1 must filter only the ORIGIN row, not the LEAD target.
                -- Previously this was on the JOIN condition which also dropped rows that
                -- would have been the DESTINATION of the preceding segment — making the
                -- last segment of every trip vanish (whenever its destination was the
                -- trip terminal at frac=1.0). Mirrors Python's `if af >= 1.0: continue`.
            ) AS _
            WHERE d_id IS NOT NULL AND o_frac < 1
        ) AS _
        GROUP BY o_id, d_id, vehicle_name, geom
    ),
    trackpoints_near_stops AS (
        SELECT
            ss.gid as stop_id,
            ss.gid as gtfs_id,
            onboard_instance_id as instance_id,
            timestamp as time
        FROM
            raw.trackpoints t
            JOIN (
                SELECT * FROM raw.onboard_instances raw_oi
                WHERE raw_oi.valid = 'true' AND raw_oi.status = 'finished'
            ) oi ON t.onboard_instance_id = oi.id,
            transit.stops ss
        WHERE ST_DWithin(t.geom::geography, ss.geom::geography, {trackpoint_buffer_m})
    ),
    o_timestamps AS (
        -- vehicle_name is taken from od_segments (i.e. trip_stops_sequence's
        -- declared vehicle type for each OD pair), not from the onboard
        -- instance's registered vehicle. This mirrors the Python pipeline —
        -- see refresh_sdi_derived_algorithm.py::_build_od_stats where
        -- tp_near_df is joined to od_segments[["o_id","vehicle_name"]].
        -- Historically this CTE joined raw.onboard_instances → transit.vehicles
        -- and used v."name" as o_vehicle_name; that produced tier-1 rows under
        -- vehicle_name combos Python would reject (because the onboard
        -- instance's registered vehicle can differ from the vehicle_name
        -- declared on the trip).
        SELECT
            seg.o_id as o_id,
            min(trackpoints_near_stops.gtfs_id) as o_gtfs_id,
            instance_id as o_instance_id,
            max(time) as o_time,
            seg.vehicle_name as o_vehicle_name
        FROM trackpoints_near_stops, od_segments as seg, raw.onboard_instances oi
        WHERE trackpoints_near_stops.stop_id::text = seg.o_id
          AND oi.id = instance_id
        GROUP BY instance_id, seg.o_id, seg.vehicle_name
    ),
    d_timestamps AS (
        SELECT
            seg.d_id as d_id,
            min(trackpoints_near_stops.gtfs_id) as d_gtfs_id,
            instance_id as d_instance_id,
            max(time) as d_time,
            seg.vehicle_name as d_vehicle_name
        FROM trackpoints_near_stops, od_segments as seg, raw.onboard_instances oi
        WHERE trackpoints_near_stops.stop_id::text = seg.d_id
          AND oi.id = instance_id
        GROUP BY instance_id, seg.d_id, seg.vehicle_name
    ),
    od_timestamps AS (
        -- Python also requires the vehicle_name to match between the o and d
        -- sides of the same instance (see _build_od_stats: merge on
        -- ["instance_id", "vehicle_name"]). Enforce that here too.
        SELECT
            o_id, d_id,
            o_gtfs_id as from_id,
            d_gtfs_id as to_id,
            o_time, d_time,
            o_instance_id as instance_id,
            o_vehicle_name as vehicle_name
        FROM o_timestamps
        JOIN d_timestamps ON o_instance_id = d_instance_id
                          AND o_vehicle_name = d_vehicle_name
    ),
    od_observations AS (
        SELECT
            tmp.o_id, tmp.d_id, tmp.vehicle_name, tmp.from_id, tmp.to_id,
            tmp.time_diff, i.gid AS interval_id, i.name AS interval_name,
            i.start_time AS interval_start
        FROM (
            SELECT *, EXTRACT(EPOCH FROM (d_time - o_time)) AS time_diff
            FROM od_timestamps
        ) AS tmp,
            transit.intervals i
        WHERE (tmp.time_diff > 0)
          AND (tmp.o_time::time BETWEEN i.start_time AND i.end_time)
          AND (i.active)
    ),
    -- Tier 1 — CALCULATED durations. GROUP BY uses {vehicle_group_expr}:
    --   distinguish on  → obs.vehicle_name      (one average per vehicle type)
    --   distinguish off → '_pooled_'::text      (pooled across vehicle types)
    -- vehicle_name is projected via MIN() so the column exists in both modes;
    -- downstream joins use {{vehicle_join_condition_*}} which is empty in pool mode.
    tier1_calculated AS (
        SELECT
            obs.o_id,
            obs.d_id,
            obs.interval_id,
            min(obs.interval_name) as interval_name,
            min(obs.vehicle_name) as vehicle_name,
            min(obs.from_id) as from_id,
            min(obs.to_id) as to_id,
            min(obs.interval_start) as interval_start,
            avg(obs.time_diff)::int as duration,
            count(*) as n_samples,
            'calculated'::text as calc_method
        FROM od_observations obs
        GROUP BY obs.o_id, obs.d_id, obs.interval_id, {vehicle_group_expr}
    ),
    -- ========================================================================
    -- PHASE 2: Build the full candidate set (od_segments × intervals) and layer
    -- on tiered estimations. All estimation tiers (2-6) draw from
    -- tier1_calculated ONLY (no cascade), so behavior does not depend on CTE
    -- evaluation order.
    -- ========================================================================
    intervals_active AS (
        SELECT gid AS interval_id,
               name AS interval_name,
               to_char(start_time, 'HH24:MI:SS') AS interval_start
        FROM transit.intervals WHERE active = true
    ),
    -- Every candidate OD row: one per (segment × interval).
    candidates AS (
        SELECT seg.o_id, seg.d_id, seg.vehicle_name, seg.dist, seg.geom,
               iv.interval_id, iv.interval_name, iv.interval_start
        FROM od_segments seg CROSS JOIN intervals_active iv
    ),
    -- Tier 1 speed per (o_id, d_id, vehicle_name, interval_id). Used by tiers 2-6.
    -- In pool mode, tier1_calculated.vehicle_name is '_pooled_' -sentinel is dropped
    -- from the join condition, so one pooled row matches every vehicle's segment.
    --
    -- NOTE: `od_segments` may contain multiple rows for the same (o_id,d_id,vehicle_name)
    -- if different trips produced slightly-different geometric substrings of that
    -- segment (same is true in the Python side). Since the per-segment `dist` is
    -- already averaged inside `od_segments`, the distinct geom rows carry the same
    -- dist value. We pre-aggregate to one row per key here so downstream joins on
    -- `base` don't multiply candidate rows.
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
    -- Apply Tier 1 to candidates; null duration/calc_method for rows without a match.
    base AS (
        SELECT c.o_id, c.d_id, c.vehicle_name, c.dist, c.geom,
               c.interval_id, c.interval_name, c.interval_start,
               t1.from_id, t1.to_id,
               t1.duration AS t1_duration,
               t1.n_samples AS t1_n_samples,
               t1.calc_method AS t1_calc_method,
               t1w.speed_mps AS t1_speed_mps
        FROM candidates c
        LEFT JOIN tier1_calculated t1
          ON c.o_id = t1.o_id AND c.d_id = t1.d_id AND c.interval_id = t1.interval_id
          {vehicle_join_condition_c_t1}
        LEFT JOIN tier1_with_speed t1w
          ON c.o_id = t1w.o_id AND c.d_id = t1w.d_id AND c.interval_id = t1w.interval_id
          AND c.vehicle_name = t1w.vehicle_name
    ),
    -- ---- Trip-sequence context for Tier 2 ----
    -- Each (o_id, d_id, vehicle_name) segment may belong to multiple trips.
    -- To mirror Python _fill_tier2_trip_neighbors — which walks each trip
    -- independently and uses setdefault() (first-writer-wins across trips) —
    -- we must expose a candidate row on EVERY trip that contains its segment,
    -- not a single representative trip. A previous version used a
    -- seg_trip_one CTE picking the lowest t_id; that missed cases where a
    -- candidate's segment and its tier-1 neighbor are on the same trip but
    -- the neighbor's "representative" trip was a different one.
    tss_segments AS (
        SELECT t_id,
               stop_id::text AS o_id,
               LEAD(stop_id, 1) OVER (PARTITION BY t_id ORDER BY stop_sequence)::text AS d_id,
               vehicle_name,
               stop_sequence
        FROM transit.trip_stops_sequence
    ),
    -- Expand base: one row per (candidate × trip containing its segment).
    -- For rows whose segment appears on no trip (unlikely, but possible),
    -- the LEFT JOIN keeps them with NULL t_id so the tier-3/4/5/6 fallbacks
    -- can still fire downstream.
    with_trip_ctx_all AS (
        SELECT b.*, st.t_id, st.stop_sequence
        FROM base b
        LEFT JOIN tss_segments st
          ON b.o_id = st.o_id AND b.d_id = st.d_id AND b.vehicle_name = st.vehicle_name
         AND st.d_id IS NOT NULL
    ),
    -- ---- Tier 2: running-fill across same trip + interval + vehicle ----
    -- Use the array_agg(FILTER)[idx] trick as a PG 12-compatible substitute for
    -- LAST_VALUE / FIRST_VALUE ... IGNORE NULLS (which is PG 16+).
    --
    -- prev_speed: the last non-null tier-1 speed at or before the current row,
    --             within the same (t_id, interval_id, vehicle_name) partition.
    -- next_speed: the first non-null tier-1 speed at or after the current row.
    -- prev_seq/next_seq: the stop_sequence of those neighbors, used for the
    --                   inverse-distance weighting.
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
    -- Per-trip tier-2 candidate: inverse-sequence-distance-weighted avg of
    -- neighbor speeds (or single-sided where only one side has a match).
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
        WHERE t1_duration IS NULL           -- only fill rows actually missing tier 1
          AND t_id IS NOT NULL              -- only rows that landed on a trip
    ),
    -- Collapse to one tier-2 speed per (o,d,veh,iv). Python uses setdefault()
    -- → first-writer-wins; here that's DISTINCT ON ordered by t_id ASC, which
    -- is deterministic and consistent with the lowest-trip-id winning.
    t2_speed AS (
        SELECT DISTINCT ON (o_id, d_id, vehicle_name, interval_id)
               o_id, d_id, vehicle_name, interval_id, t2_speed
        FROM t2_per_trip
        WHERE t2_speed IS NOT NULL
        ORDER BY o_id, d_id, vehicle_name, interval_id, t_id
    ),
    -- ---- Tier 3: same segment (o,d[,vehicle]), averaged across other intervals ----
    -- The grouping keys depend on distinguish_speeds_by_vehicle. The
    -- {tier3_group_keys} placeholder expands to "o_id, d_id, vehicle_name"
    -- in distinguish mode, or "o_id, d_id" in pool mode. The corresponding
    -- join condition {tier3_join_condition_b_t3} matches.
    tier3_speed AS (
        SELECT {tier3_group_keys},
               AVG(speed_mps) AS t3_speed
        FROM tier1_with_speed
        GROUP BY {tier3_group_keys}
    ),
    -- ---- Tier 4: same trip + vehicle + interval avg speed ----
    -- Python walks each trip this segment belongs to (seg_to_trips), fetches
    -- the per-(trip, vehicle, interval) avg speed for each, and averages
    -- those. We mirror by (a) exploding each tier-1 segment across all trips
    -- it appears on, (b) computing per-(trip, vehicle, interval) avg,
    -- (c) re-aggregating back per-(segment, interval) by averaging across
    -- the trips this segment belongs to.
    t1_trip_link AS (
        SELECT t1w.o_id, t1w.d_id, t1w.vehicle_name, t1w.interval_id,
               t1w.speed_mps, st.t_id
        FROM tier1_with_speed t1w
        JOIN tss_segments st
          ON t1w.o_id = st.o_id AND t1w.d_id = st.d_id AND t1w.vehicle_name = st.vehicle_name
        WHERE st.d_id IS NOT NULL
    ),
    tier4_trip_speed AS (
        -- per-(trip, vehicle, interval) average speed
        SELECT t_id, vehicle_name, interval_id, AVG(speed_mps) AS trip_speed
        FROM t1_trip_link
        GROUP BY t_id, vehicle_name, interval_id
    ),
    seg_all_trips AS (
        -- each segment → list of trips it appears on
        SELECT DISTINCT o_id, d_id, vehicle_name, t_id
        FROM tss_segments
        WHERE d_id IS NOT NULL
    ),
    tier4_speed AS (
        -- average the per-trip speeds across all trips this segment belongs to
        SELECT sat.o_id, sat.d_id, sat.vehicle_name, tts.interval_id,
               AVG(tts.trip_speed) AS t4_speed
        FROM seg_all_trips sat
        JOIN tier4_trip_speed tts
          ON sat.t_id = tts.t_id AND sat.vehicle_name = tts.vehicle_name
        GROUP BY sat.o_id, sat.d_id, sat.vehicle_name, tts.interval_id
    ),
    -- ---- Tier 5: same vehicle + interval citywide ----
    tier5_speed AS (
        SELECT vehicle_name, interval_id, AVG(speed_mps) AS t5_speed
        FROM tier1_with_speed
        GROUP BY vehicle_name, interval_id
    ),
    -- ---- Tier 6: per-interval global, with overall-global fallback ----
    tier6_per_interval AS (
        SELECT interval_id, AVG(speed_mps) AS t6_speed_iv
        FROM tier1_with_speed GROUP BY interval_id
    ),
    tier6_global AS (
        SELECT AVG(speed_mps) AS t6_speed_all FROM tier1_with_speed
    ),
    -- Gather candidate speeds from every tier (each is nullable).
    -- We read from `base` (1 row per candidate) and left-join each tier's
    -- per-key speed. Previously assembled read from `with_neighbors`, but
    -- that CTE is now trip-exploded (see with_trip_ctx_all comments).
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
    -- Final: pick the duration and calc_method by tier priority.
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
        row_number() over () as gid,
        f.o_id, f.d_id,
        f.interval_id, f.interval_name, f.interval_start,
        f.vehicle_name,
        f.dist,
        f.duration_final as duration,
        (f.dist / NULLIF(f.duration_final, 0)) * 3.6 as speed,
        f.calc_method_final as calc_method,
        COALESCE(f.t1_n_samples, 0) as n_samples,
        f.geom::geometry(LINESTRING, 4326)
    FROM final f
    WHERE f.duration_final IS NOT NULL
); -- semicolon required by run_ddl_sql in script4plugin.py
