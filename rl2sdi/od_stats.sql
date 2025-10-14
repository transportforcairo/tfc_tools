DROP MATERIALIZED VIEW IF EXISTS transit.od_stats;
CREATE MATERIALIZED VIEW transit.od_stats AS (
    WITH RECURSIVE
    od_segments AS(
        SELECT DISTINCT ON
        (   o_id,
            d_id,
            vehicle_name,
            geom)
            o_id,
            d_id,
            vehicle_name,
            avg(dist) as dist,
            geom,
            count(*)
        FROM (
                SELECT
                    *
                FROM
                    (
                        SELECT
                            ts.stop_id as o_id,
                            LEAD(ts.stop_id, 1) OVER ( PARTITION BY ts.t_id ORDER BY ts.t_id, ts.stop_sequence ) as d_id,
                            ts.vehicle_name,
                            LEAD(ts.distance, 1) OVER ( PARTITION BY ts.t_id ORDER BY ts.t_id, ts.stop_sequence ) - ts.distance as dist, 
                            st_linesubstring( trip.geom, ts.distance_frac, LEAD(ts.distance_frac, 1) OVER ( PARTITION BY ts.t_id ORDER BY ts.t_id, ts.stop_sequence )) as geom
                        FROM
                            transit.trip_stops_sequence ts JOIN transit.trips_view trip ON ts.t_id = trip.gid AND ts.distance_frac < 1
                    ) AS _
                WHERE
                    d_id IS NOT NULL
            ) AS _
        GROUP BY
            o_id,
            d_id,
            vehicle_name,
            geom
    ),
    trackpoints_near_stops AS (
        SELECT 
            ss.gid as stop_id,
			ss.gid as gtfs_id,
            onboard_instance_id as instance_id,
            timestamp as time
        FROM
            raw.trackpoints t JOIN (select * from raw.onboard_instances raw_oi where raw_oi.valid = 'true' and raw_oi.status = 'finished') oi ON t.onboard_instance_id = oi.id,
            transit.stops ss
        WHERE
            ST_DWithin(t.geom::geography, ss.geom::geography, 30)
    ),
    o_timestamps AS(
        SELECT 
            seg.o_id as o_id,
			min(trackpoints_near_stops.gtfs_id) as o_gtfs_id, -- Aggregate fn. used only to not need to group redundantly by gtfs_id, since o_id and gtfs_id are both unique per stop.
            instance_id as o_instance_id, 
            max(time) as o_time,
            v."name" as o_vehicle_name
        FROM trackpoints_near_stops, od_segments as seg, raw.onboard_instances oi , transit.vehicles v 
        -- WHERE trackpoints_near_stops.stop_id::text = seg.o_id and oi.id  = instance_id and oi.vehicle_id = v.gid 
        WHERE trackpoints_near_stops.stop_id::text = seg.o_id and oi.id  = instance_id and oi.vehicle_id::text = v.gid::text
        GROUP BY instance_id, seg.o_id, v."name" 
    ),
    d_timestamps AS(
        SELECT 
            seg.d_id as d_id,
			min(trackpoints_near_stops.gtfs_id) as d_gtfs_id,
            instance_id as d_instance_id, 
            max(time) as d_time,
            v."name" as d_vehicle_name
        FROM trackpoints_near_stops, od_segments as seg, raw.onboard_instances oi , transit.vehicles v 
        -- WHERE trackpoints_near_stops.stop_id::text = seg.d_id and oi.id  = instance_id and oi.vehicle_id = v.gid 
        WHERE trackpoints_near_stops.stop_id::text = seg.d_id and oi.id  = instance_id and oi.vehicle_id::text = v.gid::text
        GROUP BY instance_id, seg.d_id, v."name" 
    ),
    od_timestamps AS (
        SELECT 
            o_id,
            d_id,
			o_gtfs_id as from_id,
            d_gtfs_id as to_id,
            o_time,
            d_time,
            o_instance_id as instance_id,
            o_vehicle_name as vehicle_name
        FROM o_timestamps JOIN d_timestamps ON o_instance_id = d_instance_id
    ),
    avg_durations AS (
        SELECT
            row_number() over () as gid,
            o_id,
            d_id,
            i.gid interval_id,
            vehicle_name,
            min(tmp.from_id) as from_id,
            min(tmp.to_id) as to_id,
            min(i.start_time) interval_start,
            avg(time_diff)::int duration
        FROM
        
            ( SELECT *,  EXTRACT (EPOCH FROM (d_time - o_time)) as time_diff FROM od_timestamps ) AS tmp,
            transit.intervals i
        WHERE
            (time_diff > 0) AND (o_time::time BETWEEN i.start_time and i.end_time) AND (i.active)
        GROUP BY
            o_id,
            d_id,
            i.gid,
            vehicle_name
    )
    SELECT
        row_number() over () as gid,
		o_id,
        d_id,
        interval_id,
        interval_start,
        vehicle_name,
        od_segments.dist as dist,
        duration,
        (od_segments.dist/duration)*3.6 as speed,	-- Convert m/s to km/h (1 m/s = 3.6 km/h)
        od_segments.geom::geometry(LINESTRING, 4326)
    FROM        
        od_segments JOIN avg_durations USING (o_id, d_id, vehicle_name)
        ); -- WE ADDED THE SEMICOLON TO WORK WITH PSYCOPG2 FUNCTION run_ddl_sql IN script4plugin.py

