DROP MATERIALIZED VIEW IF EXISTS transit.trip_stops_sequence;
CREATE MATERIALIZED VIEW transit.trip_stops_sequence AS(

with recursive stop_distance_along_trip as 
(SELECT 
    row_number() OVER () AS gid,
    t.gid as t_id,
    t.observer_id  as observer_trip_id,
    s.gid::text AS stop_id,
--    s.stop_id as gtfs_id,
    s.stop_name as stop_name,
    st_linelocatepoint(t.geom, s.geom::geometry) * st_length(t.geom::geography) AS distance,
    st_linelocatepoint(t.geom, s.geom::geometry) AS distance_frac,
    v.name::text AS vehicle_name
  FROM   
    	transit.stops s
    	join transit.trips t on st_dwithin(t.geom::geography, s.geom::geography, 1::real)
	    LEFT JOIN transit.agencies t1 ON t.agency_id = t1.gid
	    LEFT JOIN transit.vehicles v ON t1.vehicle_id = v.gid
),
enriched_pairs as (
  select *, 
  distance - lag(distance, 1) over (order by observer_trip_id, distance_frac) as distance_from_prev,
  row_number() OVER (PARTITION BY t_id ORDER BY distance_frac) AS stop_sequence

  from stop_distance_along_trip
)
select * from enriched_pairs where distance_from_prev >= 100 

);
--------------------------------------
CREATE INDEX IF NOT EXISTS trip_stops_sequence_idx
  ON transit.trip_stops_sequence
  (t_id);