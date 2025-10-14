SELECT
        t0.gid gid, 
        t0.observer_route_id route_id, 
        t0.route_type route_type,
        t0.service_id service_id,
        CONCAT(t1.common_name, ' ', t0.agency_serial) route_short,
        
        CASE
            WHEN LOWER(origin_terminal.name)>LOWER(dest_terminal.name) THEN CONCAT(dest_terminal.name,' - ',origin_terminal.name)
            ELSE CONCAT(origin_terminal.name,' - ',dest_terminal.name)
        END route_long,
        
        t0.observer_id  observer_id, 
        t0.direction_id direction_id, 
        t0.o_id o_id, 
        t0.d_id d_id, 
        t0.geom geom,
        origin_terminal.name origin, 
        dest_terminal.name destination,
        t1.agency_id agency_id, 
        v.name vehicle_name, 
        v.passenger_capacity passenger_capacity,
        trip_short, 
        fare,
        len_m/1000 len_km
    FROM 
        transit.trips as t0
        LEFT JOIN transit.agencies t1 ON t0.agency_id = t1.gid
        LEFT JOIN transit.vehicles v ON t1.vehicle_id = v.gid
        LEFT JOIN transit.terminals origin_terminal ON t0.o_id = origin_terminal.gid
        LEFT JOIN transit.terminals dest_terminal ON t0.d_id = dest_terminal.gid
        ,
        ST_LENGTH(t0.geom::geography) len_m
        ,
        CONCAT(origin_terminal.name, '-', dest_terminal.name) trip_short;