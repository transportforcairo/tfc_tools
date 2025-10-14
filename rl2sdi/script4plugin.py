def run_migration(observer_db_params, sdi_db_params, observer_project_id, feedback):
    # %%
    import geopandas as gpd
    import pandas as pd
    import psycopg2
    from sqlalchemy import create_engine, text
    from sqlalchemy import Table, Column, Integer, String, MetaData, ForeignKey, Sequence, DateTime
    from sqlalchemy import inspect
    from sqlalchemy.schema import CreateSchema # TRIAL TO CREATE SCHEMAS DIFFERENTLY
    from geoalchemy2 import Geometry, WKTElement
    import matplotlib.pyplot as plt
    import numpy as np
    from shapely import wkt, wkb, geometry
    import urllib.parse

    import os
    def sql_path(filename):
        return os.path.join(os.path.dirname(__file__), filename)

    # %%
    import sqlalchemy
    sqlalchemy.__version__

    # %%
    OBSERVER_PROJECT_ID = observer_project_id #MODIFIED THIS

    # %%
    def get_pg_connection_url(dbname, user, password, host, port):
        safe_password = urllib.parse.quote_plus(password)
        return f"postgresql+psycopg2://{user}:{safe_password}@{host}:{port}/{dbname}"

    observer_db_engine = create_engine(
        get_pg_connection_url(**observer_db_params), pool_recycle=3600 #MODIFIED VARIABLE NAME HERE
    )
    observer_db_con = observer_db_engine.connect()

    feedback.pushInfo("RouteLab database is connected")

    # %%
    sdi_db_engine = create_engine(
        get_pg_connection_url(**sdi_db_params), pool_recycle=3600 #MODIFIED VARIABLE NAME HERE
    )
    city_db_con = sdi_db_engine.connect()

    feedback.pushInfo("City database is connected")

    # %%
    # NEW FUNCTION FOR EXECUTING SQL WITH PSYCOPG2 AND NOT SQLALCHEMY
    def run_ddl_sql(sql, connection):
        try:
            raw_connection = connection.connection
            cursor = raw_connection.cursor()
            cursor.execute(sql)
            raw_connection.commit()
            cursor.close()
        except Exception as e:
            raise RuntimeError(f"Failed to execute DDL: {e}")

    # %%
    trips = gpd.read_postgis(
        con=observer_db_con,
        sql=f"""
        select
            *
        from
            v_trips_ext vte
        where
            project_id = '{OBSERVER_PROJECT_ID}'
            and deleted_at is null
            and geometry is not null
        """,
        geom_col="geometry",
    )

    feedback.pushInfo("City SQL is finished")

    # %%
    observer_setting_id = observer_db_con.execute(
        text(f"""select id from settings where project_id='{OBSERVER_PROJECT_ID}' and deleted_at is null""")
    ).first()[0]

    # %%
    # observer_setting_id

    feedback.pushInfo("observer_setting_id step is finished")

    # %%

    run_ddl_sql("CREATE EXTENSION IF NOT EXISTS postgis;", city_db_con)

    feedback.pushInfo("CREATE EXTENSION IF NOT EXISTS postgis is finished")

    # %%

    ''' trial 2 using psycopg - THIS ONE WORKED!! '''
    # import psycopg2

    def create_schemas_directly(db_params, feedback):
        try:
            conn = psycopg2.connect(
                host=db_params["host"],
                port=db_params["port"],
                dbname=db_params["dbname"],
                user=db_params["user"],
                password=db_params["password"]
            )
            cur = conn.cursor()
            cur.execute("CREATE SCHEMA IF NOT EXISTS transit;")
            cur.execute("CREATE SCHEMA IF NOT EXISTS raw;")
            conn.commit()
            cur.close()
            conn.close()
            feedback.pushInfo("✅ Schemas 'transit' and 'raw' created successfully via psycopg2.")
        except Exception as e:
            feedback.reportError(f"⚠️ Failed to create schemas via psycopg2: {e}")

    create_schemas_directly(sdi_db_params, feedback)

    feedback.pushInfo("creating transit and raw schemas is finished")

    # NEW LINE FROM GITHUB UPDATE

    # I added CASCADE here to avoid an error of dependencies. 
    # NOTE: CASCDE was not part of the original script logic
    run_ddl_sql("""
        drop table if exists transit.vehicles CASCADE; 
        create table if not exists transit.vehicles (
            gid   SERIAL primary key,
            name VARCHAR,
            passenger_capacity  INT
        );
    """, city_db_con)

    # %%
    # NEW FROM GITHUB UPDATE
    # Observer schema doesn't have any of the fields
    
    # Create agencies table
    agencies = pd.read_sql(
        con=observer_db_con,
        sql=text(f"select * from agencies where setting_id='{observer_setting_id}' and deleted_at is null"),
    )

    # there was a part specific to [69] about types of agencies that I did not add here

    run_ddl_sql("""
        DROP TABLE if exists transit.agencies;
        create table if not exists transit.agencies (
            gid SERIAL primary key,
            agency_id text,
            agency_name text,
            agency_url text,
            agency_timezone text,
            common_name text,
            has_serial bool,
            vehicle_id int4,
            CONSTRAINT agencies_fk FOREIGN KEY (vehicle_id) REFERENCES transit.vehicles(gid)
        );
        """, city_db_con)
    
    vehicle_names = (
    agencies["name"]
    .dropna()
    .drop_duplicates()
    .tolist()
    )

    # Insert vehicles and capture gid<->name mapping
    insert_vehicle_sql = sqlalchemy.text("""
    INSERT INTO transit.vehicles (name, passenger_capacity)
    VALUES (:name, :capacity)
    RETURNING gid, name
    """)

    vehicle_gid_by_name = {}
    with city_db_con.begin():
        for nm in vehicle_names:
            row = city_db_con.execute(insert_vehicle_sql, {"name": nm, "capacity": None}).mappings().one()
            vehicle_gid_by_name[row["name"]] = row["gid"]

    # 2) Prepare payload for transit.agencies
    mapped = agencies.copy()
    mapped["has_serial"] = mapped["serial"].notna()

    payload = []
    for _, r in (
        mapped[["id", "name", "has_serial"]]
        .drop_duplicates(subset=["id"])
        .iterrows()
    ):
        payload.append({
            "agency_id":        r["name"],
            "agency_name":      r["name"],
            "agency_url":       None,
            "agency_timezone":  None,
            "common_name":      r["name"],
            "has_serial":       bool(r["has_serial"]),
            "vehicle_id":       vehicle_gid_by_name.get(r["name"])  # link by name
        })

    # 3) Insert agencies with FK to vehicles
    insert_agencies_sql = sqlalchemy.text("""
        INSERT INTO transit.agencies
            (agency_id, agency_name, agency_url, agency_timezone, common_name, has_serial, vehicle_id)
        VALUES
            (:agency_id, :agency_name, :agency_url, :agency_timezone, :common_name, :has_serial, :vehicle_id)
    """)

    city_db_con.execute(insert_agencies_sql, payload)


    feedback.pushInfo("transit.agencies is finished")
    feedback.pushInfo("creating vehicles is finished")



    # %%
    # NEW FROM GITHUB UPDATE
    pd.read_sql(con=city_db_con, sql=text("select * from transit.agencies"))

    # %%
    '''THE NEXT PART IS FULLY FROM THE NEW GITHUB'''
    
    ## Terminals
    terminals = gpd.read_postgis(
        con=observer_db_con,
        sql=f"select * from terminals where setting_id='{observer_setting_id}' and deleted_at is null",
        geom_col="geometry",
    )

    # %%
    # Transform terminals
    valid_terminals = (
        terminals.copy()
        .query("status == 'accepted'")[["id", "name", "geometry"]]
        .rename({"geometry": "geom", "id": "observer_id"}, axis=1)
        .dropna()
        .set_geometry("geom")
    )
    valid_terminals = valid_terminals[~valid_terminals.name.str.lower().str.contains("test")]

    # %%
    # )  # NOTE: Always wrap raw SQL strings passed to .execute() with text(...) to deal with multi-line strings when working with sqlalchemy
   
    run_ddl_sql("DROP table if exists transit.terminals CASCADE;", city_db_con)
   
    run_ddl_sql("""
        create table if not exists transit.terminals (
            gid SERIAL primary key,
            name text,
            name_ar text,
            geom GEOMETRY(GEOMETRY, 4326),
            observer_id text
        );
        """, city_db_con)
    valid_terminals.to_postgis(
        name="terminals",
        schema="transit",
        con=city_db_con,
        if_exists="append",
        dtype={"geom": Geometry("GEOMETRY", srid=4326)},
    )

    feedback.pushInfo("transit.terminals is finished")

    # %%
    # Reading the intervals from Observer
    project_intervals = pd.read_sql(con=observer_db_con, sql=f"""
    select i.id, "start" ,"end", i.days, name, i.deleted_at 
    from 
        intervals i 
            join frequency_setting_intervals fsi on i.id = fsi.interval_id 
            join frequency_settings fs on fsi.frequency_setting_id = fs.id
    where 
        fs.setting_id='{observer_setting_id}'
    """)

    run_ddl_sql("""
        drop table if exists transit.intervals CASCADE;
        create table if not exists transit.intervals (
            gid serial primary key,
            start_time time,
            end_time time,
            observer_id text,
            name varchar,
            active boolean
        );
        """, city_db_con)

    project_intervals["active"] = project_intervals["deleted_at"].isna()

    payload = []
    for _, r in (
        project_intervals[["id", "name", "start", "end", "active"]]
        .iterrows()
    ):
        payload.append({
            "start_time":        r["start"],
            "end_time":      r["end"],
            "observer_id":       r["id"],
            "name":  r["name"],
            "active":      bool(r["active"])
        })

    # 3) Insert agencies with FK to vehicles
    insert_intervals_sql = sqlalchemy.text("""
        INSERT INTO transit.intervals
            (start_time, end_time, observer_id, name, active)
        VALUES
            (:start_time, :end_time, :observer_id, :name, :active)
    """)

    city_db_con.execute(insert_intervals_sql, payload)
    
    feedback.pushInfo("transit.intervals is finished")

    # %%
    # this part was relocated above as per the new version of the originaln script.ipynb
    sdi_intervals = pd.read_sql(
        con=city_db_con, sql="""select * from transit.intervals"""
    )

    sdi_terminals = gpd.read_postgis(
        con=city_db_con, sql="""select * from transit.terminals"""
    )
    sdi_agencies = pd.read_sql(
        con=city_db_con, sql="""select * from transit.agencies"""
    )

    observer_trips = gpd.read_postgis(
        con=observer_db_con,
        sql=f"""select *, ST_POINT(0,0) from v_trips_ext where setting_id='{observer_setting_id}' and deleted_at is null""",
        geom_col="geometry",
    )

    observer_trips = (
        observer_trips.merge(
            sdi_terminals[["observer_id", "gid"]],
            left_on="origin_id",
            right_on="observer_id",
        )
        .rename({"gid": "o_id"}, axis=1)
        .drop(["observer_id"], axis=1)
        .merge(
            sdi_terminals[["observer_id", "gid"]],
            left_on="destination_id",
            right_on="observer_id",
        )
        .rename({"gid": "d_id"}, axis=1)
        .drop(["observer_id"], axis=1)
    )

    observer_trips = (
        observer_trips.drop(["agency_id"], axis=1)
        .merge(
            sdi_agencies[["agency_id", "gid"]].rename(columns={"agency_id": "agency"}), left_on="agency", right_on="agency"
        )
        .rename({"gid": "agency_id"}, axis=1)
    )

    # %%
    '''the large section below was relocated from down to up here'''

    observer_trips = observer_trips.rename(
        {
            "direction": "direction_id",
            "geometry": "geom",
            "id": "observer_id",
            "route_id": "observer_route_id",
            "bus_number": "agency_serial",
        },
        axis=1,
    ).set_geometry("geom")
    observer_trips = observer_trips[
        [
            "o_id",
            "d_id",
            "agency_id",
            "geom",
            "direction_id",
            "observer_id",
            "observer_route_id",
            "agency_serial",
            "fare",
        ]
    ]

    observer_trips["fare"] = pd.to_numeric(observer_trips["fare"], errors="ignore")
    print(observer_trips['fare'].isna().sum())
    observer_trips.fare = observer_trips.fare.fillna(np.ceil(observer_trips.fare.mean()))

    observer_trips["service_id"] = "Ground_Daily"

    # %%
    frequency_instances = pd.read_sql(
        con=observer_db_con,
        sql=f"select * from v_frequency_instances_ext where setting_id='{observer_setting_id}' and deleted_at is null",
    )

    # TODO: Apply the headway on both outbound and inbound trips of a route.
    # In the final GTFS, make sure that every trip exported has a frequency entry

    valid_frequency_instances = (
        frequency_instances[
            ~frequency_instances.origin.str.lower().str.contains("test")
            & ~frequency_instances.destination.str.lower().str.contains("test")
        ]
        .query("status == 'finished'")[["trip_id", "interval", "avg_headway_sec"]]
        .groupby(["trip_id", "interval"])
        .agg({"avg_headway_sec": "mean"})
        .reset_index()
    )

    valid_frequency_instances.query("trip_id == 'Cy2rV_BRB3PfJEgiazMtI'")

    valid_frequency_instances.avg_headway_sec.isna().sum()

    intervals = sdi_intervals[['gid', 'observer_id']].rename({"gid":"interval_gid", "observer_id":"interval"}, axis=1)

    sdi_intervals

    trips = observer_trips.drop(['geom'], axis=1)

    trips = observer_trips.drop(['geom'], axis=1)
    trips["gid"] = range(1, len(trips) + 1)

    trips = trips.merge(intervals, how='cross')

    trips = trips.merge(valid_frequency_instances[["trip_id", "interval", "avg_headway_sec"]], left_on=['observer_id', 'interval'], right_on=['trip_id', 'interval'], how='left')

    trips['avg_headway_sec'] = trips.groupby(['observer_route_id', 'interval'], group_keys=True)['avg_headway_sec'].apply(
        lambda df: df.fillna(method='ffill').fillna(method='bfill')).reset_index(drop=True, level=[0,1])

    trips['avg_headway_agency_interval'] = trips.groupby(['agency_id', 'interval'])[
        'avg_headway_sec'].transform('mean')


    trips['ratio_headway'] = trips['avg_headway_sec'] / trips['avg_headway_agency_interval']

    len(trips[trips['trip_id'].isna()].observer_id.unique())

    trips['avg_ratio_headway'] = trips.groupby(['observer_route_id'])['ratio_headway'].transform('mean')

    trips["headway_estimation_method"] = np.where(trips['avg_headway_sec'].isna(), "from_similar_agency_interval", "from_own_freq_surveys")
    trips['final_headway'] = np.floor(trips['avg_headway_sec'].fillna(
        trips['avg_ratio_headway'] * trips['avg_headway_agency_interval']))

    (trips['avg_ratio_headway']).isna().sum()

    sdi_trips = observer_trips[['observer_id', 'geom']].merge(trips, on='observer_id').drop(['interval'],axis=1).drop_duplicates(subset=['observer_id'])

    # %%
    # Create trips table.
    # Add only the rows that have headway value

    run_ddl_sql("DROP table if exists transit.trips CASCADE;", city_db_con)

    run_ddl_sql("""
        create table if not exists transit.trips (
            gid serial primary key,
            o_id integer references transit.terminals(gid),
            d_id integer references transit.terminals(gid),
            agency_id integer,
            direction_id integer,
            observer_id text,
            observer_route_id text,
            fare real,
            agency_serial text,
            route_type integer DEFAULT 3,
            service_id text,
            geom geometry(LINESTRING, 4326)
        );
        """, city_db_con)


    # %%

    sdi_trips.sort_values("gid")[[
        "o_id",
        "d_id",
        "agency_id",
        "direction_id",
        "observer_id",
        "observer_route_id",
        "fare",
        "agency_serial",
        "service_id",
        'geom'
    ]].to_crs(4326).to_postgis(
        name="trips",
        schema="transit",
        con=city_db_con,
        if_exists="append",
        index=False
    )

    trips_intervals = trips[["gid", "final_headway", "interval_gid", "headway_estimation_method"]].rename(
        {"gid": "trip_id", "final_headway": "headway_secs", "interval_gid": "interval_id"}, axis=1
    )

    # %%

    run_ddl_sql("""
        drop table if exists transit.trips_intervals;
        create table if not exists transit.trips_intervals (
            gid SERIAL primary key,
            trip_id integer references transit.trips(gid),
            interval_id integer references transit.intervals(gid),
            headway_secs integer,
            headway_estimation_method text
        );
        """, city_db_con)

    trips_intervals.to_sql(
        name="trips_intervals",
        schema="transit",
        con=city_db_con,
        if_exists="append",
        index=False,
    )

    trips.sort_values(['observer_route_id', 'interval_gid'])[['observer_id', 'interval_gid', 'avg_headway_sec']]

    # %%
    ## Apply update trips view
    with open(sql_path("updated_trips_view.sql")) as f:
        trips_view_sql_query = f.read()
        run_ddl_sql(f"""
                DROP MATERIALIZED VIEW IF EXISTS transit.trips_view cascade;
                CREATE MATERIALIZED VIEW transit.trips_view AS
                {trips_view_sql_query};
                """, city_db_con)



    # %% [markdown]
    # ## Raw data
    '''LARGE SECTION BELOW WAS RELOCATED FROM THE BOTTOM TO UP HERE.
    THIS IS THE LAST SECTION IN script.ipynb'''

    # %%
    sdi_agencies = pd.read_sql(
        con=city_db_con, sql="""select * from transit.agencies"""
    )

    # %%
    onboard_instances = gpd.read_postgis(
        con=observer_db_con,
        sql=f"""select * from v_onboard_instances_ext where project_id='{OBSERVER_PROJECT_ID}' and deleted_at is null and geometry is not null and ST_IsValid(geometry)""",
        geom_col="geometry",
    ).assign(
    matched_geometry = lambda df: df.matched_geometry.apply(wkb.loads))

    onboard_instances = onboard_instances.drop(['agency_id'], axis=1).merge(
        sdi_agencies[["agency_id", "vehicle_id"]],
        left_on="agency",
        right_on="agency_id",
    ).drop(["agency_id"], axis=1)

    onboard_instances = onboard_instances.reset_index(names="gid")

    onboard_instances.to_postgis(
        name="onboard_instances",
        schema="raw",
        con=city_db_con,
    #     FIXME
        if_exists="replace",
        dtype={"gid": Integer, "geom": Geometry("GEOMETRY", srid=4326)},
    )

    run_ddl_sql("""
        SELECT addgeometrycolumn ('raw', 'onboard_instances','matched_geom',4326, 'LINESTRING', 2);
        UPDATE raw.onboard_instances SET matched_geom = st_geomfromtext(matched_geometry, 4326);
        ALTER TABLE raw.onboard_instances DROP COLUMN matched_geometry;
        ALTER TABLE raw.onboard_instances RENAME COLUMN matched_geom TO matched_geometry;
    """, city_db_con)

    # %% [markdown]
    # ### Onboard Stops

    # %%
    raw_stops = gpd.read_postgis(
        con=observer_db_con,
        sql=f"""
            WITH r as (select id as route_id from routes where setting_id='{observer_setting_id}'),
            t as (select id as trip_id from trips join r on trips.route_id = r.route_id),
            oi as (select id as oi_id, status, valid from onboard_instances join t on onboard_instances.trip_id = t.trip_id)
            select * from onboard_instance_stops tp join oi on tp.onboard_instance_id = oi.oi_id
            where deleted_at is null
            and geometry is not null
            and ST_IsValid(geometry)
            """,
        geom_col="geometry",
    )

    # %%
    raw_stops = raw_stops.rename(columns={'status':'parent_onboard_instance_status', 'valid':'parent_onboard_instance_valid'})

    # %%
    raw_stops = raw_stops.assign(
        board_male=lambda df: df.board_categorized.apply(lambda b: b['male'] if b else None),
        board_female=lambda df: df.board_categorized.apply(lambda b: b['female'] if b else None),
        alight_male=lambda df: df.alight_categorized.apply(lambda b: b['male'] if b else None),
        alight_female=lambda df: df.alight_categorized.apply(lambda b: b['female'] if b else None)
    )

    # raw_stops

    # %%
    # Columns: gid, name, geom, board, alight, observer_id, onboard_instance_observer_id, h3_index (calculated by Observer)

    run_ddl_sql("DROP table if exists raw.stops;", city_db_con)

    run_ddl_sql("""
        create table if not exists raw.stops (
            gid SERIAL primary key,
            name text,
            board integer,
            alight integer,
            board_male integer,
            alight_male integer,
            board_female integer,
            alight_female integer,
            observer_id text,
            onboard_instance_observer_id text,
            h3_index text,
            parent_onboard_instance_status text,
            parent_onboard_instance_valid boolean,
            geom GEOMETRY(GEOMETRY, 4326)
        );
        """, city_db_con) 

    raw_stops[
        ["name", "geometry", "onboard_instance_id", "id", "board", "alight",  "board_male", "alight_male",  "board_female", "alight_female", "h3_index", "parent_onboard_instance_status", "parent_onboard_instance_valid"]
    ].rename_geometry("geom").rename(
        columns={"onboard_instance_id": "onboard_instance_observer_id", "id": "observer_id"}
    ).to_postgis(
        name="stops",
        schema="raw",
        con=city_db_con,
        if_exists="append",
        dtype={"geom": Geometry("GEOMETRY", srid=4326)},
    )

    # %%
    # Columns: gid, name, geom, board, alight, observer_id, onboard_instance_observer_id, h3_index (calculated by Observer)

    run_ddl_sql("DROP table if exists raw.stops_created_at;", city_db_con) 
    

    run_ddl_sql( """
        create table if not exists raw.stops_created_at (
            observer_id text,
            created_at timestamptz
        );
        """, city_db_con) 
    
    raw_stops[
        ["id", "created_at"]
    ].rename(
        columns={ "id": "observer_id"}
    ).to_sql(
        name="stops_created_at",
        schema="raw",
        con=city_db_con,
        if_exists="replace",
        index=False
    )
    feedback.pushInfo("raw.stops is finished")

    # %% [markdown]
    # ### Onboard Trackpoints

    # %%
    trackpoints = gpd.read_postgis(
        con=observer_db_con,
        sql=f"""
            WITH r as (select id as route_id from routes where setting_id='{observer_setting_id}'),
            t as (select id as trip_id from trips join r on trips.route_id = r.route_id),
            oi as (select id as oi_id, status, valid from onboard_instances join t on onboard_instances.trip_id = t.trip_id)
            select * from onboard_instance_track_points tp join oi on tp.onboard_instance_id = oi.oi_id
            where deleted_at is null
            and geometry is not null
            and ST_IsValid(geometry)
            """,
        geom_col="geometry",
    )

    # %%
    trackpoints = trackpoints.rename(columns={'status':'onboard_instance_status', 'valid':'onboard_instance_valid'})

    # %%
    # Columns in city_cairo raw.trackpoints table: timestamp, geom, onboard_instance_id

    run_ddl_sql("DROP table if exists raw.trackpoints;", city_db_con) 

    run_ddl_sql("""
        create table if not exists raw.trackpoints (
            gid SERIAL primary key,
            timestamp timestamptz(0),
            geom GEOMETRY(GEOMETRY, 4326),
            onboard_instance_id text,
            onboard_instance_status text,
            onboard_instance_valid boolean
        );
        """, city_db_con) 


    trackpoints[["timestamp", "geometry", "onboard_instance_id", "onboard_instance_status", "onboard_instance_valid"]].rename_geometry(
        "geom"
    ).to_postgis(
        name="trackpoints",
        schema="raw",
        con=city_db_con,
        if_exists="append",
        dtype={"geom": Geometry("GEOMETRY", srid=4326)},
    )
    feedback.pushInfo("raw.trackpoints is finished")

    # %%
    # city_db_con.execute(text("DROP table if exists raw.indentification_instances"))
    run_ddl_sql("DROP table if exists raw.indentification_instances;", city_db_con)

    identification_instances = gpd.read_postgis(
        con=observer_db_con,
        sql=f"""select * from v_terminal_trips_ext where project_id='{OBSERVER_PROJECT_ID}' and deleted_at is null and geometry is not null and ST_IsValid(geometry)""",
        geom_col="geometry",
    )

    identification_instances = identification_instances.reset_index(names="gid")
    identification_instances.to_postgis(
        name="identification_instances",
        schema="raw",
        con=city_db_con,
        if_exists="append",
        dtype={"gid": Integer, "geom": Geometry("GEOMETRY", srid=4326)},
    )
    feedback.pushInfo("raw.identification is finished")

    # %%
    run_ddl_sql("DROP table if exists raw.frequency_instances;", city_db_con)

    frequency_instances = gpd.read_postgis(
        con=observer_db_con,
        sql=f"""select * from v_frequency_instances_ext where project_id='{OBSERVER_PROJECT_ID}' and deleted_at is null and geometry is not null and ST_IsValid(geometry)""",
        geom_col="geometry",
    )

    frequency_instances = frequency_instances.reset_index(names="gid")
    frequency_instances.to_postgis(
        name="frequency_instances",
        schema="raw",
        con=city_db_con,
        if_exists="append",
        dtype={"gid": Integer, "geom": Geometry("GEOMETRY", srid=4326)},
    )
    feedback.pushInfo("raw.frequency_instances is finished")

    # %%
    # TODO: LCTF continue from here. Waiting on stops digitization to be done.
    # UPDATE: Stops digitization automated by running this sql below
    with open(sql_path("create_processed_stops.sql")) as f:
        create_processed_stops_sql_query = f.read()
        run_ddl_sql(create_processed_stops_sql_query, city_db_con)


    with open(sql_path("trips_stops_sequence.sql")) as f:
            trips_stops_sequence_sql_query = f.read()
            run_ddl_sql(trips_stops_sequence_sql_query, city_db_con) # we used run_ddl_sql here bec. trips_stops_sequence.sql has CREATE in it
    with open(sql_path("od_stats.sql")) as f:
        od_stats_sql = f.read()
        # city_db_con.execute(text(od_stats_sql))
        run_ddl_sql(od_stats_sql, city_db_con)
    feedback.pushInfo("od_stats is finished")


    # %%
    # Reassign representative geometry to trips using updated algorithm (Frechet distance)
    city_db_con.execute(text( # I ADDED text() HERE
            f"""
            with recursive

    

            valid_oi as (
    select
        *
    from
        raw.onboard_instances oi
    where
        oi.status = 'finished'
        and oi.valid = 'True'
        and oi.matched_geometry is not null
            ),
            onboard_survey_sibling_pairs_initial as (
    select
        row_number() over () as idx,
        oi_1.trip_id,
        array[oi_1.id,
        oi_2.id] as id_pair
        --	st_frechetdistance(oi_1.matched_geometry, oi_2.matched_geometry)
    from
        valid_oi oi_1
    join valid_oi oi_2 on
        oi_1.trip_id = oi_2.trip_id
        and oi_1.id != oi_2.id
            ),
            onboard_survey_sibling_pairs_unnested as (
    select
        *,
        unnest(id_pair) as id_pair_unnested
    from
        onboard_survey_sibling_pairs_initial
    order by
        idx,
        id_pair_unnested
            ),
            onboard_survey_sibling_pairs_distinct as (
    select
        distinct trip_id,
        array_agg(id_pair_unnested) as id_pair
    from
        onboard_survey_sibling_pairs_unnested
    group by
        idx,
        trip_id
    order by
        trip_id
            ),
            frechet_dist as (
    select
        oi1.trip_id,
        p.id_pair,
        oi1.id as id1,
        oi2.id as id2,
        --	oi1.matched_geometry as geom1,
        --	oi2.matched_geometry as geom2,

    

        st_frechetdistance(st_transform(oi1.matched_geometry,
        3857),
        st_transform(oi2.matched_geometry,
        3857) ) as frechet_dist
    from
        onboard_survey_sibling_pairs_distinct p
    join valid_oi oi1 on
        p.id_pair[1] = oi1.id
    join valid_oi oi2 on
        p.id_pair[2] = oi2.id

    

            ),

    

            frechet_dist_avg as (
    select
        trip_id,
        unnest(id_pair) as oi_id,
        avg(frechet_dist) as frechet_dist_avg
    from
        frechet_dist
    group by
        oi_id,
        trip_id
    order by
        trip_id,
        frechet_dist_avg asc

    

            ),

    

            final_trip_geoms as (
    select
        distinct on
        (f.trip_id) f.trip_id,
        valid_oi.matched_geometry
    from
        frechet_dist_avg f
    join valid_oi on
        f.oi_id = valid_oi.id
    )update
        transit.trips t
    set
        geom = matched_geometry
    from
        final_trip_geoms g
    where
        t.observer_id = g.trip_id
        and g.matched_geometry is not null
                
                """)
        )

    try:
        ...
        feedback.pushInfo("Migration finished.")
    finally:
        try:
            observer_db_con.close()
        except Exception:
            pass
        try:
            city_db_con.close()
        except Exception:
            pass
        try:
            observer_db_engine.dispose()
        except Exception:
            pass
        try:
            sdi_db_engine.dispose()
        except Exception:
            pass

    '''LARGE SECTION ABOVE WAS RELOCATED FROM THE BOTTOM TO UP HERE.
    THIS IS THE LAST SECTION IN script.ipynb'''
