import os
import pandas as pd
import geopandas as gpd
from sqlalchemy import create_engine
import urllib.parse
from qgis.core import QgsProviderRegistry
from qgis.core import QgsProviderRegistry, QgsDataSourceUri

# NOTE: If you got an error on QGIS saying: ModuleNotFoundError: No module named 'sqlalchemy'
# the solution is to do the following: 
# Open OSGeo4W Shell 
# Run the following command inside the OSGeo4W shell: python3 -m pip install sqlalchemy
# and You can confirm it worked by running: python3 -c "import sqlalchemy; print(sqlalchemy.__version__)"


def download_db_data(conn_name, output_dir):
    """
    Connect to the SDI PostgreSQL database and export all relevant transit tables
    to CSV and GeoJSON files for GTFS construction.

    Parameters:
        conn_name (str): Name of the PostgreSQL connection in QGIS
        output_dir (str): Directory to save exported files
    """

    # 🔍 Look up the PostgreSQL connection from QGIS
    provider = QgsProviderRegistry.instance().providerMetadata("postgres")
    connection = provider.findConnection(conn_name)

    if connection is None:
        raise ValueError(f"Could not find PostgreSQL connection named '{conn_name}' in QGIS.")

    # uri = connection.uri()
    uri_str = connection.uri()
    uri = QgsDataSourceUri(uri_str)  # ✅ Parse the URI string properly
    db_name = uri.database()
    db_user = uri.username()
    db_password = uri.password()
    db_host = uri.host()
    db_port = uri.port()

    # Create output folder
    os.makedirs(output_dir, exist_ok=True)

    # 👇 encode the password to be URL-safe 
    # (This ensures special characters like @, &, %, : are treated correctly by the URL parser.)
    safe_password = urllib.parse.quote_plus(db_password)

    # Connect to the database using SQLAlchemy (safe connection string)
    conn_str = f"postgresql://{db_user}:{safe_password}@{db_host}:{db_port}/{db_name}"
    engine = create_engine(conn_str)
    conn = engine.connect()
    print("✅ Connected to the database.")

    # -------------------- agency -------------------- #
    agency_df = pd.read_sql_table("agencies", con=engine, schema="transit")
    vehicles_df = pd.read_sql_table("vehicles", con=engine, schema="transit")
    vehicles_df = vehicles_df.drop(columns="gid", errors="ignore") # this line is to avoid having "gid" in vehicles and agency, and the errors="ignore" was during plugin modifications
    agency_df = agency_df.merge(vehicles_df, left_on="vehicle_id", right_index=True, how="left")
    agency_df.rename(columns={"name": "vehicle_name"}, inplace=True)
    agency_df.to_csv(os.path.join(output_dir, "agency.csv"), index=False)

    pickup_dropoff = agency_df[["agency_id", "has_serial"]].copy()
    pickup_dropoff["continuous_pickup"] = pickup_dropoff["has_serial"].apply(lambda x: 1 if x else 0)
    pickup_dropoff["continuous_drop_off"] = pickup_dropoff["has_serial"].apply(lambda x: 1 if x else 0)
    pickup_dropoff.drop(columns=["has_serial"], inplace=True)

    # -------------------- stops -------------------- #
    stops_query = "SELECT * FROM transit.stops"
    stops_gdf = gpd.read_postgis(stops_query, con=conn)
    stops_gdf.to_file(os.path.join(output_dir, "stops.geojson"), driver="GeoJSON")

    # -------------------- trips -------------------- #
    trips_query = "SELECT * FROM transit.trips_view WHERE geom IS NOT NULL"
    trips_gdf = gpd.read_postgis(trips_query, con=conn)
    trips_gdf = trips_gdf.merge(pickup_dropoff, on="agency_id", how="left")
    trips_gdf.to_file(os.path.join(output_dir, "trips.geojson"), driver="GeoJSON")

    # -------------------- terminals -------------------- #
    terminals_query = "SELECT * FROM transit.terminals"
    terminals_gdf = gpd.read_postgis(terminals_query, con=conn)
    terminals_gdf.to_file(os.path.join(output_dir, "terminals.geojson"), driver="GeoJSON")

    # -------------------- intervals -------------------- #
    intervals_query = "SELECT * FROM transit.intervals WHERE active = TRUE"
    intervals_df = pd.read_sql_query(intervals_query, con=engine)
    intervals_df.to_csv(os.path.join(output_dir, "intervals.csv"), index=False)

    # -------------------- frequencies -------------------- #
    frequencies_df = pd.read_sql_table("trips_intervals", con=engine, schema="transit")
    frequencies_df = frequencies_df[frequencies_df["interval_id"].isin(intervals_df["gid"])]
    frequencies_df.to_csv(os.path.join(output_dir, "frequencies.csv"), index=False)

    # -------------------- trip stop sequence -------------------- #
    tss_df = pd.read_sql_table("trip_stops_sequence", con=engine, schema="transit")
    tss_df = tss_df.drop(columns=["gid"], errors="ignore")
    tss_df.to_csv(os.path.join(output_dir, "trip_stop_sequence.csv"), index=False)

    # -------------------- travel time OD stats -------------------- #
    travel_df = pd.read_sql_table("od_stats", con=engine, schema="transit")
    travel_df = travel_df[["o_id", "d_id", "interval_id", "interval_start", "duration", "vehicle_name"]]
    travel_df.to_csv(os.path.join(output_dir, "travel_times_trackpoints.csv"), index=False)

    conn.close()
    print("✅ All data successfully downloaded and saved.")
