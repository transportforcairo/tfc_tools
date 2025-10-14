'''
Main Steps
1. Load data:
trips.geojson from earlier (each row = a route geometry)
terminals.geojson (commented out)
intervals.csv to handle potential duplication (though unused here)
2. Filter duplicates:
Remove any duplicate trip_ids.
3. Assign trip headsign:
Use the destination column → trip_headsign.
4. Assign GTFS fields:
shape_id = trip_id + "_Shape"
Preserve: route_id, service_id, trip_headsign, direction_id, shape_id, trip_id
5. Drop geometry and export:
Save to trips.txt for GTFS
Save trips_with_intervals.csv for internal use
'''
import os
import pandas as pd
import geopandas as gpd


def generate(data_dir, data_raw_dir):
    """
    Generate GTFS trips.txt file from trips.geojson and prepare trips_with_intervals.csv.

    Parameters:
        data_dir (str): Output GTFS folder
        data_raw_dir (str): Input folder with geojson and intervals
    """
    print("🚌 Generating trips.txt...")

    trips_path = os.path.join(data_raw_dir, "trips.geojson")
    intervals_path = os.path.join(data_raw_dir, "intervals.csv")
    terminals_path = os.path.join(data_raw_dir, "terminals.geojson")

    # Load trips and remove duplicates
    trips_gdf = gpd.read_file(trips_path)
    trips_gdf["trip_id"] = trips_gdf["observer_id"].astype(str)  # Create trip_id from observer_id

    # Drop duplicates
    trips_gdf = trips_gdf.drop_duplicates(subset="trip_id").copy()

    # Assign trip_headsign from destination
    trips_gdf["trip_headsign"] = trips_gdf["destination"]

    # Assign shape_id
    trips_gdf["shape_id"] = trips_gdf["trip_id"].astype(str) + "_Shape"

    # Select required GTFS columns
    trips_df = trips_gdf[[
        "route_id", "service_id", "trip_headsign",
        "direction_id", "shape_id", "trip_id"
    ]].copy()

    # Drop geometry and raw ID column
    trips_df = trips_df.drop(columns=["geometry","observer_id"], errors="ignore")

    # Save to GTFS output
    trips_txt = os.path.join(data_dir, "trips.txt")
    trips_df.to_csv(trips_txt, index=False)

    # Save a CSV for internal use (optional)
    trips_csv = os.path.join(data_raw_dir, "trips_with_intervals.csv")
    trips_df.to_csv(trips_csv, index=False)

    print("✅ trips.txt and trips_with_intervals.csv written.")
