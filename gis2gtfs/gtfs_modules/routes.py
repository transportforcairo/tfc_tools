'''
Main Steps:
1. Load trips data (trips.geojson):
Contains route-level attributes (e.g., route_id, agency_id, route_type, etc.)
2. Deduplicate:
Retain only one record per route_id (remove duplicates caused by origin/destination trip direction)
3. (Commented-out section): Intended logic to map numeric agency_id values to actual agency strings via agency.csv. Not used currently.
4. Add route_short_name:
Rename route_short to route_short_name
5. Include conditional logic:
If continuous_dropoff_pickup == TRUE, add continuous_pickup and continuous_drop_off columns.
Otherwise, exclude them.
6. Export routes.txt in GTFS format.
'''
import os
import pandas as pd
import geopandas as gpd


def generate(data_dir, data_raw_dir, continuous_dropoff_pickup=True):
    """
    Generate GTFS routes.txt file.

    Parameters:
        data_dir (str): Output GTFS folder
        data_raw_dir (str): Folder containing raw trips.geojson
        continuous_dropoff_pickup (bool): Whether to include continuous pickup/dropoff columns
    """
    print("📥 Loading trip data from trips.geojson...")

    trips_path = os.path.join(data_raw_dir, "trips.geojson")
    trips_gdf = gpd.read_file(trips_path)

    # Remove duplicate route_ids (keep first occurrence)
    routes_df = trips_gdf.drop_duplicates(subset=["route_id"]).copy()

    # Drop geometry
    routes_df = routes_df.drop(columns="geometry", errors="ignore")

    # Rename route_short to route_short_name
    routes_df.rename(columns={"route_short": "route_short_name"}, inplace=True)

    # Select columns
    if continuous_dropoff_pickup:
        required_columns = [
            "route_id", "agency_id", "route_long", "route_short_name",
            "route_type", "continuous_pickup", "continuous_drop_off"
        ]
    else:
        required_columns = [
            "route_id", "agency_id", "route_long", "route_short_name", "route_type"
        ]

    routes_df = routes_df[required_columns].rename(columns={"route_long": "route_long_name"})

    output_file = os.path.join(data_dir, "routes.txt")
    routes_df.to_csv(output_file, index=False)
    print("✅ routes.txt written.")
