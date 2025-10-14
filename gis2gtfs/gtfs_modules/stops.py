'''
Main Steps
1. Load stops.geojson — created earlier from the DB.
2. Filter:
Remove any features with empty geometries using st_is_empty().
# 3. Select and rename:
# Keep only gid and stop_name, and rename gid → stop_id.
# 4. Extract coordinates:
# 5. stop_lat and stop_lon extracted from point geometry.
'''
import os
import pandas as pd
import geopandas as gpd


def generate(data_dir, data_raw_dir):
    """
    Generate GTFS stops.txt file from stops.geojson.

    Parameters:
        data_dir (str): Output folder for GTFS
        data_raw_dir (str): Folder containing raw stops.geojson
    """
    print("📍 Generating stops.txt...")

    input_file = os.path.join(data_raw_dir, "stops.geojson")
    output_file = os.path.join(data_dir, "stops.txt")

    # Read as GeoDataFrame
    stops_gdf = gpd.read_file(input_file)

    # Remove empty geometries
    stops_gdf = stops_gdf[~stops_gdf.geometry.is_empty].copy()

    # Extract coordinates BEFORE dropping geometry
    # stops_gdf["stop_lat"] = stops_gdf.geometry.y
    # stops_gdf["stop_lon"] = stops_gdf.geometry.x

    # Select relevant columns and rename
    # stops_gdf = stops_gdf[["gid", "stop_name", "stop_lat", "stop_lon"]]
    stops_gdf = stops_gdf[["stop_id", "stop_name", "stop_lat", "stop_lon"]]
    # stops_gdf = stops_gdf.rename(columns={"gid": "stop_id"})

    # Save
    stops_gdf.to_csv(output_file, index=False)
    print("✅ stops.txt written.")
