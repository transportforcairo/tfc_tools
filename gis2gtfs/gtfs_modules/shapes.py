'''
Main Steps:
1. Load geometry:
Reads trips.geojson, which contains LINESTRING geometries per trip
2. Create shape_id:
Combines trip_id with "_Shape" suffix
3. Explode geometries:
Uses the equivalent of R's st_cast(..., "POINT") to convert each linestring into a sequence of points
4. Add shape attributes:
shape_pt_sequence: point order
shape_pt_lat / shape_pt_lon: extracted from geometry
5. Drop geometry and export as shapes.txt
'''
import os
import pandas as pd
import geopandas as gpd
from shapely.geometry import LineString


def generate(data_dir, data_raw_dir):
    """
    Generate GTFS shapes.txt by extracting points from LINESTRING geometries
    in trips.geojson.

    Parameters:
        data_dir (str): Output directory for GTFS files
        data_raw_dir (str): Directory containing raw geojson files
    """
    print("🧭 Generating shapes.txt from trip geometries...")

    trips_path = os.path.join(data_raw_dir, "trips.geojson")
    output_file = os.path.join(data_dir, "shapes.txt")

    # Load trips.geojson
    trips_gdf = gpd.read_file(trips_path)

    # Ensure trip_id exists (from observer_id)
    trips_gdf["trip_id"] = trips_gdf["observer_id"].astype(str)

    # Create shape_id
    trips_gdf["shape_id"] = trips_gdf["trip_id"] + "_Shape"

    # --- Explode each linestring into vertices (manual st_cast equivalent) ---
    shape_records = []

    for _, row in trips_gdf.iterrows():
        shape_id = row["shape_id"]
        geom = row["geometry"]

        # Only handle LineString geometries
        if isinstance(geom, LineString):
            coords = list(geom.coords)
        else:
            continue  # skip non-LineString geometries

        for i, (lon, lat) in enumerate(coords, start=1): # this starts the counter at 1 not the default 0
            shape_records.append({
                "shape_id": shape_id,
                "shape_pt_lat": lat,
                "shape_pt_lon": lon,
                "shape_pt_sequence": i
            })

    # Create final DataFrame and save
    shapes_df = pd.DataFrame(shape_records)
    shapes_df.to_csv(output_file, index=False)

    print("✅ shapes.txt written.")
