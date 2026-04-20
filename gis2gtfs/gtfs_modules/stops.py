'''
Main Steps
1. Load stops.geojson — created earlier from the DB.
2. Filter out empty geometries.
3. Derive stop_lat and stop_lon FROM THE POINT GEOMETRY ONLY.
   The point geometry is the single source of truth. Any stop_lat/stop_lon
   columns that may happen to exist on the input are ignored, because they
   can silently desync from the geometry when stops are moved in QGIS.
4. Emit GTFS-compliant stops.txt with stop_id, stop_name, stop_lat, stop_lon.
'''
import os
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

    stops_gdf = gpd.read_file(input_file)

    # Remove empty/null geometries — they cannot produce valid GTFS coordinates.
    before = len(stops_gdf)
    stops_gdf = stops_gdf[stops_gdf.geometry.notna() & ~stops_gdf.geometry.is_empty].copy()
    dropped = before - len(stops_gdf)
    if dropped:
        print(f"⚠️  Dropped {dropped} stop(s) with null/empty geometry.")

    # Ensure CRS is WGS84 (EPSG:4326) for GTFS lat/lon coordinates.
    # GeoJSON per RFC 7946 is always 4326; this is defensive in case the loader
    # returns a different CRS or None.
    if stops_gdf.crs is None:
        stops_gdf = stops_gdf.set_crs("EPSG:4326")
    elif str(stops_gdf.crs).upper() not in ("EPSG:4326", "OGC:CRS84"):
        stops_gdf = stops_gdf.to_crs("EPSG:4326")

    # Derive stop_lat/stop_lon STRICTLY from the point geometry.
    # Never trust any pre-existing stop_lat/stop_lon columns — they are known
    # to desync from the geometry when stops are edited visually in QGIS.
    stops_gdf = stops_gdf.drop(columns=["stop_lat", "stop_lon"], errors="ignore")
    stops_gdf["stop_lat"] = stops_gdf.geometry.y
    stops_gdf["stop_lon"] = stops_gdf.geometry.x

    # Select the required GTFS columns and write.
    stops_gdf = stops_gdf[["stop_id", "stop_name", "stop_lat", "stop_lon"]]
    stops_gdf.to_csv(output_file, index=False, encoding='utf-8')
    print("✅ stops.txt written.")
