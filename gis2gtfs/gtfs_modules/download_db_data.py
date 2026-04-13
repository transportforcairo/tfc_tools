import os
import pandas as pd
import geopandas as gpd

# NOTE: QGIS loads plugins as packages under the plugin folder.
# Use relative imports so module resolution works regardless of sys.path.
from ...tfc_tools_common.sdi_io import SDISource, read_df, read_gdf

def download_db_data(sdi_mode, conn_name, gpkg_path, output_dir):
    """
    Export all relevant SDI tables to CSV and GeoJSON files for GTFS construction.
    Supports either:
      - Postgres (QGIS connection name), or
      - GeoPackage (tfc_rl_sdi_gpkg schema)

    Parameters:
        sdi_mode (str): 'postgres' or 'gpkg'
        conn_name (str|None): Name of the PostgreSQL connection in QGIS (postgres mode)
        gpkg_path (str|None): Path to GeoPackage (gpkg mode)
        output_dir (str): Directory to save exported files
    """

    source = SDISource(
        mode=sdi_mode,
        conn_name=conn_name,
        gpkg_path=gpkg_path,
    )
    source.validate()

    os.makedirs(output_dir, exist_ok=True)

    # -------------------- agency -------------------- #
    agency_df = read_df(source, "transit.agencies")
    vehicles_df = read_df(source, "transit.vehicles")

    # Robust join on vehicles.gid (not DataFrame index)
    if "gid" in vehicles_df.columns:
        agency_df = agency_df.merge(
            vehicles_df[[c for c in ["gid", "name", "passenger_capacity"] if c in vehicles_df.columns]],
            left_on="vehicle_id",
            right_on="gid",
            how="left",
            suffixes=("", "_veh"),
        )
        agency_df.drop(columns=["gid"], inplace=True, errors="ignore")
    else:
        # fallback (should not happen in TfC schema)
        agency_df = agency_df.merge(vehicles_df, left_on="vehicle_id", right_index=True, how="left")

    if "name" in agency_df.columns:
        agency_df.rename(columns={"name": "vehicle_name"}, inplace=True)
    elif "name_veh" in agency_df.columns:
        agency_df.rename(columns={"name_veh": "vehicle_name"}, inplace=True)

    agency_df.to_csv(os.path.join(output_dir, "agency.csv"), index=False)

    pickup_dropoff = agency_df[["agency_id", "has_serial"]].copy()
    pickup_dropoff["continuous_pickup"] = pickup_dropoff["has_serial"].apply(lambda x: 1 if x else 0)
    pickup_dropoff["continuous_drop_off"] = pickup_dropoff["has_serial"].apply(lambda x: 1 if x else 0)
    pickup_dropoff.drop(columns=["has_serial"], inplace=True)

    # -------------------- stops -------------------- #
    stops_gdf = read_gdf(source, "transit.stops")
    stops_gdf.to_file(os.path.join(output_dir, "stops.geojson"), driver="GeoJSON")

    # -------------------- trips -------------------- #
    trips_gdf = read_gdf(source, "transit.trips_view", where="geom IS NOT NULL" if source.mode == "postgres" else None)
    if source.mode == "gpkg":
        trips_gdf = trips_gdf[~trips_gdf.geometry.isna()].copy()
    trips_gdf = trips_gdf.merge(pickup_dropoff, on="agency_id", how="left")
    trips_gdf.to_file(os.path.join(output_dir, "trips.geojson"), driver="GeoJSON")

    # -------------------- terminals -------------------- #
    terminals_gdf = read_gdf(source, "transit.terminals")
    terminals_gdf.to_file(os.path.join(output_dir, "terminals.geojson"), driver="GeoJSON")

    # -------------------- intervals -------------------- #
    intervals_df = read_df(source, "transit.intervals")
    if "active" in intervals_df.columns:
        intervals_df = intervals_df[intervals_df["active"].isin([True, 1, "1", "t", "T", "true", "TRUE"])]
    intervals_df.to_csv(os.path.join(output_dir, "intervals.csv"), index=False)

    # -------------------- frequencies -------------------- #
    frequencies_df = read_df(source, "transit.trips_intervals")
    frequencies_df = frequencies_df[frequencies_df["interval_id"].isin(intervals_df["gid"])]
    frequencies_df.to_csv(os.path.join(output_dir, "frequencies.csv"), index=False)

    # -------------------- trip stop sequence -------------------- #
    tss_df = read_df(source, "transit.trip_stops_sequence")
    tss_df = tss_df.drop(columns=["gid"], errors="ignore")
    tss_df.to_csv(os.path.join(output_dir, "trip_stop_sequence.csv"), index=False)

    # -------------------- travel time OD stats -------------------- #
    travel_df = read_df(source, "transit.od_stats")
    travel_df = travel_df[["o_id", "d_id", "interval_id", "interval_start", "duration", "vehicle_name"]]
    travel_df.to_csv(os.path.join(output_dir, "travel_times_trackpoints.csv"), index=False)

    print("✅ All data successfully downloaded and saved.")
