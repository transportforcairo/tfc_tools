'''
Main Steps:
1. Load source files:
frequencies.csv: contains trip_id, interval_id, headway_secs
trips.geojson: contains geometries and gid (used to retrieve trip_id)
intervals.csv: contains interval start and end times keyed by gid
2. Join tables:
frequencies ← joined with trips on trip_id = gid to get trip_id
Then joined with intervals on interval_id = gid to get start_time, end_time
3. Select and rename:
Keep: trip_id, start_time, end_time, headway_secs
Ensure column names match GTFS
4. Output:
Save as frequencies.txt in GTFS format
'''
import os
import pandas as pd
import geopandas as gpd


def generate(data_dir, data_raw_dir):
    """
    Generate GTFS frequencies.txt by joining trips and intervals to frequency definitions.

    Parameters:
        data_dir (str): Output folder for GTFS files
        data_raw_dir (str): Input folder containing raw CSV and GeoJSON data
    """
    print("📥 Loading raw frequency data...")

    frequencies_path = os.path.join(data_raw_dir, "frequencies.csv")
    trips_path = os.path.join(data_raw_dir, "trips.geojson")
    intervals_path = os.path.join(data_raw_dir, "intervals.csv")
    output_file = os.path.join(data_dir, "frequencies.txt")

    # Load inputs
    freq_df = pd.read_csv(frequencies_path)
    freq_df.rename(columns={"trip_id": "trip_id_freq"}, inplace=True) # rename field to avoid conflict with trips.geojson field of trip_id
    trips_gdf = gpd.read_file(trips_path)[["gid","observer_id"]].copy() # copy the field "gid" and "observer_id" then rename it as "trip_id"
    trips_gdf["trip_id"] = trips_gdf["observer_id"].astype(str)  # same as R: stringr::str_c(gid)
    intervals_df = pd.read_csv(intervals_path)

    # Drop geometries for joining
    trips_df = trips_gdf.drop(columns=["geometry","observer_id"], errors="ignore")

    # Join trips using: frequencies.trip_id_freq = trips.gid
    freq_df = freq_df.merge(trips_df, left_on="trip_id_freq", right_on="gid", how="inner")
    
    # Join intervals: frequencies.interval_id = intervals.gid
    freq_df = freq_df.merge(intervals_df, left_on="interval_id", right_on="gid", how="left")

    # Select
    freq_df = freq_df[["trip_id", "start_time", "end_time", "headway_secs"]]

    # # 🚨 TEMP PATCH: Drop rows with missing headway_secs temporarily (avoid error during conversion)
    # if freq_df["headway_secs"].isnull().any():
    #     print("⚠️ Warning: Dropping rows with missing headway_secs (temporary fix)")
    #     freq_df = freq_df.dropna(subset=["headway_secs"])
    # # 🚨 END OF TEMP PATCH   

    # # 🚨 FUTURE PATCH: Raise error if any missing headway_secs values are found
    # if freq_df["headway_secs"].isnull().any():
    #     raise ValueError("❌ Error: transit.frequencies table has missing values in 'headway_secs'. Please fix the SDI data.")
    # # 🚨 END OF FUTURE PATCH

    freq_df["headway_secs"] = freq_df["headway_secs"].astype(int) # Write headway_secs as integer (no decimals)
    
    # Save GTFS output
    freq_df.to_csv(output_file, index=False)
    print("✅ frequencies.txt written.")
