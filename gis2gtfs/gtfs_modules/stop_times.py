'''
This script builds the stop_times.txt GTFS file, which specifies:
- Sequence of stops for each trip
- Arrival/departure times at each stop

Main Steps
🔹 1. Load Required Data
intervals.csv: list of time intervals (start time, ID)
travel_times_trackpoints_filled_na.csv: OD-pair travel times
trip_stop_sequence.csv: sequences of stops for each trip
agency.csv: provides vehicle_name, passenger_capacity (used in joining)
🔹 2. Duplicate trips for all time intervals
For each stop in a trip, replicate it across intervals
Assign one interval ID and start time to each trip copy (currently: the first interval only)
🔹 3. Add travel time between stops
Add a from_id column (previous stop)
Join OD travel time by matching: from_id, stop_id, interval_id, and vehicle_name
🔹 4. Fix or remove missing travel times
Fill missing durations using static speed (10 m/s) × distance
Discard trips with multiple missing times (disabled in current version)
🔹 5. Calculate arrival and departure times
Compute cumulative travel time
Add dwell time per stop
Convert everything to HH:MM:SS
Round to the nearest second
🔹 6. Export stop_times.txt
'''
import os
import pandas as pd
import numpy as np
from datetime import timedelta


def generate(data_dir, data_raw_dir, dwell_time_sec=15):
    """
    Generate GTFS stop_times.txt using trip_stop_sequence, travel_times, and intervals.

    Parameters:
        data_dir (str): Output folder for GTFS files
        data_raw_dir (str): Folder containing raw inputs
        dwell_time_sec (int): Default dwell time in seconds
    """
    print("🕓 Generating stop_times.txt...")

    # --- Load data
    intervals = pd.read_csv(os.path.join(data_raw_dir, "intervals.csv"))
    travel_times = pd.read_csv(os.path.join(data_raw_dir, "travel_times_trackpoints_filled_na.csv"))
    stop_seq = pd.read_csv(os.path.join(data_raw_dir, "trip_stop_sequence.csv"))
    agency = pd.read_csv(os.path.join(data_raw_dir, "agency.csv"))[["gid", "agency_id", "vehicle_name"]]
    output_path = os.path.join(data_dir, "stop_times.txt")

    # ✅ Fix for schema mismatch: observer_trip_id → trip_id
    if "trip_id" not in stop_seq.columns and "observer_trip_id" in stop_seq.columns:
        stop_seq = stop_seq.rename(columns={"observer_trip_id": "trip_id"})

    # --- TEMP: Attach static interval_start from first row of intervals.csv (start_time field)
    # stop_seq["interval_start"] = str(intervals.loc[0, "start_time"])
    static_start_time = str(intervals.loc[0, "start_time"])
    # Normalize to HH:MM:SS
    if len(static_start_time.split(":")) == 2:
        static_start_time += ":00"
    stop_seq["interval_start_custom"] = static_start_time

    stop_seq["interval_id"] = str(intervals.loc[0, "gid"])

    # ensure trip_id is str
    stop_seq["trip_id"] = stop_seq["trip_id"].astype(str)

    # --- Add from_id and segment length
    stop_seq = stop_seq.sort_values(by=["trip_id", "stop_sequence"])
    stop_seq["from_id"] = stop_seq.groupby("trip_id")["stop_id"].shift(1)
    stop_seq["length"] = stop_seq.groupby("trip_id")["distance"].diff().fillna(0)

    # --- Join travel times
    travel_times["interval_start"] = travel_times["interval_start"].astype(str)
    travel_times["o_id"] = travel_times["o_id"].astype(str)
    travel_times["d_id"] = travel_times["d_id"].astype(str)

    stop_seq = stop_seq.astype({
        "stop_id": str, "from_id": str, "interval_id": float, "vehicle_name": str
    })

    merged = pd.merge(
        stop_seq,
        travel_times,
        how="left",
        left_on=["from_id", "stop_id", "interval_id", "vehicle_name"],
        right_on=["o_id", "d_id", "interval_id", "vehicle_name"]
    )

    # --- Fill missing durations with static speed (10 m/s)
    merged["duration"] = merged.apply(
        lambda row: np.ceil(row["length"] / 10) if pd.isna(row["duration"]) and row["stop_sequence"] != 1 else row["duration"],
        axis=1
    ).fillna(0)

    # --- Calculate cumulative time and arrival/departure
    def to_hms(seconds):
        return (timedelta(seconds=int(seconds)))

    merged["cumulative_duration"] = merged.groupby("trip_id")["duration"].cumsum()

    # Parse interval_start as timedelta
    base_time = pd.to_timedelta(merged["interval_start_custom"]) #This ensures that the version of interval_start used for calculating arrival_time isn't overwritten during the merge.

    merged["arrival_time"] = base_time + pd.to_timedelta(merged["cumulative_duration"], unit='s') + \
                             pd.to_timedelta((merged["stop_sequence"] - 1) * dwell_time_sec, unit='s')
    merged["departure_time"] = merged["arrival_time"] + pd.to_timedelta(dwell_time_sec, unit='s')

    # Format times as HH:MM:SS
    merged["arrival_time"] = merged["arrival_time"].dt.components.apply(
        lambda row: f"{int(row.hours):02}:{int(row.minutes):02}:{int(row.seconds):02}", axis=1
    )
    merged["departure_time"] = merged["departure_time"].dt.components.apply(
        lambda row: f"{int(row.hours):02}:{int(row.minutes):02}:{int(row.seconds):02}", axis=1
    )

    # --- Final output
    stop_times = merged[["trip_id", "stop_id", "stop_sequence", "arrival_time", "departure_time"]]
    stop_times = stop_times.sort_values(by=["trip_id", "stop_sequence"])

    stop_times.to_csv(output_path, index=False)

    print("✅ stop_times.txt written.")
