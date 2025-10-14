import os
import pandas as pd


def estimate_missing_values(travel_times: pd.DataFrame, vehicle: int = 0, od_pairs_buffer: pd.DataFrame = None) -> pd.DataFrame:
    """
    Fill missing travel time durations by computing OD-pair means.
    If vehicle == 1, disaggregate by vehicle_name.
    """
    travel_times["id"] = travel_times["o_id"].astype(str) + "-" + travel_times["d_id"].astype(str)

    if vehicle == 1:
        # Create all combinations of OD pair + interval + vehicle
        od_combinations = (
            travel_times[["id", "interval_id", "vehicle_name"]]
            .drop_duplicates()
            .merge(
                travel_times[["id", "o_id", "d_id"]].drop_duplicates(),
                on="id"
            )
            .drop(columns=["id"])
        )

        travel_times = (
            od_combinations
            .merge(travel_times, on=["o_id", "d_id", "interval_id", "vehicle_name"], how="left")
        )

        travel_times["duration"] = travel_times.groupby(["o_id", "d_id", "vehicle_name"])["duration"] \
            .transform(lambda x: x.fillna(x.mean()))

    else:
        # Create all combinations of OD pair + interval (ignore vehicle type)
        od_combinations = (
            travel_times[["id", "interval_id"]]
            .drop_duplicates()
            .merge(
                travel_times[["id", "o_id", "d_id"]].drop_duplicates(),
                on="id"
            )
            .drop(columns=["id"])
        )

        travel_times = (
            od_combinations
            .merge(travel_times, on=["o_id", "d_id", "interval_id"], how="left")
        )

        if od_pairs_buffer is not None:
            travel_times = travel_times.merge(
                od_pairs_buffer[["o_id", "d_id", "length_m"]],
                on=["o_id", "d_id"],
                how="left"
            )

        travel_times["duration"] = travel_times.groupby(["o_id", "d_id"])["duration"] \
            .transform(lambda x: x.fillna(x.mean()))

    return travel_times.drop(columns=["id"])


def run_estimate_missing_times(input_csv, output_csv, vehicle=0, od_pairs_buffer_csv=None):
    """
    Main wrapper function for running the missing time estimation
    """
    print("📥 Loading travel time data...")
    travel_times = pd.read_csv(input_csv)

    od_pairs_buffer = None
    if od_pairs_buffer_csv:
        od_pairs_buffer = pd.read_csv(od_pairs_buffer_csv)

    print("🧠 Estimating missing values...")
    travel_times_filled = estimate_missing_values(travel_times, vehicle=vehicle, od_pairs_buffer=od_pairs_buffer)

    print("💾 Saving filled travel times...")
    travel_times_filled.to_csv(output_csv, index=False)
    print("✅ Done: saved to", output_csv)
