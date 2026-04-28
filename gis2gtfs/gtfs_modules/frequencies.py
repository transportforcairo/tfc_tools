'''
Build frequencies.txt with one row per (expanded_trip_id, interval). The
expanded trip_id carries the per-interval suffix, so (trip_id, start_time)
is unique by construction and the GTFS duplicate_key error goes away.

Inputs
- frequencies.csv: raw transit.trips_intervals table, has (trip_id, interval_id, headway_secs)
- trips.geojson: used to map the transit.trips.gid -> observer_id (= base trip_id)
- intervals.csv: used for (interval start/end) and for the suffix name

The merge chain here is identical to the legacy file; the only change is
that the final trip_id is '<base>_<suffix>' instead of just '<base>'.
'''
import os
import pandas as pd
import geopandas as gpd

from ._interval_expansion import load_intervals_and_suffixes, canonical_hms, expanded_trip_id


def generate(data_dir, data_raw_dir):
    '''Generate frequencies.txt with expanded per-interval trip_ids.'''
    print("Loading raw frequency data...")

    freq_path = os.path.join(data_raw_dir, "frequencies.csv")
    trips_path = os.path.join(data_raw_dir, "trips.geojson")
    output_file = os.path.join(data_dir, "frequencies.txt")

    freq_df = pd.read_csv(freq_path, encoding="utf-8")
    # Rename the incoming FK (which refers to transit.trips.gid) so it
    # doesn't shadow the real GTFS trip_id column we're going to build.
    freq_df = freq_df.rename(columns={"trip_id": "trip_gid"})

    trips_gdf = gpd.read_file(trips_path)[["gid", "observer_id"]].copy()
    trips_gdf["base_trip_id"] = trips_gdf["observer_id"].astype(str)

    # Active intervals + suffix map.
    iv, _suffix_by_gid = load_intervals_and_suffixes(data_raw_dir)
    active_ids = set(iv["gid"].astype(str))
    freq_df["interval_id"] = freq_df["interval_id"].astype(str)
    freq_df = freq_df[freq_df["interval_id"].isin(active_ids)].copy()

    # Attach base trip_id (transit.trips.gid -> observer_id).
    freq_df = freq_df.merge(
        trips_gdf[["gid", "base_trip_id"]],
        left_on="trip_gid", right_on="gid", how="inner",
    )

    # Attach interval timing + suffix.
    iv_keep = iv[["gid", "start_time", "end_time", "suffix"]].copy()
    iv_keep["gid"] = iv_keep["gid"].astype(str)
    freq_df = freq_df.merge(
        iv_keep, left_on="interval_id", right_on="gid",
        how="left", suffixes=("", "_iv"),
    )

    # Build the expanded GTFS trip_id.
    freq_df["trip_id"] = [
        expanded_trip_id(bt, suf)
        for bt, suf in zip(freq_df["base_trip_id"], freq_df["suffix"])
    ]

    # Canonical HH:MM:SS for GTFS.
    freq_df["start_time"] = freq_df["start_time"].apply(canonical_hms)
    freq_df["end_time"] = freq_df["end_time"].apply(canonical_hms)

    freq_df = freq_df[["trip_id", "start_time", "end_time", "headway_secs"]].copy()

    # Guardrail: (trip_id, start_time) must be unique in frequencies.txt.
    dup = freq_df.duplicated(subset=["trip_id", "start_time"], keep=False)
    if dup.any():
        sample = freq_df.loc[dup].head().to_dict("records")
        raise ValueError(
            "frequencies.txt has duplicate (trip_id, start_time) rows even "
            f"after interval expansion. Sample: {sample}"
        )

    freq_df["headway_secs"] = freq_df["headway_secs"].astype(int)
    freq_df.to_csv(output_file, index=False, encoding="utf-8")
    print(f"frequencies.txt written ({len(freq_df)} rows).")
