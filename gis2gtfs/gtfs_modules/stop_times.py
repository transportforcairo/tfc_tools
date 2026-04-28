'''
Build stop_times.txt with per-interval travel times.

Change vs the previous version
------------------------------
Previously this module pinned every stop-row to:
    interval_id       = intervals.loc[0, 'gid']             # line 66
    interval_start    = intervals.loc[0, 'start_time']      # line 57
...which meant the merge against travel_times only ever matched morning-peak
rows, every trip's arrival_time was anchored to 06:30:00, and every other
interval in travel_times_trackpoints_filled_na.csv was ignored. On top of
that, trip_id in the output was the BASE trip_id, so any per-interval
variation could not have coexisted anyway.

This version:
1. Cross-joins trip_stop_sequence x active intervals (one stop-row per
   interval) so every row carries its OWN interval_id + interval_start.
2. Composes the EXPANDED trip_id (<base>_<suffix>) using the same helper
   as trips.py and frequencies.py, so all three stay consistent.
3. Merges per-interval OD durations -- interval 1 rows match interval 1
   OD stats, interval 3 rows match interval 3, etc.
4. Uses each interval's interval_start as the time base for that trip's
   stop_times (not a global anchor).
'''
import os
import pandas as pd
import numpy as np
from datetime import timedelta

from ._interval_expansion import (
    load_intervals_and_suffixes,
    canonical_hms,
    expanded_trip_id,
)


def generate(data_dir, data_raw_dir, dwell_time_sec=15):
    '''Generate GTFS stop_times.txt using per-interval OD travel times.'''
    print("Generating stop_times.txt (per-interval)...")

    intervals_df, _suffix_by_gid = load_intervals_and_suffixes(data_raw_dir)
    travel_times = pd.read_csv(
        os.path.join(data_raw_dir, "travel_times_trackpoints_filled_na.csv"),
        encoding="utf-8",
    )
    stop_seq = pd.read_csv(
        os.path.join(data_raw_dir, "trip_stop_sequence.csv"), encoding="utf-8"
    )
    output_path = os.path.join(data_dir, "stop_times.txt")

    # Schema-mismatch fix: some exports use 'observer_trip_id' instead of
    # 'trip_id'. Rename so the rest of this function treats the column as
    # the BASE trip id (pre-suffix).
    if "trip_id" not in stop_seq.columns and "observer_trip_id" in stop_seq.columns:
        stop_seq = stop_seq.rename(columns={"observer_trip_id": "trip_id"})
    stop_seq = stop_seq.rename(columns={"trip_id": "base_trip_id"})
    stop_seq["base_trip_id"] = stop_seq["base_trip_id"].astype(str)

    # Cross-join stops x intervals. Every (trip, stop) becomes N rows, one
    # per active interval, each carrying its own interval_id / interval_start
    # / suffix. This is the fix: the old code collapsed this to N=1.
    iv_keep = intervals_df[["gid", "start_time", "suffix"]].copy()
    iv_keep = iv_keep.rename(columns={"gid": "interval_id", "start_time": "interval_start"})
    iv_keep["interval_id"] = iv_keep["interval_id"].astype(str)
    iv_keep["interval_start"] = iv_keep["interval_start"].apply(canonical_hms)

    expanded = stop_seq.merge(iv_keep, how="cross")

    # Compose the expanded trip_id used across trips.txt / frequencies.txt.
    expanded["trip_id"] = [
        expanded_trip_id(bt, suf)
        for bt, suf in zip(expanded["base_trip_id"], expanded["suffix"])
    ]

    # from_id + segment length computed PER expanded trip. The group key has
    # to be the expanded trip_id, not the base, so the shift/diff don't leak
    # across intervals.
    expanded = expanded.sort_values(by=["trip_id", "stop_sequence"]).reset_index(drop=True)
    # IMPORTANT: shift() on an int column promotes to float (because NaN needs
    # it), and later astype(str) would then emit "136.0" instead of "136",
    # which would never match travel_times.o_id = "136". Cast stop_id to str
    # BEFORE shifting so from_id is already a proper text column.
    expanded["stop_id"] = expanded["stop_id"].astype(str)
    expanded["from_id"] = expanded.groupby("trip_id")["stop_id"].shift(1)
    expanded["length"] = expanded.groupby("trip_id")["distance"].diff().fillna(0)

    # Merge keys all coerced to str on both sides. from_id already is str
    # (it inherited from stop_id); interval_id and vehicle_name are set
    # explicitly here.
    expanded = expanded.astype({
        "interval_id": str, "vehicle_name": str,
    })
    tt = travel_times.copy()
    tt["interval_id"] = tt["interval_id"].astype(str)
    tt["o_id"] = tt["o_id"].astype(str)
    tt["d_id"] = tt["d_id"].astype(str)
    tt["vehicle_name"] = tt["vehicle_name"].astype(str)
    # The OD table can carry multiple rows per (o_id, d_id, interval_id,
    # vehicle_name) when the same stop pair participates in multiple trips
    # with slightly different geometries. Duration is effectively identical
    # across those duplicates. Deduplicate so each stop_times row gets
    # exactly one merge match -- otherwise a single segment fans out to N
    # rows in the output.
    tt = tt.drop_duplicates(subset=["o_id", "d_id", "interval_id", "vehicle_name"])

    merged = pd.merge(
        expanded,
        tt[["o_id", "d_id", "interval_id", "vehicle_name", "duration"]],
        how="left",
        left_on=["from_id", "stop_id", "interval_id", "vehicle_name"],
        right_on=["o_id", "d_id", "interval_id", "vehicle_name"],
    )

    # Fallback: static 10 m/s where the OD table has no matching row.
    # stop_sequence == 1 is legitimately null (it's the origin).
    missing_mask = merged["duration"].isna() & (merged["stop_sequence"] != 1)
    fallback_count = int(missing_mask.sum())
    merged.loc[missing_mask, "duration"] = np.ceil(merged.loc[missing_mask, "length"] / 10)
    merged["duration"] = merged["duration"].fillna(0)
    if fallback_count:
        print(f"  [info] {fallback_count} segment(s) used static 10 m/s fallback.")

    # Arrival / departure built per expanded trip using ITS own interval_start.
    merged["cumulative_duration"] = merged.groupby("trip_id")["duration"].cumsum()
    base_time = pd.to_timedelta(merged["interval_start"])
    merged["arrival_time"] = (
        base_time
        + pd.to_timedelta(merged["cumulative_duration"], unit="s")
        + pd.to_timedelta((merged["stop_sequence"] - 1) * dwell_time_sec, unit="s")
    )
    merged["departure_time"] = merged["arrival_time"] + pd.to_timedelta(dwell_time_sec, unit="s")

    def _fmt(td_series):
        # GTFS stop_times allow hours >= 24 to represent service that
        # continues past midnight on the same service day (e.g. a trip
        # that departs at 23:55 and ends at 24:10 must NOT wrap to 00:10,
        # otherwise the validator flags `start_and_end_range_out_of_order`
        # and `stop_time_with_arrival_before_previous_departure_time`).
        #
        # Timedelta.components.hours wraps at 24 (the rest goes into the
        # `days` field), so building HH:MM:SS from components silently
        # drops the day count. Use total_seconds() instead so hours can
        # legitimately exceed 23.
        total_secs = td_series.dt.total_seconds().astype("int64")
        hours, rem = divmod(total_secs, 3600)
        minutes, seconds = divmod(rem, 60)
        return pd.Series(
            [f"{h:02d}:{m:02d}:{s:02d}" for h, m, s in zip(hours, minutes, seconds)],
            index=td_series.index,
        )
    merged["arrival_time"] = _fmt(merged["arrival_time"])
    merged["departure_time"] = _fmt(merged["departure_time"])

    stop_times = merged[
        ["trip_id", "stop_id", "stop_sequence", "arrival_time", "departure_time"]
    ].sort_values(by=["trip_id", "stop_sequence"])

    stop_times.to_csv(output_path, index=False, encoding="utf-8")
    print(f"stop_times.txt written ({len(stop_times)} rows).")
