'''
Generate calendar.txt dynamically from the day-of-week flags carried on
each interval in transit_intervals (mon, tue, wed, thu, fri, sat, sun).

How it works
------------
- Read intervals.csv and filter to active rows (same rule as everywhere
  else in the pipeline).
- For each interval, compute its 7-bit day pattern and a canonical
  service_id (e.g. 'svc_weekday', 'svc_weekend', 'svc_daily', or the
  fallback 'svc_<pattern>' for irregular patterns).
- Emit one calendar.txt row per distinct service_id actually referenced
  by the intervals. No hardcoded 'Ground_Daily' anymore.

Behaviour on older exports
--------------------------
If intervals.csv lacks the mon..sun columns (legacy exports produced
before the refresh plugin carried those fields), this module falls back
to the previous behaviour: emit one row for `fallback_service_id` that
runs every day. That keeps the module backwards-compatible while new
exports automatically benefit from the accurate calendar.
'''
import os
import pandas as pd

from ._interval_expansion import (
    DAY_COLS,
    attach_day_patterns,
    distinct_service_patterns,
)


def generate(data_dir, data_raw_dir, start_date, end_date,
             fallback_service_id="Ground_Daily"):
    """Generate GTFS calendar.txt.

    Parameters
    ----------
    data_dir : str
        Output folder where calendar.txt is written.
    data_raw_dir : str
        Folder containing intervals.csv (the raw export produced earlier
        in the pipeline by download_db_data).
    start_date, end_date : int or str
        YYYYMMDD values for calendar.start_date / end_date, applied to
        every service row.
    fallback_service_id : str
        Service id used when intervals.csv has no day columns (legacy
        behaviour). In that case the row runs every day of the week.
    """
    output_file = os.path.join(data_dir, "calendar.txt")
    intervals_path = os.path.join(data_raw_dir, "intervals.csv")
    iv = pd.read_csv(intervals_path, encoding="utf-8")

    # Honour the active flag like every other module does.
    if "active" in iv.columns:
        iv = iv[iv["active"].astype(str).str.lower().isin(
            ["true", "1", "t", "yes"]
        )].copy()

    # Detect whether the export includes day flags.
    has_day_cols = all(c in iv.columns for c in DAY_COLS)

    if not has_day_cols:
        # Legacy fallback: one row that runs every day, keyed on the
        # caller-provided fallback_service_id. This keeps existing
        # downstream tooling working on older exports.
        calendar_df = pd.DataFrame([{
            "service_id": fallback_service_id,
            "monday": 1, "tuesday": 1, "wednesday": 1, "thursday": 1,
            "friday": 1, "saturday": 1, "sunday": 1,
            "start_date": start_date,
            "end_date": end_date,
        }])
        calendar_df.to_csv(output_file, index=False, encoding="utf-8")
        print(
            f"calendar.txt written (legacy single-service mode: "
            f"'{fallback_service_id}', runs every day)."
        )
        return

    # New behaviour: derive one service per distinct day pattern.
    iv_p = attach_day_patterns(iv)
    services = distinct_service_patterns(iv_p)

    rows = []
    for _, row in services.iterrows():
        pattern = row["day_pattern"]
        rows.append({
            "service_id": row["service_id"],
            "monday":    int(pattern[0]),
            "tuesday":   int(pattern[1]),
            "wednesday": int(pattern[2]),
            "thursday":  int(pattern[3]),
            "friday":    int(pattern[4]),
            "saturday":  int(pattern[5]),
            "sunday":    int(pattern[6]),
            "start_date": start_date,
            "end_date": end_date,
        })

    calendar_df = pd.DataFrame(rows)
    calendar_df.to_csv(output_file, index=False, encoding="utf-8")
    print(
        f"calendar.txt written ({len(calendar_df)} service(s) derived "
        f"from intervals.mon-sun)."
    )
