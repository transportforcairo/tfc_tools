'''
Shared helpers used by trips.py, frequencies.py, and stop_times.py so they
all produce the SAME per-interval trip_id suffix tokens. Keeping this
logic in one place means trips.txt, frequencies.txt, and stop_times.txt
are guaranteed to agree -- they use the exact same suffix function.

Why a per-interval suffix at all?
---------------------------------
GTFS requires (trip_id, start_time) to be unique in frequencies.txt.
When the operator defines multiple intervals that share a start_time (or
overlap), each base trip produces colliding rows. The canonical GTFS fix
is to split the base trip into one GTFS trip per interval, appending a
suffix derived from the interval's name. This module implements the
suffix logic and the uniqueness guard.
'''
import os
import re
import pandas as pd


# Whitespace becomes '_'; runs of '_' collapse; leading/trailing '_' stripped.
# "WEEKDAY_ Morning peak" -> "WEEKDAY_Morning_peak".
_WS_RE = re.compile(r"\s+")
_MULTI_UNDERSCORE = re.compile(r"_+")


def _normalise_interval_name(name, gid):
    """Turn a free-form interval name into a trip_id-safe suffix token.

    Falls back to 'int<gid>' if the name is missing/empty so every interval
    still ends up with a deterministic, unique suffix.
    """
    if name is None or (isinstance(name, float) and pd.isna(name)):
        return f"int{gid}"
    token = _WS_RE.sub("_", str(name).strip())
    token = _MULTI_UNDERSCORE.sub("_", token).strip("_")
    return token or f"int{gid}"


def canonical_hms(value):
    """Normalise 'HH:MM', 'HH:MM:SS', or 'HH:MM:SS.ffffff' to 'HH:MM:SS'.

    Safe to call on NaN: returns NaN unchanged. Exists because Postgres
    `time` values round-trip through datetime with microseconds, and GTFS
    rejects anything but HH:MM:SS.
    """
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return value
    s = str(value).strip()
    if "." in s:
        s = s.split(".", 1)[0]
    parts = s.split(":")
    if len(parts) == 2:
        parts = parts + ["00"]
    h, m, sec = (int(x) for x in parts)
    return f"{h:02d}:{m:02d}:{sec:02d}"


def build_interval_suffix_map(intervals_df):
    """Return {interval_gid (str): suffix (str)} with unique suffixes.

    If two interval names normalise to the same token we disambiguate
    deterministically by appending '_<gid>' so every suffix is unique.
    Also asserts (start_time, end_time) is unique in the intervals table
    -- otherwise a frequencies row with a given (start, end) could map to
    multiple intervals, which downstream joins could not resolve.
    """
    if "gid" not in intervals_df.columns:
        raise ValueError("intervals is missing the 'gid' column")

    # First pass: candidate suffixes.
    candidates = {}
    for _, row in intervals_df.iterrows():
        gid = row["gid"]
        name = row.get("name") if "name" in intervals_df.columns else None
        candidates[str(gid)] = _normalise_interval_name(name, gid)

    # Detect collisions and disambiguate.
    seen = {}
    for gid, tok in candidates.items():
        seen.setdefault(tok, []).append(gid)
    resolved = {}
    for tok, gids in seen.items():
        if len(gids) == 1:
            resolved[gids[0]] = tok
        else:
            for gid in gids:
                resolved[gid] = f"{tok}_{gid}"
    return resolved


def load_intervals_and_suffixes(data_raw_dir):
    '''Read intervals.csv, filter to active rows, and return
    (intervals_df_with_suffix_column, suffix_by_gid_dict).

    If the day columns (mon..sun) are present, the returned dataframe
    also carries 'day_pattern' and 'service_id' columns derived from
    them. If they are missing (legacy intervals.csv), both columns are
    filled with None so callers can fall back to the legacy behaviour
    (hardcoded 'Ground_Daily').
    '''
    iv = pd.read_csv(os.path.join(data_raw_dir, "intervals.csv"), encoding="utf-8")
    if "active" in iv.columns:
        iv = iv[iv["active"].astype(str).str.lower().isin(
            ["true", "1", "t", "yes"]
        )].copy()
    suffix_by_gid = build_interval_suffix_map(iv)
    iv["suffix"] = iv["gid"].astype(str).map(suffix_by_gid)

    has_day_cols = all(c in iv.columns for c in DAY_COLS)
    if has_day_cols:
        iv = attach_day_patterns(iv)
    else:
        iv["day_pattern"] = None
        iv["service_id"] = None

    return iv, suffix_by_gid


def expanded_trip_id(base_trip_id, suffix):
    '''Compose the GTFS trip_id from the base trip id and an interval suffix.'''
    return f"{base_trip_id}_{suffix}"


# =========================================================================
# Day-pattern / service_id derivation
# =========================================================================
#
# Each row of transit_intervals carries seven boolean columns (mon..sun)
# saying which weekdays the interval runs on. That information must reach
# calendar.txt so consumers know when the service actually operates.
#
# A "day pattern" is just the 7-bit vector (mon, tue, wed, thu, fri, sat,
# sun) encoded as a 7-character string of '1'/'0'. Distinct patterns map
# to distinct service_ids. Intervals that share a pattern share a
# service_id.
#
# We emit nice names for the common patterns (daily, weekday, weekend,
# saturday, sunday, weekday+sat, weekday+sun) and fall back to
# 'svc_<pattern>' (e.g. 'svc_1010100') for irregular ones so the service
# name is still deterministic.

DAY_COLS = ("mon", "tue", "wed", "thu", "fri", "sat", "sun")

_PATTERN_NICE_NAMES = {
    "1111111": "svc_daily",
    "1111100": "svc_weekday",
    "0000011": "svc_weekend",
    "0000010": "svc_saturday",
    "0000001": "svc_sunday",
    "1111110": "svc_weekday_plus_sat",
    "1111101": "svc_weekday_plus_sun",
}


def _as_bool_int(value):
    '''Coerce a boolean-ish cell to 0/1.

    transit_intervals stores booleans as INTEGER 0/1 in SQLite, but
    intervals.csv may round-trip them as strings. Be lenient.
    '''
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return 0
    s = str(value).strip().lower()
    return 1 if s in ("1", "true", "t", "yes") else 0


def day_pattern_of(row):
    '''Return the 7-char '0'/'1' day-pattern string for an intervals row.

    row must be a pandas Series (or dict-like) exposing the columns
    mon..sun.
    '''
    missing = [c for c in DAY_COLS if c not in row]
    if missing:
        raise ValueError(
            f"intervals row is missing day column(s): {missing}. "
            "Expected columns: mon, tue, wed, thu, fri, sat, sun."
        )
    bits = [_as_bool_int(row[c]) for c in DAY_COLS]
    if sum(bits) == 0:
        raise ValueError(
            f"intervals row {row.get('gid', '?')} has no active days "
            "(mon..sun all zero). Refusing to emit a calendar row for a "
            "service that never runs."
        )
    return "".join(str(b) for b in bits)


def service_id_for_pattern(pattern):
    '''Map a 7-bit pattern to a deterministic service_id.

    Well-known patterns get readable names; everything else becomes
    'svc_<pattern>' so the service_id is still stable across runs.
    '''
    return _PATTERN_NICE_NAMES.get(pattern, f"svc_{pattern}")


def attach_day_patterns(intervals_df):
    '''Return a copy of intervals_df with two added columns:
    'day_pattern' (the 7-bit string) and 'service_id' (the canonical name).

    Raises ValueError if the day columns are missing -- the caller is
    expected to have confirmed they exist.
    '''
    out = intervals_df.copy()
    out["day_pattern"] = out.apply(day_pattern_of, axis=1)
    out["service_id"] = out["day_pattern"].map(service_id_for_pattern)
    return out


def distinct_service_patterns(intervals_df_with_patterns):
    '''Return a DataFrame with one row per distinct (service_id, day_pattern)
    pair found in the intervals. Used by calendar.py to emit one
    calendar.txt row per service actually referenced.
    '''
    return (
        intervals_df_with_patterns[["service_id", "day_pattern"]]
        .drop_duplicates()
        .sort_values("service_id")
        .reset_index(drop=True)
    )
