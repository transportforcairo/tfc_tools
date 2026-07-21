#!/usr/bin/env python3
"""
gtfs_attributes.py — pre-flight validator for GTFS-bound source attributes.

Part of the "edit GIS layers → re-export GTFS" UX (see
``GTFS_Edit_Model_Specification_v1.docx``, section 5). The Refresh SDI derived
layers tool calls this BEFORE it rebuilds the derived layers, so that a user who
edited the source layers in QGIS is told about GTFS-breaking values (bad
route_type, invalid timezone, duplicate keys, orphaned foreign keys, ambiguous
intervals) up front instead of discovering them as a cryptic GTFS validator
error after export.

Scope
-----
This checks the *source* (authored) layers only — the columns a user edits. It
does NOT re-implement the geometric derivation (that is the Refresh tool's job)
and it does NOT touch the derived layers.

Design
------
Depends on pandas only (already a plugin dependency). No QGIS, shapely or fiona
imports, so it can be unit-tested in a bare Python environment. Geometry is
irrelevant to these checks and is ignored if present.

Severity
--------
* ``ERROR`` — will produce an invalid or misleading GTFS feed (duplicate keys,
  invalid enum, orphaned foreign key, ambiguous intervals).
* ``WARN``  — questionable but exportable (blank names, odd fare, missing URL).

The caller decides what to do with them. The Refresh tool reports every issue as
processing feedback and, only in *strict* mode, aborts when any ERROR is present.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Optional

import pandas as pd

ERROR = "ERROR"
WARN = "WARN"

# GTFS route_type: standard set + the extended-route-types range (100–1899).
_STD_ROUTE_TYPES = {0, 1, 2, 3, 4, 5, 6, 7, 11, 12}

_URL_RE = re.compile(r"^https?://.+", re.IGNORECASE)
_TIME_RE = re.compile(r"^\d{1,3}:[0-5]\d(:[0-5]\d)?$")
# Suffix normalisation, identical to gis2gtfs/_interval_expansion.py.
_WS_RE = re.compile(r"\s+")
_MULTI_US = re.compile(r"_+")

_TRUE = {"1", "true", "t", "yes", "y"}
_FALSE = {"0", "false", "f", "no", "n", "", "nan", "none"}


@dataclass
class Issue:
    severity: str
    layer: str
    field: str
    message: str
    count: int = 0
    sample: Optional[list] = None

    def format(self) -> str:
        loc = f"{self.layer}.{self.field}" if self.field else self.layer
        head = f"[{self.severity}] {loc}: {self.message}"
        if self.count:
            head += f" ({self.count} row(s))"
        if self.sample:
            head += f" e.g. {self.sample[:5]}"
        return head


@dataclass
class Result:
    issues: list = field(default_factory=list)
    skipped: list = field(default_factory=list)  # human notes about missing tables

    def add(self, severity, layer, fieldname, message, count=0, sample=None):
        self.issues.append(Issue(severity, layer, fieldname, message, count, sample))

    @property
    def errors(self):
        return [i for i in self.issues if i.severity == ERROR]

    @property
    def warnings(self):
        return [i for i in self.issues if i.severity == WARN]

    @property
    def has_errors(self):
        return bool(self.errors)

    def to_lines(self):
        return [i.format() for i in self.issues]


# --------------------------------------------------------------------------- #
# Coercion helpers
# --------------------------------------------------------------------------- #
def _col(df, *names):
    """Return the first present column name (case-insensitive) or None."""
    if df is None:
        return None
    lower = {c.lower(): c for c in df.columns}
    for n in names:
        if n.lower() in lower:
            return lower[n.lower()]
    return None


def _s(df, colname):
    """String series, stripped, with NaN preserved as ''."""
    return df[colname].map(lambda v: "" if v is None or (isinstance(v, float) and pd.isna(v)) else str(v).strip())


def _is_blank(series):
    return series.map(lambda v: v == "")


def _boolish_bad(series):
    """Rows whose value is not recognisably boolean."""
    def bad(v):
        if v is None or (isinstance(v, float) and pd.isna(v)):
            return False
        return str(v).strip().lower() not in (_TRUE | _FALSE)
    return series.map(bad)


def _active_mask(df):
    ac = _col(df, "active")
    if ac is None:
        return pd.Series([True] * len(df), index=df.index)
    return df[ac].map(lambda v: str(v).strip().lower() in _TRUE)


def _norm_suffix(name):
    tok = _WS_RE.sub("_", str(name).strip())
    return _MULTI_US.sub("_", tok).strip("_")


# --------------------------------------------------------------------------- #
# Per-layer checks
# --------------------------------------------------------------------------- #
def _check_stops(r, stops):
    if stops is None or stops.empty:
        return
    sid = _col(stops, "stop_id")
    if sid:
        s = _s(stops, sid)
        blank = _is_blank(s)
        if blank.any():
            r.add(ERROR, "transit_stops", "stop_id", "null/blank stop_id", int(blank.sum()))
        dup = s[~blank].duplicated(keep=False)
        if dup.any():
            r.add(ERROR, "transit_stops", "stop_id", "duplicate stop_id",
                  int(dup.sum()), sorted(set(s[~blank][dup]))[:5])
    name = _col(stops, "stop_name", "name")
    if name:
        blank = _is_blank(_s(stops, name))
        if blank.any():
            r.add(WARN, "transit_stops", "stop_name", "blank stop_name", int(blank.sum()))


def _check_terminals(r, terminals):
    if terminals is None or terminals.empty:
        return
    name = _col(terminals, "name")
    if name:
        blank = _is_blank(_s(terminals, name))
        if blank.any():
            r.add(WARN, "transit_terminals", "name",
                  "blank terminal name (feeds route_long_name / headsign)", int(blank.sum()))


def _check_agencies(r, agencies):
    if agencies is None or agencies.empty:
        return
    aid = _col(agencies, "agency_id")
    if aid:
        s = _s(agencies, aid)
        blank = _is_blank(s)
        if blank.any():
            r.add(ERROR, "transit_agencies", "agency_id", "null/blank agency_id", int(blank.sum()))
        dup = s[~blank].duplicated(keep=False)
        if dup.any():
            r.add(ERROR, "transit_agencies", "agency_id", "duplicate agency_id",
                  int(dup.sum()), sorted(set(s[~blank][dup]))[:5])
    aname = _col(agencies, "agency_name", "common_name", "name")
    if aname:
        blank = _is_blank(_s(agencies, aname))
        if blank.any():
            r.add(WARN, "transit_agencies", "agency_name", "blank agency_name", int(blank.sum()))
    url = _col(agencies, "agency_url")
    if url:
        s = _s(agencies, url)
        bad = s.map(lambda v: v != "" and not _URL_RE.match(v))
        if bad.any():
            r.add(WARN, "transit_agencies", "agency_url", "malformed URL (expected http/https)",
                  int(bad.sum()), list(s[bad])[:5])
    tz = _col(agencies, "agency_timezone")
    if tz:
        s = _s(agencies, tz)
        valid = _tz_set()
        if valid is None:
            r.add(WARN, "transit_agencies", "agency_timezone",
                  "timezone validation skipped (zoneinfo unavailable)")
        else:
            bad = s.map(lambda v: v != "" and v not in valid)
            if bad.any():
                r.add(ERROR, "transit_agencies", "agency_timezone", "invalid IANA timezone",
                      int(bad.sum()), sorted(set(s[bad]))[:5])
    hs = _col(agencies, "has_serial")
    if hs and _boolish_bad(agencies[hs]).any():
        r.add(WARN, "transit_agencies", "has_serial", "non-boolean value",
              int(_boolish_bad(agencies[hs]).sum()))


def _check_vehicles(r, vehicles):
    if vehicles is None or vehicles.empty:
        return
    name = _col(vehicles, "name")
    if name:
        blank = _is_blank(_s(vehicles, name))
        if blank.any():
            r.add(WARN, "transit_vehicles", "name", "blank vehicle name", int(blank.sum()))
    cap = _col(vehicles, "passenger_capacity")
    if cap:
        num = pd.to_numeric(vehicles[cap], errors="coerce")
        bad = (num.notna() & (num < 0)) | ((vehicles[cap].notna()) & num.isna())
        if bad.any():
            r.add(WARN, "transit_vehicles", "passenger_capacity",
                  "non-integer or negative capacity", int(bad.sum()))


def _check_trips_view(r, tv):
    if tv is None or tv.empty:
        return
    rid = _col(tv, "route_id", "observer_route_id")
    if rid:
        blank = _is_blank(_s(tv, rid))
        if blank.any():
            r.add(ERROR, "transit_trips_view", "route_id", "null/blank route_id", int(blank.sum()))
    rt = _col(tv, "route_type")
    if rt:
        num = pd.to_numeric(tv[rt], errors="coerce")
        def _bad_rt(v):
            if pd.isna(v):
                return False
            iv = int(v)
            return not (iv in _STD_ROUTE_TYPES or (100 <= iv <= 1899))
        bad = num.map(_bad_rt)
        if bad.any():
            r.add(ERROR, "transit_trips_view", "route_type", "invalid GTFS route_type",
                  int(bad.sum()), sorted(set(tv[rt][bad].tolist()))[:5])
    did = _col(tv, "direction_id")
    if did:
        num = pd.to_numeric(tv[did], errors="coerce")
        bad = num.notna() & ~num.isin([0, 1])
        if bad.any():
            r.add(ERROR, "transit_trips_view", "direction_id", "direction_id must be 0 or 1",
                  int(bad.sum()), sorted(set(tv[did][bad].tolist()))[:5])
    fare = _col(tv, "fare")
    if fare:
        num = pd.to_numeric(tv[fare], errors="coerce")
        bad = (num.notna() & (num < 0)) | (tv[fare].notna() & num.isna())
        if bad.any():
            r.add(WARN, "transit_trips_view", "fare", "non-numeric or negative fare", int(bad.sum()))


def _check_intervals(r, intervals):
    if intervals is None or intervals.empty:
        r.add(ERROR, "transit_intervals", "", "no intervals present; feed cannot be built")
        return
    active = _active_mask(intervals)
    if not active.any():
        r.add(ERROR, "transit_intervals", "active", "no active intervals; feed will be empty")

    for tcol in ("start_time", "end_time"):
        c = _col(intervals, tcol)
        if c:
            s = _s(intervals, c)
            bad = s.map(lambda v: v != "" and not _TIME_RE.match(v.split(".")[0]))
            if bad.any():
                r.add(ERROR, "transit_intervals", tcol, "unparseable time (expected HH:MM[:SS])",
                      int(bad.sum()), list(s[bad])[:5])

    for d in ("mon", "tue", "wed", "thu", "fri", "sat", "sun"):
        c = _col(intervals, d)
        if c and _boolish_bad(intervals[c]).any():
            r.add(WARN, "transit_intervals", d, "non-boolean day flag",
                  int(_boolish_bad(intervals[c]).sum()))

    name = _col(intervals, "name", "interval_name")
    if name:
        blank = _is_blank(_s(intervals, name))
        if blank.any():
            r.add(WARN, "transit_intervals", "name", "blank interval name", int(blank.sum()))

    # §10 — (start_time, end_time) must be unique among ACTIVE intervals.
    sc, ec = _col(intervals, "start_time"), _col(intervals, "end_time")
    if sc and ec:
        act = intervals[active]
        pairs = list(zip(_s(act, sc), _s(act, ec)))
        seen, dups = set(), set()
        for p in pairs:
            (dups if p in seen else seen).add(p)
        if dups:
            r.add(ERROR, "transit_intervals", "start_time,end_time",
                  "active intervals share the same (start,end); frequencies joins are ambiguous",
                  len(dups), [f"{a}-{b}" for a, b in list(dups)[:5]])

    # §10 — advisory: distinct names that normalise to the same trip_id suffix.
    if name:
        act = intervals[active]
        tok_map = {}
        for nm in _s(act, name):
            if nm:
                tok_map.setdefault(_norm_suffix(nm), set()).add(nm)
        collided = {k: v for k, v in tok_map.items() if len(v) > 1}
        if collided:
            sample = [sorted(v) for v in list(collided.values())[:3]]
            r.add(WARN, "transit_intervals", "name",
                  "distinct interval names normalise to the same trip_id suffix "
                  "(auto-disambiguated by gid, but rename for cleaner ids)",
                  len(collided), sample)


def _check_trips_intervals(r, ti, intervals, trips_view):
    if ti is None or ti.empty:
        return
    hw = _col(ti, "headway_secs")
    if hw:
        num = pd.to_numeric(ti[hw], errors="coerce")
        bad = (num.notna() & (num <= 0)) | (ti[hw].notna() & num.isna())
        if bad.any():
            r.add(ERROR, "transit_trips_intervals", "headway_secs",
                  "headway_secs must be a positive integer", int(bad.sum()))
        if num.isna().any():
            r.add(WARN, "transit_trips_intervals", "headway_secs",
                  "null headway_secs (row will be dropped or need a fallback)", int(num.isna().sum()))
    tcol, icol = _col(ti, "trip_id"), _col(ti, "interval_id")
    if tcol and icol:
        pair_dup = ti.duplicated(subset=[tcol, icol], keep=False)
        if pair_dup.any():
            r.add(ERROR, "transit_trips_intervals", "trip_id,interval_id",
                  "duplicate (trip_id, interval_id) frequency rows", int(pair_dup.sum()))


# --------------------------------------------------------------------------- #
# Cross-layer referential integrity (§9)
# --------------------------------------------------------------------------- #
def _orphans(child, child_col, parent, *parent_cols):
    """Return child values whose key is absent from all parent_cols."""
    if child is None or parent is None:
        return None
    cc = _col(child, child_col)
    if cc is None:
        return None
    valid = set()
    for pc in parent_cols:
        col = _col(parent, pc)
        if col is not None:
            valid |= set(_s(parent, col))
    cs = _s(child, cc)
    miss = cs.map(lambda v: v != "" and v not in valid)
    return cs[miss]


def _check_referential(r, T):
    ag, veh, iv, ti, tv, term = (T.get(k) for k in
                                 ("agencies", "vehicles", "intervals", "trips_intervals",
                                  "trips_view", "terminals"))

    o = _orphans(ag, "vehicle_id", veh, "gid")
    if o is not None and len(o):
        r.add(ERROR, "transit_agencies", "vehicle_id",
              "vehicle_id not found in transit_vehicles.gid", len(o), sorted(set(o))[:5])

    o = _orphans(tv, "agency_id", ag, "agency_id", "gid")
    if o is not None and len(o):
        r.add(ERROR, "transit_trips_view", "agency_id",
              "agency_id not found in transit_agencies", len(o), sorted(set(o))[:5])

    for k in ("o_id", "d_id"):
        o = _orphans(tv, k, term, "gid")
        if o is not None and len(o):
            r.add(WARN, "transit_trips_view", k,
                  "terminal reference not found in transit_terminals.gid", len(o), sorted(set(o))[:5])

    o = _orphans(ti, "interval_id", iv, "gid")
    if o is not None and len(o):
        r.add(ERROR, "transit_trips_intervals", "interval_id",
              "interval_id not found in transit_intervals.gid", len(o), sorted(set(o))[:5])

    o = _orphans(ti, "trip_id", tv, "gid")
    if o is not None and len(o):
        r.add(ERROR, "transit_trips_intervals", "trip_id",
              "trip_id not found in transit_trips_view.gid", len(o), sorted(set(o))[:5])


# --------------------------------------------------------------------------- #
# Optional-dep helpers
# --------------------------------------------------------------------------- #
_TZ_CACHE = False


def _tz_set():
    """Return the set of valid IANA timezones, or None if unavailable."""
    global _TZ_CACHE
    if _TZ_CACHE is False:
        try:
            from zoneinfo import available_timezones
            _TZ_CACHE = set(available_timezones())
        except Exception:
            _TZ_CACHE = None
    return _TZ_CACHE


# --------------------------------------------------------------------------- #
# Public entry point
# --------------------------------------------------------------------------- #
def validate_gtfs_attributes(tables: dict, *, strict: bool = False) -> Result:
    """Validate GTFS-bound source attributes.

    Parameters
    ----------
    tables : dict[str, pandas.DataFrame | None]
        Logical name → DataFrame. Recognised keys:
        ``agencies, vehicles, intervals, trips_intervals, stops, terminals,
        trips_view``. Missing / None entries are skipped (noted in
        ``Result.skipped``), never raised — older exports lack some tables.
    strict : bool
        Only affects nothing here; it is echoed back for the caller's
        convenience. Severity classification is independent of strictness.

    Returns
    -------
    Result
    """
    r = Result()
    for key in ("agencies", "vehicles", "intervals", "trips_intervals",
                "stops", "terminals", "trips_view"):
        if tables.get(key) is None:
            r.skipped.append(key)

    _check_stops(r, tables.get("stops"))
    _check_terminals(r, tables.get("terminals"))
    _check_agencies(r, tables.get("agencies"))
    _check_vehicles(r, tables.get("vehicles"))
    _check_trips_view(r, tables.get("trips_view"))
    _check_intervals(r, tables.get("intervals"))
    _check_trips_intervals(r, tables.get("trips_intervals"),
                           tables.get("intervals"), tables.get("trips_view"))
    _check_referential(r, tables)
    return r
