"""tfc_tools_common.layer_styles

Embed default QGIS attribute-form definitions INSIDE a GeoPackage so the
"edit GIS layers → re-export GTFS" UX travels with the data, independently of
whatever .qgs/.qgz project the user builds.

Mechanism
---------
QGIS stores per-layer styles (symbology AND the attribute-form config: per-field
read-only flags, constraints, widgets) in a standard ``layer_styles`` table
inside the GeoPackage. When a layer that has a default style row there
(``useAsDefault=1``) is added to ANY project, QGIS auto-applies it. So instead of
shipping loose .qml files the user must manage, our export/refresh tools write
the default form directly into the gpkg.

What we embed
-------------
* ``transit_trips_view`` — the derived columns are set read-only so users can't
  waste effort hand-editing values the Refresh tool authoritatively recomputes
  (see GTFS_Edit_Model_Specification, section 6). Users change the true inputs
  (terminals, agency common_name/serial) instead.
* Source key fields — not-null + unique constraints on ``transit_stops.stop_id``
  and ``transit_agencies.agency_id`` so bad keys can't be committed via the form.

This is a UI guardrail, not enforcement: a determined user can toggle a field
editable or edit via SQL. That's acceptable — the real guarantee is the Refresh
recompute (overwrites bypassed edits) plus the validation pre-flight. This just
keeps honest users from editing fields that won't stick.

The function is best-effort: any failure is reported and swallowed so it can
never break an export or refresh. It must run inside QGIS (needs QgsVectorLayer).
"""

from __future__ import annotations

import os
import contextlib

# Derived columns on transit_trips_view that the Refresh tool recomputes.
# Locking them prevents silent overwrite-on-refresh of hand edits.
LOCKED_TRIPS_VIEW_COLS = [
    "route_short", "route_long", "trip_short", "len_km",
    "origin", "destination", "vehicle_name", "passenger_capacity",
]

# Single-field key constraints: layer -> {field: (not_null, unique)}.
KEY_CONSTRAINTS = {
    "transit_stops": {"stop_id": (True, True)},
    "transit_agencies": {"agency_id": (True, True)},
}

_STYLE_NAME = "tfc_edit_model"
_STYLE_DESC = ("TfC edit-model default form: read-only derived trips_view "
               "columns; not-null/unique key constraints.")


def embed_default_forms(gpkg_path: str, feedback=None) -> None:
    """Write default attribute-form styles into the gpkg's layer_styles table.

    Safe to call repeatedly; existing TfC style rows are replaced. Never raises.
    """
    def _log(msg):
        if feedback is not None:
            with contextlib.suppress(Exception):
                feedback.pushInfo(msg)

    def _warn(msg):
        if feedback is not None:
            with contextlib.suppress(Exception):
                feedback.reportError(msg)

    if not gpkg_path or not os.path.exists(gpkg_path):
        return

    try:
        from qgis.core import QgsVectorLayer, QgsFieldConstraints
    except Exception as e:  # pragma: no cover - only importable inside QGIS
        _warn(f"Could not embed edit-form styles (QGIS API unavailable): {e}")
        return

    # Only layers we actually want to style.
    targets = ["transit_trips_view"] + list(KEY_CONSTRAINTS.keys())
    styled = []

    for layer_name in targets:
        uri = f"{gpkg_path}|layername={layer_name}"
        try:
            lyr = QgsVectorLayer(uri, layer_name, "ogr")
            if not lyr.isValid():
                continue  # layer absent from this export — fine

            applied = False

            # 1) read-only derived columns (trips_view only)
            if layer_name == "transit_trips_view":
                cfg = lyr.editFormConfig()
                for fname in LOCKED_TRIPS_VIEW_COLS:
                    idx = lyr.fields().indexOf(fname)
                    if idx >= 0:
                        cfg.setReadOnly(idx, True)
                        applied = True
                lyr.setEditFormConfig(cfg)

            # 2) key constraints (hard: block commit if violated)
            for fname, (not_null, unique) in KEY_CONSTRAINTS.get(layer_name, {}).items():
                idx = lyr.fields().indexOf(fname)
                if idx < 0:
                    continue
                if not_null:
                    lyr.setFieldConstraint(
                        idx, QgsFieldConstraints.Constraint.ConstraintNotNull,
                        QgsFieldConstraints.ConstraintStrength.ConstraintStrengthHard)
                    applied = True
                if unique:
                    lyr.setFieldConstraint(
                        idx, QgsFieldConstraints.Constraint.ConstraintUnique,
                        QgsFieldConstraints.ConstraintStrength.ConstraintStrengthHard)
                    applied = True

            if not applied:
                continue

            _replace_tfc_style(lyr)
            lyr.saveStyleToDatabase(_STYLE_NAME, _STYLE_DESC, True, "")
            styled.append(layer_name)
        except Exception as e:
            _warn(f"Could not embed edit-form style for {layer_name}: {e}")

    if styled:
        _log(f"Embedded default edit-form styles for: {', '.join(styled)}")


def _replace_tfc_style(lyr) -> None:
    """Delete any prior TfC style rows for this layer so re-runs don't pile up."""
    try:
        res = lyr.listStylesInDatabase()
        # SIP binding returns (relatedCount, ids, names, descriptions, msgError)
        ids, names = None, None
        if isinstance(res, (tuple, list)):
            for part in res:
                if isinstance(part, (list, tuple)) and part and isinstance(part[0], str):
                    if ids is None:
                        ids = list(part)
                    elif names is None:
                        names = list(part)
                        break
        if ids and names:
            for sid, nm in zip(ids, names):
                if nm == _STYLE_NAME:
                    with contextlib.suppress(Exception):
                        lyr.deleteStyleFromDatabase(sid)
    except Exception:
        # Style listing/deletion is best-effort cleanup; a failure here just
        # means a duplicate style row, which is harmless.
        return
