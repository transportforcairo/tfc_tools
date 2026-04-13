"""TfC Tools - Processing widget wrapper helpers

This module contains a custom Processing widget wrapper used to *visually*
toggle (enable/disable) parameter widgets based on an enum selection.

In QGIS Processing algorithm dialogs, it is standard UX practice that
mutually-exclusive inputs are visually disabled when not applicable.

We implement this using the Processing 'widget_wrapper' metadata mechanism,
which allows attaching a custom wrapper to a parameter without building a
full custom dialog.
"""

from __future__ import annotations

from typing import Dict, List, Optional

from qgis.PyQt.QtWidgets import QComboBox

from qgis.core import QgsProcessingParameterEnum

# processing.gui.wrappers.WidgetWrapper is the supported base
# for metadata-defined wrappers in QGIS Processing dialogs.
from processing.gui.wrappers import WidgetWrapper


class SDISourceToggleWidgetWrapper(WidgetWrapper):
    """Enum wrapper which toggles (enables/disables) other parameter widgets.

    The wrapper expects an ``enable_map`` dict (passed via parameter metadata)
    with the following structure:

        {
          "0": ["PARAM_A", "PARAM_B"],
          "1": ["PARAM_C"],
        }

    Where keys are enum indices as strings, and values are lists of parameter
    names which should be enabled when that enum option is selected.

    All parameters mentioned across all lists are treated as the target set.
    Any target parameter not in the enabled list for the current index will be
    disabled.
    """

    def __init__(
        self,
        param,
        dialog,
        row: int = 0,
        col: int = 0,
        enable_map: Optional[Dict[str, List[str]]] = None,
        **kwargs,
    ):
        self._enable_map: Dict[str, List[str]] = enable_map or {}
        self._wrappers_by_name = {}
        super().__init__(param, dialog, row=row, col=col, **kwargs)

    # ---- WidgetWrapper API -------------------------------------------------

    def createWidget(self, **kwargs):
        combo = QComboBox()

        # Populate from enum options (fallback to empty list)
        opts: List[str] = []
        pdef = self.parameterDefinition()
        if isinstance(pdef, QgsProcessingParameterEnum):
            opts = [str(o) for o in (pdef.options() or [])]

        for o in opts:
            combo.addItem(o)

        combo.currentIndexChanged.connect(lambda *_: self._on_value_changed())
        return combo

    def setValue(self, value):
        # Enum parameters are represented by their index
        if value is None:
            return
        try:
            self.widget.setCurrentIndex(int(value))
        except Exception:
            pass

    def value(self):
        return int(self.widget.currentIndex())

    # ---- Cross-parameter wiring -------------------------------------------

    def postInitialize(self, wrappers):
        # Build lookup for related wrappers by parameter name
        self._wrappers_by_name = {
            w.parameterDefinition().name(): w
            for w in wrappers
            if w is not None and w.parameterDefinition() is not None
        }
        self._apply_toggle()

    # ---- Internals ---------------------------------------------------------

    def _on_value_changed(self):
        self._apply_toggle()
        # Notify Processing dialog that this parameter changed
        self.widgetValueHasChanged.emit(self)

    def _apply_toggle(self):
        # Determine which targets exist
        targets = set()
        for enabled in self._enable_map.values():
            targets.update(enabled)

        if not targets:
            return

        idx = str(self.value())
        enabled_now = set(self._enable_map.get(idx, []))

        for name in targets:
            wrapper = self._wrappers_by_name.get(name)
            if not wrapper:
                continue

            is_enabled = name in enabled_now

            # Disable/enable both label and widget for a clear visual cue
            try:
                w = wrapper.wrappedWidget()
                if w is not None:
                    w.setEnabled(is_enabled)
            except Exception:
                # Fallback for wrappers which may not expose wrappedWidget
                try:
                    if hasattr(wrapper, "widget") and wrapper.widget is not None:
                        wrapper.widget.setEnabled(is_enabled)
                except Exception:
                    pass

            try:
                lbl = wrapper.wrappedLabel()
                if lbl is not None:
                    lbl.setEnabled(is_enabled)
            except Exception:
                # Some wrappers embed labels inside widgets
                pass
