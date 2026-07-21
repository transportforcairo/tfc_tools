import os
import sys
# nosec B404 — subprocess is required to install the plugin's pinned
# dependencies into a plugin-local ./libs folder; see _pip_install below, where
# the command is built in-process and run without a shell.
import subprocess  # nosec B404
import traceback
import site
import contextlib

from qgis.PyQt.QtWidgets import QMessageBox

# Keep default empty; we read pinned deps from requirements.txt at the plugin root.
DEFAULT_PKGS = []


def plugin_root_dir():
    # __file__ is this module inside the plugin
    here = os.path.dirname(os.path.abspath(__file__))
    # go one up to plugin root
    return os.path.dirname(here)


def libs_dir():
    return os.path.join(plugin_root_dir(), "libs")


def _ensure_libs_path_first():
    import importlib

    importlib.invalidate_caches()

    ld = libs_dir()
    if ld not in sys.path:
        sys.path.insert(0, ld)
    else:
        # move to front
        sys.path.remove(ld)
        sys.path.insert(0, ld)


def _resolve_pkgs():
    # read requirements.txt if exists; else DEFAULT_PKGS
    req = os.path.join(plugin_root_dir(), "requirements.txt")
    pkgs = []
    if os.path.exists(req):
        for line in open(req, "r", encoding="utf-8").read().splitlines():
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            pkgs.append(s)
    else:
        pkgs = DEFAULT_PKGS[:]
    return pkgs


def _find_python_exe():
    exe = sys.executable
    # On Windows, sys.executable may be qgis-bin.exe; prefer python.exe next to it
    if os.name == "nt" and (exe.lower().endswith("qgis-bin.exe") or exe.lower().endswith("qgis-ltr-bin.exe")):
        cand = os.path.join(os.path.dirname(exe), "python.exe")
        if os.path.exists(cand):
            return cand
    return exe


def _version_in_target(dist_name: str, target_path: str):
    """Return installed distribution version inside target_path, else None."""
    try:
        import importlib.metadata

        for dist in importlib.metadata.distributions(path=[target_path]):
            name = (dist.metadata.get("Name") or "").strip().lower()
            if name == dist_name.strip().lower():
                return dist.version
    except Exception:
        return None
    return None


def ensure_deps(show_ui=True):
    """Install pinned deps into ./libs (plugin-local) on demand.

    Notes:
    - We intentionally do NOT validate against the global QGIS Python environment because
      it may contain unrelated dependency conflicts (pip shows scary warnings).
    - We also suppress the Windows console window for pip to avoid confusing users.
    """
    _ensure_libs_path_first()

    target = libs_dir()
    os.makedirs(target, exist_ok=True)

    # Add target to site dirs so pkg resources work
    site.addsitedir(target)

    pkgs = _resolve_pkgs()
    if not pkgs:
        return

    missing = []
    for spec in pkgs:
        mod = spec.split("==")[0].split(">=")[0].split("[")[0]
        version_pin = None
        if "==" in spec:
            version_pin = spec.split("==")[1]
        elif ">=" in spec:
            version_pin = spec.split(">=")[1]

        installed_ver = _version_in_target(mod, target)
        if installed_ver is None:
            missing.append(spec)
        elif version_pin and installed_ver != version_pin:
            missing.append(spec)

    if not missing:
        return

    py = _find_python_exe()
    cmd = [
        py,
        "-m",
        "pip",
        "install",
        "--disable-pip-version-check",
        "--no-input",
        "--no-warn-script-location",
        "--upgrade",
        "--target",
        target,
    ] + missing

    creationflags = 0
    if os.name == "nt":
        try:
            creationflags = subprocess.CREATE_NO_WINDOW
        except Exception:
            creationflags = 0

    try:
        # nosec B603 — `cmd` is built entirely in-process from sys.executable plus
        # the plugin's own pinned requirement strings (see above); no user input
        # reaches it. It is passed as an argument list with shell=False, so there
        # is no shell interpolation.
        proc = subprocess.run(  # nosec B603
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            creationflags=creationflags,
            check=False,
        )
        if proc.returncode != 0:
            raise RuntimeError(proc.stdout or "pip install failed")
    except Exception as e:
        msg = (
            "Couldn't install required Python packages for TfC Tools.\n\n"
            f"Command: {' '.join(cmd)}\n\n"
            f"Error: {e}\n\n"
            f"Output (if any):\n{getattr(e, 'args', [''])[0]}\n\n{traceback.format_exc()}\n"
            "You can install manually with the same command."
        )
        if show_ui:
            with contextlib.suppress(Exception):
                QMessageBox.critical(None, "TfC Tools: dependency install failed", msg)
        raise


def ensure_paths():
    """Ensure the plugin's ./libs directory is on sys.path (no pip install)."""
    _ensure_libs_path_first()
    target = libs_dir()
    os.makedirs(target, exist_ok=True)
    site.addsitedir(target)


def ensure_bootstrap():
    # Public one-liner: put libs on sys.path and install if needed
    _ensure_libs_path_first()
    ensure_deps(show_ui=True)


# ---------------------------------------------------------------------------
# Runtime library compatibility check
# ---------------------------------------------------------------------------
# Some libraries (pandas, geopandas, shapely, numpy) are NOT pinned in
# requirements.txt because they are heavyweight and typically already provided
# by the QGIS Python environment. The trade-off is that the user may have a
# version older or newer than what TfC Tools has been tested against.
#
# This module reports a single message at plugin load if any installed runtime
# library falls outside the tested range. It never raises; the plugin still
# loads. The check runs at most once per QGIS session.

# (lower_inclusive, upper_exclusive). Set to None on either side to skip that bound.
# Lower bounds are versions known to work. Upper bounds are the next major
# release in which breaking changes are expected.
RUNTIME_LIB_RANGES = {
    "pandas":    ("2.0", "3.0"),
    "geopandas": ("0.13", "2.0"),
    "shapely":   ("2.0", "3.0"),
    "numpy":     ("1.23", "3.0"),
}

_compat_check_done = False


def _parse_version(v):
    """Parse 'X.Y[.Z][...]' to a comparable tuple of ints.

    Tolerant of suffixes like '2.2.3.dev0+abc' or '1.26.4rc1'; trailing
    non-numeric junk is dropped. Returns () on parse failure (which compares
    less than any real version, so an unparseable version is treated as
    "below floor" — the safer side to err on)."""
    if not isinstance(v, str):
        return ()
    parts = []
    for chunk in v.split(".")[:3]:
        digits = ""
        for ch in chunk:
            if ch.isdigit():
                digits += ch
            else:
                break
        if not digits:
            break
        parts.append(int(digits))
    return tuple(parts)


def _get_installed_version(dist_name):
    """Return the installed version string for a top-level package, or None."""
    try:
        import importlib.metadata
        return importlib.metadata.version(dist_name)
    except Exception:
        return None


def check_runtime_lib_compatibility(show_ui=True):
    """Check installed versions of unpinned runtime libs against tested ranges.

    Reports a single QMessageBox if anything is out of range. Always safe to
    call — exceptions are swallowed so a broken check never prevents plugin
    load. Runs at most once per QGIS session (controlled by _compat_check_done).
    """
    global _compat_check_done
    if _compat_check_done:
        return
    _compat_check_done = True

    try:
        too_old = []   # (name, installed, floor)
        too_new = []   # (name, installed, ceiling)
        missing = []   # name

        for name, (lo, hi) in RUNTIME_LIB_RANGES.items():
            installed = _get_installed_version(name)
            if installed is None:
                missing.append(name)
                continue
            iv = _parse_version(installed)
            if not iv:
                # Unparseable; skip rather than guess.
                continue
            if lo is not None and iv < _parse_version(lo):
                too_old.append((name, installed, lo))
            if hi is not None and iv >= _parse_version(hi):
                too_new.append((name, installed, hi))

        if not (too_old or too_new or missing):
            return

        lines = []
        if missing:
            lines.append(
                "Missing libraries (the plugin will not run without these):\n  - "
                + "\n  - ".join(missing)
            )
        if too_old:
            lines.append(
                "Below the tested floor (older than what TfC Tools has been tested with;\n"
                "the plugin may still work but some features could fail):\n  - "
                + "\n  - ".join(f"{n} {v} (tested with >= {floor})" for n, v, floor in too_old)
            )
        if too_new:
            lines.append(
                "Above the tested ceiling (newer than what TfC Tools has been tested with;\n"
                "watch for silently incorrect results — for example, merges returning empty):\n  - "
                + "\n  - ".join(f"{n} {v} (tested with < {ceil})" for n, v, ceil in too_new)
            )
        msg = (
            "TfC Tools detected runtime library versions outside its tested range.\n\n"
            + "\n\n".join(lines)
            + "\n\nThe plugin will continue to load. If you encounter unexpected "
              "behaviour, consider aligning your QGIS Python environment to the "
              "tested versions, or report the issue to TfC."
        )

        # Log to QGIS message panel (always)
        with contextlib.suppress(Exception):
            from qgis.core import QgsMessageLog, Qgis
            QgsMessageLog.logMessage(msg, "TfC Tools", Qgis.MessageLevel.Warning)

        # Modal once per session (only if UI requested and a Qt app exists)
        if show_ui:
            with contextlib.suppress(Exception):
                QMessageBox.warning(None, "TfC Tools: library version check", msg)
    except Exception:
        # Never let the version check itself break plugin load.
        return
