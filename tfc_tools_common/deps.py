import os
import sys
import subprocess
import traceback
import site

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
        proc = subprocess.run(
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
            try:
                QMessageBox.critical(None, "TfC Tools: dependency install failed", msg)
            except Exception:
                pass
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
