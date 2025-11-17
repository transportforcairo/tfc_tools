
import os, sys, subprocess, traceback, site
from qgis.PyQt.QtWidgets import QMessageBox

DEFAULT_PKGS = [
    
]

def plugin_root_dir():
    # __file__ is this module inside the plugin
    here = os.path.dirname(os.path.abspath(__file__))
    # go one up to plugin root
    return os.path.dirname(here)

def libs_dir():
    return os.path.join(plugin_root_dir(), "libs")

def _ensure_libs_path_first():
    import importlib

    # If these were imported earlier from core site-packages, drop them
    # for name in ("pandas", "numpy"):
    #     for k in list(sys.modules.keys()):
    #         if k == name or k.startswith(name + "."):
    #             sys.modules.pop(k, None)

    # Make sure Python sees the new lib directory
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
    if os.name == "nt" and ((exe.lower().endswith("qgis-bin.exe")) or (exe.lower().endswith("qgis-ltr-bin.exe"))):
        cand = os.path.join(os.path.dirname(exe), "python.exe")
        if os.path.exists(cand):
            return cand
    return exe

def ensure_deps(show_ui=True):
    _ensure_libs_path_first()
    target = libs_dir()
    os.makedirs(target, exist_ok=True)

    # Add target to site dirs so pkg resources work
    site.addsitedir(target)

    pkgs = _resolve_pkgs()
    if not pkgs:
        return

    # Check if imports already succeed to avoid reinstalling every load
    missing = []
    import importlib.metadata

    for spec in pkgs:
        mod = spec.split("==")[0].split(">=")[0].split("[")[0]
        version_pin = None
        if "==" in spec:
            version_pin = spec.split("==")[1]
        elif ">=" in spec:
            version_pin = spec.split(">=")[1]

        need_install = False
        try:
            installed_ver = importlib.metadata.version(mod)
            if version_pin and installed_ver != version_pin:
                need_install = True
        except importlib.metadata.PackageNotFoundError:
            need_install = True

        if need_install:
            missing.append(spec)


    if not missing:
        return

    py = _find_python_exe()
    cmd = [py, "-m", "pip", "install", "--no-warn-script-location", "--target", target] + missing

    try:
        subprocess.check_call(cmd)
    except Exception as e:
        msg = (
            "Couldn't install required Python packages for TfC Tools.\n\n"
            f"Command: {' '.join(cmd)}\n\n"
            f"Error: {e}\n\n{traceback.format_exc()}\n"
            "You can install manually with the same command."
        )
        if show_ui:
            try:
                QMessageBox.critical(None, "TfC Tools: dependency install failed", msg)
            except Exception:
                pass
        # Re-raise so calling code can decide
        raise

def ensure_bootstrap():
    # Public one-liner: put libs on sys.path and install if needed
    _ensure_libs_path_first()
    ensure_deps(show_ui=True)
