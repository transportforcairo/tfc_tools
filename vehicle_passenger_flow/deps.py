import os, sys, subprocess, site, traceback
from qgis.PyQt.QtWidgets import QMessageBox

PKGS = [
    "sqlalchemy==2.0.43",   # pin versions you rely on
    # "yourlib==x.y.z",
]

def _find_python_exe():
    """
    Returns the actual Python executable bundled with QGIS,
    not qgis-bin.exe.
    """
    # sys.executable may point to qgis-bin.exe on Windows
    exe = sys.executable

    if exe.lower().endswith("qgis-bin.exe"):
        # Windows standalone / OSGeo4W
        exe = os.path.join(os.path.dirname(exe), "..", "apps", "Python39", "python.exe")
        exe = os.path.normpath(exe)
    elif exe.endswith("QGIS"):  # macOS app bundle
        exe = os.path.join(os.path.dirname(exe), "bin", "python3")
    # Linux: sys.executable is usually already /usr/bin/python3

    return exe


def ensure_deps():
    plugin_dir = os.path.dirname(__file__)
    target_dir = os.path.join(plugin_dir, "extlibs")
    os.makedirs(target_dir, exist_ok=True)

    # add target_dir to sys.path so later imports work
    if target_dir not in sys.path:
        sys.path.insert(0, target_dir)

    # quick check: if already installed, skip
    try:
        import sqlalchemy  # replace with a lib you need
        return
    except ImportError:
        pass

    pyexe = _find_python_exe()

    try:
        cmd = [
            pyexe, "-m", "pip", "install",
            "--no-color", "--disable-pip-version-check",
            "--no-cache-dir",
            "--target", target_dir,
        ] + PKGS

        subprocess.check_call(cmd)
    except Exception as e:
        QMessageBox.critical(
            None,
            "Plugin dependency install failed",
            "I couldn't install required Python packages.\n\n"
            f"Command: {' '.join(cmd)}\n\n"
            f"Error:\n{e}\n\n{traceback.format_exc()}\n"
            "You can also install manually with:\n"
            f"\"{pyexe}\" -m pip install --target \"{target_dir}\" " + " ".join(PKGS)
        )
        raise
