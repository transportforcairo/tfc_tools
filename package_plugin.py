#!/usr/bin/env python3
"""
package_plugin.py — build a clean plugin zip for plugins.qgis.org.

Why this exists
---------------
`pb_tool.cfg` is stale (its extra_dirs list is empty, so it would ship a plugin
that cannot import), and the Makefile deploy target is unix-only with a broken
EXTRA_DIRS loop. This script is cross-platform, has no dependencies, and is the
one place that defines what actually ships.

It excludes the QGIS Plugin Builder test scaffolding, dev-only tooling and build
artefacts — the test boilerplate is what produced the Bandit B101 (assert)
findings that the plugin reviewers flagged.

Usage
-----
    python package_plugin.py                # -> tfc_tools-<version>.zip
    python package_plugin.py --out dist/    # choose output directory
    python package_plugin.py --list         # dry run: print what would ship
"""

from __future__ import annotations

import argparse
import configparser
import os
import sys
import zipfile
from pathlib import Path

PLUGIN_NAME = "tfc_tools"
ROOT = Path(__file__).resolve().parent

# Directory names excluded anywhere in the tree.
EXCLUDE_DIRS = {
    "test", "tests", "__pycache__", ".git", ".pytest_cache",
    "scripts",          # dev shell helpers
    ".mypy_cache", ".ruff_cache",
}

# Filenames excluded anywhere in the tree (build/dev tooling ships in each
# sub-plugin folder too, so matching on name rather than path is required).
EXCLUDE_NAMES = {
    "Makefile", "make.bat", "pb_tool.cfg", "pylintrc",
    "package_plugin.py", "plugin_upload.py",
    ".gitignore", ".gitattributes",
}

# Exact repo-relative paths excluded.
EXCLUDE_PATHS = {
    # Dev-only validation harness: a standalone CLI, never imported by the
    # plugin. Keeping it out of the zip also removes its SQL-building code
    # from anything the reviewers scan.
    "validation/compare_stops.py",
    "validation/sousse_baseline.json",
    # Working documents that belong in the repo, not in the plugin.
    "GTFS_Edit_Model_Specification_v1.docx",
    "GTFS_Editing_UserGuide_Section.md",
    "TfC Tools Plugin - User Guide (updated).docx",
}

EXCLUDE_SUFFIXES = {".pyc", ".pyo", ".zip", ".ts"}  # .ts = translation source; .qm ships


def plugin_version() -> str:
    cfg = configparser.ConfigParser()
    cfg.read(ROOT / "metadata.txt", encoding="utf-8")
    return cfg.get("general", "version", fallback="0.0")


def should_skip(rel: Path) -> bool:
    if any(part in EXCLUDE_DIRS for part in rel.parts):
        return True
    if rel.name in EXCLUDE_NAMES:
        return True
    if rel.as_posix() in EXCLUDE_PATHS:
        return True
    if rel.suffix.lower() in EXCLUDE_SUFFIXES:
        return True
    return False


def collect() -> list[Path]:
    out = []
    for path in sorted(ROOT.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(ROOT)
        if should_skip(rel):
            continue
        out.append(rel)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=".", help="output directory (default: current)")
    ap.add_argument("--list", action="store_true", help="dry run; list files and exit")
    args = ap.parse_args(argv)

    files = collect()
    if args.list:
        for f in files:
            print(f.as_posix())
        print(f"\n{len(files)} file(s) would be packaged.", file=sys.stderr)
        return 0

    version = plugin_version()
    out_dir = Path(args.out).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    zip_path = out_dir / f"{PLUGIN_NAME}-{version}.zip"

    # Sanity: refuse to ship a package that is obviously broken.
    required = ["__init__.py", "metadata.txt", "tfc_tools.py"]
    missing = [r for r in required if Path(r) not in files]
    if missing:
        print(f"ERROR: required file(s) missing from package: {missing}", file=sys.stderr)
        return 1

    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
        for rel in files:
            # QGIS requires everything under a single top-level plugin folder.
            zf.write(ROOT / rel, arcname=str(Path(PLUGIN_NAME) / rel))

    size_kb = zip_path.stat().st_size / 1024
    print(f"Built {zip_path.name} — {len(files)} files, {size_kb:.0f} KB")
    print(f"  {zip_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
