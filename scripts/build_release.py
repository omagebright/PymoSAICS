#!/usr/bin/env python3
"""Build the PyMOL Plugin Manager ZIP using only the standard library."""

import argparse
import zipfile
from pathlib import Path


VERSION = "0.2.3"
ROOT = Path(__file__).resolve().parents[1]
INCLUDED_ROOT_FILES = ("README.md", "LICENSE", "CHANGELOG.md")


def included_files():
    for name in INCLUDED_ROOT_FILES:
        yield ROOT / name
    for path in sorted((ROOT / "pymosaics").rglob("*")):
        if path.is_file() and "__pycache__" not in path.parts and path.suffix not in (".pyc", ".pyo"):
            yield path


def build(destination: Path) -> Path:
    destination.mkdir(parents=True, exist_ok=True)
    archive = destination / "PymoSAICS-{}.zip".format(VERSION)
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as bundle:
        for path in included_files():
            relative = path.relative_to(ROOT).as_posix()
            info = zipfile.ZipInfo(relative, date_time=(2026, 7, 22, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.create_system = 3
            mode = 0o755 if relative.startswith("pymosaics/assets/runtimes/") else 0o644
            info.external_attr = mode << 16
            bundle.writestr(info, path.read_bytes())
    return archive


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "dist")
    args = parser.parse_args()
    print(build(args.output))


if __name__ == "__main__":
    main()
