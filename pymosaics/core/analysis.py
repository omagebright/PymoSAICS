"""Post-run discovery and lightweight MOSAICS analysis without extra packages."""

import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple


@dataclass(frozen=True)
class EnergySeries:
    path: Path
    label: str
    values: Tuple[float, ...]

    @property
    def minimum(self) -> float:
        return min(self.values)

    @property
    def maximum(self) -> float:
        return max(self.values)

    @property
    def mean(self) -> float:
        return sum(self.values) / len(self.values)


@dataclass(frozen=True)
class AcceptanceSummary:
    chain: int
    accepted: int
    attempted: int
    ratio: float


@dataclass(frozen=True)
class ProjectFile:
    path: Path
    kind: str
    text_readable: bool
    loadable_in_pymol: bool


ENERGY_PATTERNS = ("*.pot_energy", "*.inter_energy", "*potential_energy*.dat", "*inter_energy*.dat")
CHAIN_ACCEPTANCE = re.compile(r"^\s*Chain\s+(\d+)\s+(\d+)\s+(\d+)\s+([0-9]*\.?[0-9]+)\s*$")
TEXT_SUFFIXES = {
    ".txt", ".log", ".input", ".inp", ".dat", ".pdb", ".tsv", ".csv", ".rtf",
    ".bond", ".bend", ".onfo", ".vdw", ".tors", ".impr", ".tors_and_impr",
    ".json", ".md",
}
HIDDEN_PROJECT_PARTS = {"resolved"}


def read_energy_series(path: Path) -> EnergySeries:
    values: List[float] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            fields = raw.split()
            if not fields:
                continue
            try:
                number = float(fields[0])
            except ValueError:
                continue
            if math.isfinite(number):
                values.append(number)
    if not values:
        raise ValueError("no numeric energy values found in {}".format(path))
    label = path.name.replace("_", " ")
    return EnergySeries(path.resolve(), label, tuple(values))


def discover_energy_series(project: Path) -> Tuple[EnergySeries, ...]:
    found = []
    seen = set()
    for pattern in ENERGY_PATTERNS:
        for path in project.expanduser().resolve().rglob(pattern):
            if path in seen or ".pymosaics" in path.parts:
                continue
            seen.add(path)
            try:
                found.append(read_energy_series(path))
            except (OSError, ValueError):
                continue
    found.sort(key=lambda item: item.path.stat().st_mtime, reverse=True)
    return tuple(found)


def parse_acceptance_log(path: Path) -> Tuple[AcceptanceSummary, ...]:
    summaries = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            match = CHAIN_ACCEPTANCE.match(raw)
            if match:
                summaries.append(
                    AcceptanceSummary(
                        int(match.group(1)),
                        int(match.group(2)),
                        int(match.group(3)),
                        float(match.group(4)),
                    )
                )
    return tuple(summaries)


def latest_log(project: Path) -> Optional[Path]:
    logs = list(project.expanduser().resolve().rglob("*.log"))
    logs.extend((project.expanduser().resolve() / ".pymosaics" / "logs").glob("*.log"))
    logs = [path for path in set(logs) if path.is_file()]
    return max(logs, key=lambda path: path.stat().st_mtime) if logs else None


def discover_pdb_outputs(project: Path, exclude: Optional[Path] = None) -> Tuple[Path, ...]:
    root = project.expanduser().resolve()
    excluded = {root / "structure.pdb"}
    if exclude is not None and exclude.exists():
        excluded.add(exclude.resolve())
    paths = []
    for path in root.rglob("*.pdb"):
        resolved = path.resolve()
        if not path.is_file() or ".pymosaics" in path.parts or resolved in excluded:
            continue
        paths.append(resolved)
    paths.sort(key=lambda item: item.stat().st_mtime, reverse=True)
    return tuple(paths)


def classify_project_file(path: Path) -> ProjectFile:
    suffix = path.suffix.lower()
    name = path.name.lower()
    if suffix == ".pdb":
        kind = "PDB structure / trajectory"
    elif suffix == ".log":
        kind = "Run log"
    elif suffix in (".input", ".inp"):
        kind = "MOSAICS input"
    elif "energy" in name:
        kind = "Energy output"
    elif suffix in (".dat", ".tsv", ".csv"):
        kind = "Data output"
    else:
        kind = "Text file" if suffix in TEXT_SUFFIXES else "Output file"
    return ProjectFile(
        path=path.resolve(),
        kind=kind,
        text_readable=suffix in TEXT_SUFFIXES,
        loadable_in_pymol=suffix == ".pdb",
    )


def discover_project_files(project: Path) -> Tuple[ProjectFile, ...]:
    """Return user-facing project files, including PymoSAICS run logs."""

    root = project.expanduser().resolve()
    if not root.is_dir():
        return ()
    files = []
    for path in root.rglob("*"):
        if not path.is_file() or path.name.startswith("."):
            continue
        relative = path.relative_to(root)
        if ".pymosaics" in relative.parts and "logs" not in relative.parts:
            continue
        if any(part in HIDDEN_PROJECT_PARTS for part in relative.parts):
            continue
        files.append(classify_project_file(path))
    files.sort(key=lambda item: item.path.stat().st_mtime, reverse=True)
    return tuple(files)


def read_text_file(path: Path, maximum_bytes: int = 10 * 1024 * 1024) -> str:
    path = path.expanduser().resolve()
    if not path.is_file():
        raise ValueError("file does not exist: {}".format(path))
    if path.stat().st_size > maximum_bytes:
        raise ValueError("file is larger than the 10 MB text-viewer limit: {}".format(path.name))
    data = path.read_bytes()
    if b"\x00" in data[:4096]:
        raise ValueError("file appears to be binary: {}".format(path.name))
    return data.decode("utf-8", errors="replace")
