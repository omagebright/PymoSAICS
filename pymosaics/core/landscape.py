"""Structural-landscape projection from multi-model PDB trajectories."""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple


LANDMARK_ATOMS = {
    "N", "CA", "C", "O",  # protein backbone
    "P", "O5'", "C5'", "C4'", "C3'", "O3'", "C1'",  # nucleic-acid backbone/sugar
    "O5*", "C5*", "C4*", "C3*", "O3*", "C1*",  # legacy PDB names
}


@dataclass(frozen=True)
class LandscapeResult:
    coordinates: Tuple[Tuple[float, float], ...]
    frame_numbers: Tuple[int, ...]
    representative_frames: Tuple[int, ...]
    rmsd_matrix: Tuple[Tuple[float, ...], ...]
    atom_count: int


def _numpy():
    try:
        import numpy as np
    except ImportError as exc:
        raise RuntimeError(
            "Structural mapping requires NumPy. Use a PyMOL build that includes NumPy."
        ) from exc
    return np


def atom_key(line: str):
    line = line.ljust(80)
    return (line[21].strip(), line[22:26].strip(), line[26].strip(), line[12:16].strip())


def atom_element(line: str) -> str:
    line = line.ljust(80)
    value = line[76:78].strip().upper()
    return value or line[12:16].strip().lstrip("0123456789")[:1].upper()


def iter_pdb_frames(lines):
    current = []
    explicit_models = False
    for raw in lines:
        raw = raw.rstrip("\r\n")
        record = raw[:6].strip().upper()
        if record == "MODEL":
            explicit_models = True
            if current:
                yield current
                current = []
            continue
        if record == "ENDMDL":
            if current:
                yield current
                current = []
            continue
        if record in ("ATOM", "HETATM"):
            current.append(raw)
            continue
        if record == "END" and current and not explicit_models:
            yield current
            current = []
    if current:
        yield current


def read_coordinate_frames(path: Path, maximum_frames: int = 500):
    """Return sampled frames with a consistent landmark-atom ordering."""

    np = _numpy()
    path = path.expanduser().resolve()
    if not path.is_file():
        raise ValueError("trajectory does not exist: {}".format(path))
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        frame_count = sum(1 for _frame in iter_pdb_frames(handle))
    if frame_count < 2:
        raise ValueError("a structural map requires a PDB trajectory with at least two frames")
    if maximum_frames < 2:
        raise ValueError("maximum_frames must be at least two")

    if frame_count > maximum_frames:
        sampled = np.linspace(0, frame_count - 1, maximum_frames, dtype=int).tolist()
        sampled_indices = tuple(dict.fromkeys(int(value) for value in sampled))
    else:
        sampled_indices = tuple(range(frame_count))

    selected_indices = set(sampled_indices)
    raw_frames = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for index, frame in enumerate(iter_pdb_frames(handle)):
            if index in selected_indices:
                raw_frames.append(frame)
    if len(raw_frames) != len(sampled_indices):
        raise ValueError("trajectory changed while it was being read; retry the structural map")

    parsed = []
    for frame in raw_frames:
        atoms = {}
        for raw in frame:
            line = raw.ljust(80)
            if line[16].strip() not in ("", "A") or atom_element(line) == "H":
                continue
            try:
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            atoms[atom_key(line)] = xyz
        parsed.append(atoms)

    common = set(parsed[0])
    for atoms in parsed[1:]:
        common.intersection_update(atoms)
    landmarks = sorted(key for key in common if key[3].upper() in LANDMARK_ATOMS)
    keys = landmarks if len(landmarks) >= 3 else sorted(common)
    if len(keys) < 3:
        raise ValueError("trajectory frames do not share at least three usable atoms")
    coordinates = np.asarray([[atoms[key] for key in keys] for atoms in parsed], dtype=float)
    return coordinates, tuple(index + 1 for index in sampled_indices)


def kabsch_rmsd(first, second) -> float:
    np = _numpy()
    first = np.asarray(first, dtype=float)
    second = np.asarray(second, dtype=float)
    if first.shape != second.shape or first.ndim != 2 or first.shape[1] != 3:
        raise ValueError("coordinate arrays must have matching (atoms, 3) shapes")
    mobile = first - first.mean(axis=0)
    reference = second - second.mean(axis=0)
    covariance = mobile.T @ reference
    left, _, right_t = np.linalg.svd(covariance)
    correction = np.eye(3)
    correction[-1, -1] = 1.0 if np.linalg.det(left @ right_t) >= 0 else -1.0
    rotation = left @ correction @ right_t
    difference = mobile @ rotation - reference
    return float(np.sqrt(np.mean(np.sum(difference * difference, axis=1))))


def pairwise_rmsd(frames):
    np = _numpy()
    frame_count = len(frames)
    matrix = np.zeros((frame_count, frame_count), dtype=float)
    for first in range(frame_count):
        for second in range(first + 1, frame_count):
            value = kabsch_rmsd(frames[first], frames[second])
            matrix[first, second] = value
            matrix[second, first] = value
    return matrix


def classical_mds(distance_matrix):
    np = _numpy()
    distances = np.asarray(distance_matrix, dtype=float)
    count = distances.shape[0]
    if distances.shape != (count, count) or count < 2:
        raise ValueError("distance matrix must be square with at least two frames")
    centering = np.eye(count) - np.ones((count, count)) / count
    gram = -0.5 * centering @ (distances * distances) @ centering
    values, vectors = np.linalg.eigh(gram)
    order = np.argsort(values)[::-1]
    positive = [(values[index], vectors[:, index]) for index in order if values[index] > 1e-12]
    result = np.zeros((count, 2), dtype=float)
    for axis, (value, vector) in enumerate(positive[:2]):
        result[:, axis] = vector * np.sqrt(value)
    return result


def select_representatives(distance_matrix, count: int):
    """Select deterministic farthest-point medoids from the sampled frames."""

    np = _numpy()
    distances = np.asarray(distance_matrix, dtype=float)
    frame_count = distances.shape[0]
    count = max(1, min(int(count), frame_count))
    selected = [int(np.argmin(distances.mean(axis=1)))]
    while len(selected) < count:
        nearest = distances[:, selected].min(axis=1)
        nearest[selected] = -1.0
        selected.append(int(np.argmax(nearest)))
    return tuple(selected)


def build_landscape(path: Path, representatives: int = 6, maximum_frames: int = 500) -> LandscapeResult:
    frames, frame_numbers = read_coordinate_frames(path, maximum_frames=maximum_frames)
    matrix = pairwise_rmsd(frames)
    coordinates = classical_mds(matrix)
    representative_indices = select_representatives(matrix, representatives)
    return LandscapeResult(
        coordinates=tuple((float(x), float(y)) for x, y in coordinates),
        frame_numbers=frame_numbers,
        representative_frames=tuple(frame_numbers[index] for index in representative_indices),
        rmsd_matrix=tuple(tuple(float(value) for value in row) for row in matrix),
        atom_count=int(frames.shape[1]),
    )


def write_landscape_table(
    path: Path, result: LandscapeResult, energies: Optional[Tuple[float, ...]] = None
) -> Path:
    """Write the exact frame-to-map correspondence for inspection and reuse."""

    representatives = set(result.representative_frames)
    if energies is not None and len(energies) != len(result.coordinates):
        raise ValueError("energy values must match the mapped frame count")
    lines = ["frame\tcoordinate_1\tcoordinate_2\trepresentative\tenergy"]
    for index, (frame, (first, second)) in enumerate(zip(result.frame_numbers, result.coordinates)):
        lines.append(
            "{}\t{:.10g}\t{:.10g}\t{}\t{}".format(
                frame,
                first,
                second,
                "yes" if frame in representatives else "no",
                "" if energies is None else "{:.10g}".format(energies[index]),
            )
        )
    path = path.expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
