"""Aligned structural, sugar-pucker, and terminal-mobility trajectory analysis."""

import math
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

from .landscape import LANDMARK_ATOMS, atom_element, atom_key, iter_pdb_frames


@dataclass(frozen=True)
class PuckerMeasurement:
    frame: int
    chain: str
    residue: str
    residue_name: str
    phase_degrees: float
    state: str


@dataclass(frozen=True)
class ResidueChange:
    chain: str
    residue: str
    residue_name: str
    rmsd_angstrom: float


@dataclass(frozen=True)
class TerminalMobility:
    chain: str
    end: str
    residue: str
    residue_name: str
    atom_count: int
    invariant_atom_count: int
    maximum_displacement: float
    mean_maximum_displacement: float


@dataclass(frozen=True)
class TrajectoryAnalysis:
    path: Path
    frame_numbers: Tuple[int, ...]
    rmsd_to_first: Tuple[float, ...]
    puckers: Tuple[PuckerMeasurement, ...]
    residue_changes: Tuple[ResidueChange, ...]
    terminal_mobility: Tuple[TerminalMobility, ...]
    aligned_atom_count: int

    @property
    def start_to_end_rmsd(self) -> float:
        return self.rmsd_to_first[-1]

    @property
    def maximum_rmsd(self) -> float:
        return max(self.rmsd_to_first)


def _numpy():
    try:
        import numpy as np
    except ImportError as exc:
        raise RuntimeError("Trajectory analysis requires NumPy from the PyMOL distribution.") from exc
    return np


def _dihedral(first, second, third, fourth) -> float:
    np = _numpy()
    first, second, third, fourth = (
        np.asarray(value, dtype=float) for value in (first, second, third, fourth)
    )
    first_bond = -(second - first)
    axis = third - second
    last_bond = fourth - third
    length = np.linalg.norm(axis)
    if length == 0:
        raise ValueError("cannot calculate a torsion around a zero-length bond")
    axis /= length
    first_plane = first_bond - np.dot(first_bond, axis) * axis
    last_plane = last_bond - np.dot(last_bond, axis) * axis
    x_value = np.dot(first_plane, last_plane)
    y_value = np.dot(np.cross(axis, first_plane), last_plane)
    return math.degrees(math.atan2(y_value, x_value))


def glycosidic_chi(atoms: Dict[str, Sequence[float]]) -> float:
    """Return a purine or pyrimidine glycosidic chi angle in degrees.

    The IUPAC atom order is O4'-C1'-N9-C4 for purines and
    O4'-C1'-N1-C2 for pyrimidines.  Legacy ``*`` sugar atom names are
    accepted so the measurement also works on older MOSAICS output.
    """

    normalized = {name.replace("*", "'").upper(): xyz for name, xyz in atoms.items()}
    if all(name in normalized for name in ("O4'", "C1'", "N9", "C4")):
        required = ("O4'", "C1'", "N9", "C4")
    elif all(name in normalized for name in ("O4'", "C1'", "N1", "C2")):
        required = ("O4'", "C1'", "N1", "C2")
    else:
        raise ValueError("residue is missing the atoms required for a glycosidic chi angle")
    return _dihedral(*(normalized[name] for name in required))


def pseudorotation_phase(atoms: Dict[str, Sequence[float]]) -> float:
    """Return the Altona–Sundaralingam furanose phase angle in degrees."""

    normalized = {name.replace("*", "'").upper(): xyz for name, xyz in atoms.items()}
    required = ("O4'", "C1'", "C2'", "C3'", "C4'")
    missing = [name for name in required if name not in normalized]
    if missing:
        raise ValueError("sugar ring is missing {}".format(", ".join(missing)))
    oxygen, c1, c2, c3, c4 = (normalized[name] for name in required)
    # Altona--Sundaralingam numbers the five endocyclic torsions from
    # C4'-O4'-C1'-C2' around the ring.  Keeping this order is important:
    # starting at O4' instead rotates the reported phase by roughly 144 degrees
    # and can therefore exchange an A-like and B-like assignment.
    torsions = (
        _dihedral(c4, oxygen, c1, c2),
        _dihedral(oxygen, c1, c2, c3),
        _dihedral(c1, c2, c3, c4),
        _dihedral(c2, c3, c4, oxygen),
        _dihedral(c3, c4, oxygen, c1),
    )
    numerator = (torsions[4] + torsions[1]) - (torsions[3] + torsions[0])
    denominator = 2.0 * torsions[2] * (
        math.sin(math.radians(36.0)) + math.sin(math.radians(72.0))
    )
    return math.degrees(math.atan2(numerator, denominator)) % 360.0


def _circular_distance(first: float, second: float) -> float:
    return abs((first - second + 180.0) % 360.0 - 180.0)


def classify_pucker(phase_degrees: float) -> str:
    if _circular_distance(phase_degrees, 18.0) <= 45.0:
        return "A-like / C3'-endo"
    if _circular_distance(phase_degrees, 162.0) <= 45.0:
        return "B-like / C2'-endo"
    return "other"


def _is_nucleic_residue_name(name: str) -> bool:
    normalized = name.strip().upper()
    return normalized in {
        "A", "C", "G", "T", "U", "DA", "DC", "DG", "DT", "DU",
        "RA", "RC", "RG", "RU", "A3", "A5", "C3", "C5", "G3", "G5",
        "T3", "T5", "U3", "U5", "DA3", "DA5", "DC3", "DC5", "DG3",
        "DG5", "DT3", "DT5", "RA3", "RA5", "RC3", "RC5", "RG3",
        "RG5", "RU3", "RU5",
    }


def _read_atom_frames(path: Path, maximum_frames: int):
    path = path.expanduser().resolve()
    if not path.is_file():
        raise ValueError("trajectory does not exist: {}".format(path))
    if maximum_frames < 2:
        raise ValueError("maximum_frames must be at least two")
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        frame_count = sum(1 for _frame in iter_pdb_frames(handle))
    if frame_count < 2:
        raise ValueError("trajectory analysis requires at least two PDB frames")
    np = _numpy()
    if frame_count > maximum_frames:
        sampled = np.linspace(0, frame_count - 1, maximum_frames, dtype=int).tolist()
        indices = tuple(dict.fromkeys(int(value) for value in sampled))
    else:
        indices = tuple(range(frame_count))
    selected = set(indices)
    frames = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for index, frame in enumerate(iter_pdb_frames(handle)):
            if index not in selected:
                continue
            atoms = OrderedDict()
            residue_names = {}
            for raw in frame:
                line = raw.ljust(80)
                if line[16].strip() not in ("", "A") or atom_element(line) == "H":
                    continue
                try:
                    xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                except ValueError:
                    continue
                key = atom_key(line)
                atoms[key] = xyz
                residue_names[key[:3]] = line[17:20].strip()
            frames.append((atoms, residue_names))
    return tuple(frames), tuple(index + 1 for index in indices)


def analyze_structure_puckers(path: Path) -> Tuple[PuckerMeasurement, ...]:
    """Measure furanose pseudorotation in the first structure/frame of a PDB."""

    path = path.expanduser().resolve()
    if not path.is_file():
        raise ValueError("structure does not exist: {}".format(path))
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        try:
            raw_frame = next(iter_pdb_frames(handle))
        except StopIteration:
            return ()
    grouped: Dict[Tuple[str, str, str], Dict[str, Sequence[float]]] = OrderedDict()
    residue_names = {}
    for raw in raw_frame:
        line = raw.ljust(80)
        if line[16].strip() not in ("", "A"):
            continue
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        key = atom_key(line)
        grouped.setdefault(key[:3], {})[key[3]] = xyz
        residue_names[key[:3]] = line[17:20].strip()
    measurements = []
    for residue_key, residue_atoms in grouped.items():
        try:
            phase = pseudorotation_phase(residue_atoms)
        except ValueError:
            continue
        measurements.append(
            PuckerMeasurement(
                1,
                residue_key[0] or "_",
                residue_key[1] + residue_key[2],
                residue_names.get(residue_key, ""),
                phase,
                classify_pucker(phase),
            )
        )
    return tuple(measurements)


def _alignment(mobile, reference):
    np = _numpy()
    mobile = np.asarray(mobile, dtype=float)
    reference = np.asarray(reference, dtype=float)
    mobile_center = mobile.mean(axis=0)
    reference_center = reference.mean(axis=0)
    covariance = (mobile - mobile_center).T @ (reference - reference_center)
    left, _, right_t = np.linalg.svd(covariance)
    correction = np.eye(3)
    correction[-1, -1] = 1.0 if np.linalg.det(left @ right_t) >= 0 else -1.0
    rotation = left @ correction @ right_t
    return rotation, mobile_center, reference_center


def analyze_trajectory(
    path: Path,
    maximum_frames: int = 500,
    invariant_tolerance: float = 0.01,
) -> TrajectoryAnalysis:
    frames, frame_numbers = _read_atom_frames(path, maximum_frames)
    common = set(frames[0][0])
    for atoms, _names in frames[1:]:
        common.intersection_update(atoms)
    common_keys = tuple(key for key in frames[0][0] if key in common)
    landmark_keys = tuple(key for key in common_keys if key[3].upper() in LANDMARK_ATOMS)
    alignment_keys = landmark_keys if len(landmark_keys) >= 3 else common_keys
    if len(alignment_keys) < 3:
        raise ValueError("trajectory frames do not share at least three alignment atoms")

    np = _numpy()
    reference_alignment = np.asarray([frames[0][0][key] for key in alignment_keys], dtype=float)
    common_index = {key: index for index, key in enumerate(common_keys)}
    alignment_indices = [common_index[key] for key in alignment_keys]
    aligned_frames = []
    rmsd_values = []
    for atoms, _names in frames:
        mobile_alignment = np.asarray([atoms[key] for key in alignment_keys], dtype=float)
        rotation, mobile_center, reference_center = _alignment(mobile_alignment, reference_alignment)
        all_coordinates = np.asarray([atoms[key] for key in common_keys], dtype=float)
        aligned = (all_coordinates - mobile_center) @ rotation + reference_center
        aligned_frames.append(aligned)
        difference = aligned[alignment_indices] - reference_alignment
        rmsd_values.append(float(np.sqrt(np.mean(np.sum(difference * difference, axis=1)))))

    puckers: List[PuckerMeasurement] = []
    for frame_number, (atoms, residue_names) in zip(frame_numbers, frames):
        grouped: Dict[Tuple[str, str, str], Dict[str, Sequence[float]]] = OrderedDict()
        for key, xyz in atoms.items():
            grouped.setdefault(key[:3], {})[key[3]] = xyz
        for residue_key, residue_atoms in grouped.items():
            try:
                phase = pseudorotation_phase(residue_atoms)
            except ValueError:
                continue
            puckers.append(
                PuckerMeasurement(
                    frame_number,
                    residue_key[0] or "_",
                    residue_key[1] + residue_key[2],
                    residue_names.get(residue_key, ""),
                    phase,
                    classify_pucker(phase),
                )
            )

    first_names = frames[0][1]
    final_aligned = aligned_frames[-1]
    initial_aligned = aligned_frames[0]
    residue_indices: Dict[Tuple[str, str, str], List[int]] = OrderedDict()
    for index, key in enumerate(common_keys):
        residue_indices.setdefault(key[:3], []).append(index)
    changes = []
    for residue_key, indices in residue_indices.items():
        difference = final_aligned[indices] - initial_aligned[indices]
        value = float(np.sqrt(np.mean(np.sum(difference * difference, axis=1))))
        changes.append(
            ResidueChange(
                residue_key[0] or "_",
                residue_key[1] + residue_key[2],
                first_names.get(residue_key, ""),
                value,
            )
        )
    changes.sort(key=lambda item: item.rmsd_angstrom, reverse=True)

    chain_residues: Dict[str, List[Tuple[str, str, str]]] = OrderedDict()
    for residue_key in residue_indices:
        chain_residues.setdefault(residue_key[0], []).append(residue_key)
    terminal = []
    for chain, residues in chain_residues.items():
        nucleic = any(_is_nucleic_residue_name(first_names.get(key, "")) for key in residues)
        if not nucleic:
            for residue_key in residues:
                atom_names = {
                    key[3].replace("*", "'").upper()
                    for key in common_keys
                    if key[:3] == residue_key
                }
                if {"O4'", "C1'", "C2'", "C3'", "C4'"}.issubset(atom_names):
                    nucleic = True
                    break
        labels = ("5'", "3'") if nucleic else ("N", "C")
        for label, residue_key in ((labels[0], residues[0]), (labels[1], residues[-1])):
            indices = residue_indices[residue_key]
            maxima = []
            for index in indices:
                maxima.append(
                    max(
                        float(np.linalg.norm(aligned[index] - initial_aligned[index]))
                        for aligned in aligned_frames
                    )
                )
            terminal.append(
                TerminalMobility(
                    chain or "_",
                    label,
                    residue_key[1] + residue_key[2],
                    first_names.get(residue_key, ""),
                    len(indices),
                    sum(value <= invariant_tolerance for value in maxima),
                    max(maxima) if maxima else 0.0,
                    sum(maxima) / len(maxima) if maxima else 0.0,
                )
            )

    return TrajectoryAnalysis(
        path.expanduser().resolve(),
        frame_numbers,
        tuple(rmsd_values),
        tuple(puckers),
        tuple(changes),
        tuple(terminal),
        len(alignment_keys),
    )
