"""Explicit MOSAICS region-file generation."""

from dataclasses import dataclass
import math
from pathlib import Path
from typing import Sequence, Tuple


@dataclass(frozen=True)
class RegionSettings:
    residues: Tuple[str, ...]
    centers: Tuple[str, ...]
    residue_pairs: Tuple[Tuple[str, str], ...]
    dependency_type: str = "independent"
    translation_sigma: str = "0.0"
    rotation_sigma: str = "0.0"
    free_translation_sigma: str = ".5e-5"
    free_rotation_sigma: str = ".5e-6"
    pair_translation_sigma: str = ".5e-5"
    pair_rotation_sigma: str = ".5e-6"


def generate_region_file(settings: RegionSettings) -> str:
    if not settings.residues:
        raise ValueError("a region must contain at least one residue")
    if settings.dependency_type != "independent":
        raise ValueError("residue-level regions must use the supported independent dependency type")
    residue_set = set(settings.residues)
    unknown_centers = [item for item in settings.centers if item not in residue_set]
    unknown_pairs = [item for pair in settings.residue_pairs for item in pair if item not in residue_set]
    if unknown_centers or unknown_pairs:
        raise ValueError("centers and residue pairs must belong to the selected region")
    if not settings.centers:
        raise ValueError("a residue-level region must define at least one rotation center")
    if len(set(settings.residues)) != len(settings.residues):
        raise ValueError("a residue may appear only once in a region")
    normalized_pairs = {tuple(sorted(pair)) for pair in settings.residue_pairs}
    if any(first == second for first, second in settings.residue_pairs):
        raise ValueError("a residue cannot be paired with itself")
    if len(normalized_pairs) != len(settings.residue_pairs):
        raise ValueError("each residue pair may appear only once")
    paired_residues = [item for pair in settings.residue_pairs for item in pair]
    if len(set(paired_residues)) != len(paired_residues):
        raise ValueError("a residue may belong to only one residue pair")
    sigma_values = (
        settings.translation_sigma,
        settings.rotation_sigma,
        settings.free_translation_sigma,
        settings.free_rotation_sigma,
        settings.pair_translation_sigma,
        settings.pair_rotation_sigma,
    )
    try:
        parsed_sigmas = tuple(float(value) for value in sigma_values)
    except (TypeError, ValueError) as exc:
        raise ValueError("all region proposal widths must be numeric") from exc
    if any(value < 0 or not math.isfinite(value) for value in parsed_sigmas):
        raise ValueError("all region proposal widths must be non-negative finite numbers")

    pair_lines = "".join(
        "        \\residue_pair{{{},{}}}\n".format(first, second)
        for first, second in settings.residue_pairs
    )
    return (
        "__________________________region____________________________\n"
        "~region[\\element_top_type{{residue}}\n"
        "        \\dependency_type{{{}}}\n\n"
        "        \\nres{{{}}}\n"
        "        \\residues{{{}}}\n\n"
        "        \\ncenter{{{}}}\n"
        "        \\centers{{{}}}\n\n"
        "        \\nrespair{{{}}}\n"
        "{}\n"
        "        \\prop_trans_sig{{{}}}\n"
        "        \\prop_rot_sig{{{}}}\n"
        "        \\prop_trans_sig_freeres{{{}}}\n"
        "        \\prop_rot_sig_freeres{{{}}}\n"
        "        \\prop_trans_sig_respair{{{}}}\n"
        "        \\prop_rot_sig_respair{{{}}}\n"
        "]\n"
    ).format(
        settings.dependency_type,
        len(settings.residues),
        ",".join(settings.residues),
        len(settings.centers),
        ",".join(settings.centers),
        len(settings.residue_pairs),
        pair_lines.rstrip(),
        settings.translation_sigma,
        settings.rotation_sigma,
        settings.free_translation_sigma,
        settings.free_rotation_sigma,
        settings.pair_translation_sigma,
        settings.pair_rotation_sigma,
    )


def write_region_file(path: Path, settings: RegionSettings) -> Path:
    path = path.expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(generate_region_file(settings), encoding="utf-8")
    return path
