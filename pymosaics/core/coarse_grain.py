"""Transparent conversion of canonical proteins to MOSAICS KB_3pt models."""

import csv
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

from .structures import PreparedStructure


BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}
PROTEIN_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}
OUTPUT_CHAIN_IDS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
RESIDUE_ALIASES = {
    "CYX": "CYS",
    "HID": "HIS",
    "HIE": "HIS",
    "HIP": "HIS",
    "HSD": "HIS",
    "HSE": "HIS",
    "HSP": "HIS",
    "MSE": "MET",
}


def _model_number(line: str, fallback: int) -> str:
    value = line[10:14].strip() if len(line) >= 14 else ""
    return value or str(fallback)


def _selected_residues(text: str, model: str, chains: Sequence[str]):
    selected = set(chains)
    residues = OrderedDict()
    current_model = "1"
    model_counter = 0
    for raw in text.splitlines():
        record = raw[:6].strip().upper()
        if record == "MODEL":
            model_counter += 1
            current_model = _model_number(raw, model_counter)
            continue
        if record not in ("ATOM", "HETATM") or current_model != model:
            continue
        line = raw.ljust(80)
        chain = line[21].strip()
        if chain not in selected or line[16].strip() not in ("", "A"):
            continue
        key = (chain, line[22:26].strip(), line[26].strip())
        name = line[12:16].strip().upper()
        residue_name = line[17:20].strip().upper()
        normalized_residue = RESIDUE_ALIASES.get(residue_name, residue_name)
        if record == "HETATM" and normalized_residue not in PROTEIN_RESIDUES:
            continue
        element = line[76:78].strip().upper()
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError as exc:
            raise ValueError("invalid coordinates for {}:{} {}".format(chain or "_", key[1], name)) from exc
        if key not in residues:
            residues[key] = {"name": residue_name, "atoms": OrderedDict()}
        residues[key]["atoms"].setdefault(name, (xyz, element))
    if not residues:
        raise ValueError("the selected model and chains contain no protein atoms")
    return residues


def _centroid(points: Sequence[Tuple[float, float, float]]) -> Tuple[float, float, float]:
    return tuple(sum(point[index] for point in points) / len(points) for index in range(3))


def _atom_line(serial: int, atom: str, residue: str, chain: str, number: int, xyz, element: str) -> str:
    if number > 9999:
        raise ValueError("three-point residue number exceeds the PDB format limit")
    if serial > 99999:
        raise ValueError("three-point atom serial exceeds the PDB format limit")
    return (
        "ATOM  {serial:5d} {atom:>4s} {residue:>3s} {chain:1s}{number:4d}    "
        "{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}"
    ).format(
        serial=serial,
        atom=atom,
        residue=residue,
        chain=chain or "A",
        number=number,
        x=xyz[0],
        y=xyz[1],
        z=xyz[2],
        element=element,
    )


def prepare_three_point_structure(
    source: Path,
    destination: Path,
    model: str,
    chains: Sequence[str],
    closure_residues: Sequence[str] = (),
) -> PreparedStructure:
    """Write CA, carbonyl O, and geometric side-chain centroid per residue.

    Residues are renumbered from one within each chain and a visible mapping is
    emitted. ``STRIDE`` defaults to one rigid segment per chain. Callers may
    identify loop/closure residues using original selectors such as ``A:42``.
    Glycine has no side chain, so its CMA pseudoatom is placed 0.01 Å from CA,
    matching the non-coincident convention used by historical MOSAICS decks.
    """

    source = source.expanduser().resolve()
    residues = _selected_residues(
        source.read_text(encoding="utf-8", errors="replace"), model, chains
    )
    closure_set = set(closure_residues)
    per_chain: Dict[str, List[Tuple[Tuple[str, str, str], dict]]] = OrderedDict()
    for key, residue in residues.items():
        per_chain.setdefault(key[0], []).append((key, residue))

    serial = 1
    atom_lines: List[str] = []
    mapping_rows = []
    emitted_residues = 0
    stride_lines = []
    if len(per_chain) > len(OUTPUT_CHAIN_IDS):
        raise ValueError("three-point conversion supports at most 62 chains")
    output_chains = {
        input_chain: OUTPUT_CHAIN_IDS[index]
        for index, input_chain in enumerate(per_chain)
    }
    for chain, entries in per_chain.items():
        output_chain = output_chains[chain]
        stride = []
        output_number = 0
        for key, residue in entries:
            original_name = residue["name"]
            residue_name = RESIDUE_ALIASES.get(original_name, original_name)
            if residue_name not in PROTEIN_RESIDUES:
                raise ValueError(
                    "unsupported three-point residue {} at {}:{}".format(
                        original_name, chain or "_", key[1]
                    )
                )
            atoms = residue["atoms"]
            missing = [name for name in ("CA", "O") if name not in atoms]
            if missing:
                input_selector = "{}:{}{}".format(chain or "_", key[1], key[2])
                mapping_rows.append(
                    {
                        "input_selector": input_selector,
                        "output_selector": "",
                        "input_residue": original_name,
                        "output_residue": residue_name,
                        "cma_method": "omitted: missing {} required for KB_3pt".format(", ".join(missing)),
                        "status": "omitted incomplete residue",
                    }
                )
                continue
            output_number += 1
            emitted_residues += 1
            ca = atoms["CA"][0]
            oxygen = atoms["O"][0]
            side_chain = [
                xyz
                for atom_name, (xyz, element) in atoms.items()
                if atom_name not in BACKBONE_ATOMS
                and element != "H"
                and not atom_name.startswith("H")
            ]
            if side_chain:
                cma = _centroid(side_chain)
                method = "geometric centroid of {} side-chain heavy atoms".format(len(side_chain))
            elif residue_name == "GLY":
                cma = (ca[0] + 0.01, ca[1], ca[2])
                method = "glycine CA + 0.01 A x-offset"
            else:
                raise ValueError(
                    "{}:{} {} has no side-chain heavy atoms".format(chain or "_", key[1], original_name)
                )
            for atom_name, xyz, element in (("CA", ca, "C"), ("O", oxygen, "O"), ("CMA", cma, "C")):
                atom_lines.append(_atom_line(serial, atom_name, residue_name, output_chain, output_number, xyz, element))
                serial += 1
            input_selector = "{}:{}{}".format(chain or "_", key[1], key[2])
            output_selector = "{}:{}".format(output_chain, output_number)
            mapping_rows.append(
                {
                    "input_selector": input_selector,
                    "output_selector": output_selector,
                    "input_residue": original_name,
                    "output_residue": residue_name,
                    "cma_method": method,
                    "status": "converted",
                }
            )
            stride.append("C" if input_selector in closure_set else "R")
        if not stride:
            raise ValueError(
                "chain {} has no complete residues with both CA and carbonyl O".format(
                    chain or "_"
                )
            )
        stride_lines.append("STRIDE ~" + "".join(stride))

    actual_chains = "".join(output_chains[chain] for chain in per_chain)
    destination = destination.expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        "\n".join(["CBLC ~" + actual_chains] + stride_lines + atom_lines + ["END"]) + "\n",
        encoding="utf-8",
    )
    mapping_path = destination.with_name(destination.stem + ".mapping.tsv")
    with mapping_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "input_selector", "output_selector", "input_residue", "output_residue",
                "cma_method", "status",
            ),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(mapping_rows)
    return PreparedStructure(destination, mapping_path, "CBLC ~" + actual_chains, len(atom_lines), emitted_residues)


def generate_chain_regions(three_point_pdb: Path) -> str:
    """Return one explicit segment-region per chain as a conservative starter."""

    residues: Dict[str, List[int]] = OrderedDict()
    seen = set()
    for raw in three_point_pdb.expanduser().read_text(encoding="utf-8", errors="replace").splitlines():
        if raw[:6].strip().upper() != "ATOM":
            continue
        line = raw.ljust(80)
        key = (line[21].strip() or "A", line[22:26].strip(), line[26].strip())
        if key in seen:
            continue
        seen.add(key)
        residues.setdefault(key[0], []).append(int(key[1]))
    if not residues:
        raise ValueError("three-point PDB contains no residues")
    blocks = []
    for chain, numbers in residues.items():
        first, last = numbers[0], numbers[-1]
        center = numbers[(len(numbers) - 1) // 2]
        blocks.append(
            "~region[\n"
            "    \\element_top_type{{segment}}\n"
            "    \\dependency_type{{independent}}\n"
            "    \\nseg{{1}}\n"
            "    \\ncenter{{1}}\n"
            "    \\segments_firstres{{{}:{}}}\n"
            "    \\segments_lastres{{{}:{}}}\n"
            "    \\segments_baseres{{{}:{}}}\n"
            "    \\centers{{{}:{}}}\n"
            "    \\prop_trans_sig{{1.2e-05}}\n"
            "    \\prop_rot_sig{{2.375e-05}}\n"
            "    \\prop_trans_sig_freeres{{0.0}}\n"
            "    \\prop_rot_sig_freeres{{0.0}}\n"
            "]\n".format(chain, first, chain, last, chain, center, chain, center)
        )
    return "\n".join(blocks)


def write_chain_regions(path: Path, three_point_pdb: Path) -> Path:
    path = path.expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(generate_chain_regions(three_point_pdb), encoding="utf-8")
    return path
