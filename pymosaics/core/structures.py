"""PDB inspection, RCSB download, chain/model selection, and MOSAICS naming."""

import csv
import math
import os
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
from urllib.request import Request, urlopen


DNA_INTERNAL = {"A": "ADD", "C": "CYD", "G": "GUD", "T": "THD"}
DNA_TERMINAL_5 = {"A": "AD5", "C": "CD5", "G": "GD5", "T": "TD5"}
DNA_TERMINAL_3 = {"A": "AD3", "C": "CD3", "G": "GD3", "T": "TD3"}
RNA_INTERNAL = {"A": "ADE", "C": "CYT", "G": "GUA", "U": "URA"}
RNA_TERMINAL_5 = {"A": "AR5", "C": "CR5", "G": "GR5", "U": "UR5"}
RNA_TERMINAL_3 = {"A": "AR3", "C": "CR3", "G": "GR3", "U": "UR3"}
MOSAICS_DNA = {
    **{value: key for key, value in DNA_INTERNAL.items()},
    **{value: key for key, value in DNA_TERMINAL_5.items()},
    **{value: key for key, value in DNA_TERMINAL_3.items()},
}
MOSAICS_RNA = {
    **{value: key for key, value in RNA_INTERNAL.items()},
    **{value: key for key, value in RNA_TERMINAL_5.items()},
    **{value: key for key, value in RNA_TERMINAL_3.items()},
}
THYMINE_ATOMS = {"C7": "C5M", "H71": "H51", "H72": "H52", "H73": "H53"}
RCSB_ID = re.compile(r"^[A-Za-z0-9]{4}$")


@dataclass(frozen=True)
class ResidueIdentifier:
    chain: str
    number: str
    insertion_code: str
    name: str

    @property
    def selector(self) -> str:
        suffix = self.insertion_code if self.insertion_code else ""
        return "{}:{}{}".format(self.chain or "_", self.number, suffix)


@dataclass(frozen=True)
class PDBMetadata:
    models: Tuple[str, ...]
    chains_by_model: Dict[str, Tuple[str, ...]]
    residues_by_model: Dict[str, Tuple[ResidueIdentifier, ...]]


@dataclass(frozen=True)
class PreparedStructure:
    pdb_path: Path
    mapping_path: Optional[Path]
    header_line: str
    atom_count: int
    residue_count: int


@dataclass(frozen=True)
class DisulfideCandidate:
    first: ResidueIdentifier
    second: ResidueIdentifier
    distance_angstrom: float
    evidence: str

    @property
    def key(self) -> Tuple[str, str]:
        return (self.first.selector, self.second.selector)


def _model_number(line: str, fallback: int) -> str:
    value = line[10:14].strip() if len(line) >= 14 else ""
    return value or str(fallback)


def inspect_pdb_text(text: str) -> PDBMetadata:
    models: List[str] = []
    chains: Dict[str, List[str]] = {}
    residues: Dict[str, List[ResidueIdentifier]] = {}
    residue_keys: Dict[str, set] = {}
    current_model = "1"
    model_counter = 0

    for raw in text.splitlines():
        record = raw[:6].strip().upper()
        if record == "MODEL":
            model_counter += 1
            current_model = _model_number(raw, model_counter)
            continue
        if record not in ("ATOM", "HETATM"):
            continue
        if current_model not in models:
            models.append(current_model)
            chains[current_model] = []
            residues[current_model] = []
            residue_keys[current_model] = set()
        line = raw.ljust(80)
        chain = line[21].strip()
        if chain not in chains[current_model]:
            chains[current_model].append(chain)
        key = (chain, line[22:26].strip(), line[26].strip())
        if key not in residue_keys[current_model]:
            residue_keys[current_model].add(key)
            residues[current_model].append(
                ResidueIdentifier(chain, key[1], key[2], line[17:20].strip())
            )

    if not models:
        raise ValueError("PDB contains no ATOM or HETATM records")
    return PDBMetadata(
        tuple(models),
        {key: tuple(value) for key, value in chains.items()},
        {key: tuple(value) for key, value in residues.items()},
    )


def inspect_pdb(path: Path) -> PDBMetadata:
    return inspect_pdb_text(path.expanduser().read_text(encoding="utf-8", errors="replace"))


def detect_disulfides_text(text: str, model: str = "1", cutoff: float = 2.5) -> Tuple[DisulfideCandidate, ...]:
    """Detect cysteine pairs using the same 2.5 Å SG cutoff as MOSAICS."""

    sg_atoms = []
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
        if line[12:16].strip().upper() != "SG" or line[17:20].strip().upper() not in ("CYS", "CYX"):
            continue
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        residue = ResidueIdentifier(
            line[21].strip(), line[22:26].strip(), line[26].strip(), line[17:20].strip()
        )
        sg_atoms.append((residue, xyz))

    candidates = []
    for index, (first, first_xyz) in enumerate(sg_atoms):
        for second, second_xyz in sg_atoms[index + 1 :]:
            distance = math.sqrt(sum((a - b) ** 2 for a, b in zip(first_xyz, second_xyz)))
            if distance <= cutoff:
                candidates.append(DisulfideCandidate(first, second, distance, "SG–SG geometry"))
    candidates.sort(key=lambda item: item.distance_angstrom)
    return tuple(candidates)


def detect_disulfides(path: Path, model: str = "1", cutoff: float = 2.5) -> Tuple[DisulfideCandidate, ...]:
    return detect_disulfides_text(
        path.expanduser().read_text(encoding="utf-8", errors="replace"), model=model, cutoff=cutoff
    )


def unambiguous_disulfide_keys(
    candidates: Sequence[DisulfideCandidate],
) -> Tuple[Tuple[str, str], ...]:
    """Return only isolated one-to-one sulfur pairs safe to select automatically."""

    occurrence_count: Dict[str, int] = {}
    for candidate in candidates:
        for selector in candidate.key:
            occurrence_count[selector] = occurrence_count.get(selector, 0) + 1
    return tuple(
        candidate.key
        for candidate in candidates
        if all(occurrence_count[selector] == 1 for selector in candidate.key)
    )


def fetch_rcsb_pdb(pdb_id: str, destination: Path, timeout: int = 20) -> Path:
    identifier = pdb_id.strip().upper()
    if not RCSB_ID.match(identifier):
        raise ValueError("RCSB identifier must contain exactly four letters or digits")
    request = Request(
        "https://files.rcsb.org/download/{}.pdb".format(identifier),
        headers={"User-Agent": "PymoSAICS/0.2 (+https://github.com/omagebright/PymoSAICS)"},
    )
    with urlopen(request, timeout=timeout) as response:
        content = response.read(20 * 1024 * 1024 + 1)
    if len(content) > 20 * 1024 * 1024:
        raise ValueError("RCSB response exceeds the 20 MB safety limit")
    text = content.decode("utf-8", errors="replace")
    inspect_pdb_text(text)
    destination = destination.expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=str(destination.parent), delete=False) as handle:
            handle.write(text)
            temporary = handle.name
        os.replace(temporary, destination)
    finally:
        if temporary and Path(temporary).exists():
            Path(temporary).unlink()
    return destination


def _selected_atom_lines(text: str, model: str, selected_chains: Sequence[str]) -> List[str]:
    output: List[str] = []
    current_model = "1"
    model_counter = 0
    chains = set(selected_chains)
    for raw in text.splitlines():
        record = raw[:6].strip().upper()
        if record == "MODEL":
            model_counter += 1
            current_model = _model_number(raw, model_counter)
            continue
        if record not in ("ATOM", "HETATM") or current_model != model:
            continue
        line = raw.ljust(80)
        if line[21].strip() not in chains:
            continue
        altloc = line[16].strip()
        if altloc not in ("", "A"):
            continue
        output.append(line.rstrip())
    if not output:
        raise ValueError("the selected model and chains contain no atoms")
    return output


def _classify_residue(name: str) -> Tuple[str, str]:
    value = name.strip().upper()
    if value in MOSAICS_DNA:
        return "dna", MOSAICS_DNA[value]
    if value in MOSAICS_RNA:
        return "rna", MOSAICS_RNA[value]
    stem = value[:-1] if value.endswith(("3", "5")) else value
    if stem in ("DA", "DC", "DG", "DT"):
        return "dna", stem[1]
    if stem in ("A", "C", "G", "U"):
        return "rna", stem
    if stem in ("RA", "RC", "RG", "RU"):
        return "rna", stem[1]
    raise ValueError("unsupported nucleic-acid residue name: {}".format(name))


def _target_residue(family: str, base: str, profile: str, position: str) -> str:
    if profile == "standard":
        return (DNA_INTERNAL if family == "dna" else RNA_INTERNAL)[base]
    if profile != "terminal":
        raise ValueError("unsupported topology profile: {}".format(profile))
    if position == "single":
        raise ValueError(
            "a one-residue chain cannot use the terminal profile because MOSAICS has separate 5-prime and 3-prime templates"
        )
    if position == "first":
        return (DNA_TERMINAL_5 if family == "dna" else RNA_TERMINAL_5)[base]
    if position == "last":
        return (DNA_TERMINAL_3 if family == "dna" else RNA_TERMINAL_3)[base]
    return (DNA_INTERNAL if family == "dna" else RNA_INTERNAL)[base]


def _atom_name(family: str, base: str, name: str) -> str:
    value = name.strip()
    if value == "OP1":
        return "O1P"
    if value == "OP2":
        return "O2P"
    if family == "rna" and value == "H2'":
        return "H2''"
    if family == "rna" and value == "HO2'":
        return "H2'"
    if family == "dna" and base == "T":
        return THYMINE_ATOMS.get(value, value)
    return value


def _format_atom(name: str) -> str:
    return ("{:<4}" if len(name) >= 4 or (name and name[0].isdigit()) else " {:<3}").format(name)[:4]


def _rename_nucleic_atoms(lines: List[str], profile: str) -> Tuple[List[str], List[Dict[str, str]]]:
    residue_order: Dict[str, List[Tuple[str, str, str]]] = {}
    residue_names: Dict[Tuple[str, str, str], str] = {}
    for raw in lines:
        line = raw.ljust(80)
        key = (line[21].strip(), line[22:26].strip(), line[26].strip())
        residue_names[key] = line[17:20].strip()
        residue_order.setdefault(key[0], [])
        if key not in residue_order[key[0]]:
            residue_order[key[0]].append(key)

    residue_targets: Dict[Tuple[str, str, str], Tuple[str, str, str, str, bool]] = {}
    for chain, keys in residue_order.items():
        for index, key in enumerate(keys):
            position = "single" if len(keys) == 1 else "first" if index == 0 else "last" if index == len(keys) - 1 else "internal"
            family, base = _classify_residue(residue_names[key])
            input_name = residue_names[key].strip().upper()
            already_mosaics = input_name in MOSAICS_DNA or input_name in MOSAICS_RNA
            residue_targets[key] = (
                family,
                base,
                position,
                _target_residue(family, base, profile, position),
                already_mosaics,
            )

    output: List[str] = []
    mapping: List[Dict[str, str]] = []
    for raw in lines:
        line = raw.ljust(80)
        key = (line[21].strip(), line[22:26].strip(), line[26].strip())
        family, base, position, residue, already_mosaics = residue_targets[key]
        before_atom = line[12:16].strip()
        after_atom = before_atom if already_mosaics else _atom_name(family, base, before_atom)
        renamed = line[:12] + _format_atom(after_atom) + line[16:17] + "{:>3}".format(residue) + line[20:]
        output.append(renamed.rstrip())
        mapping.append(
            {
                "chain": key[0],
                "residue_number": key[1],
                "position": position,
                "family": family,
                "input_residue": residue_names[key],
                "output_residue": residue,
                "input_atom": before_atom,
                "output_atom": after_atom,
            }
        )
    return output, mapping


def _write_mapping(path: Path, rows: List[Dict[str, str]]) -> None:
    columns = (
        "chain",
        "residue_number",
        "position",
        "family",
        "input_residue",
        "output_residue",
        "input_atom",
        "output_atom",
    )
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def prepare_structure(
    source: Path,
    destination: Path,
    model: str,
    chains: Sequence[str],
    chemistry: str,
    topology_profile: str,
    header_mode: str,
    disulfide_pairs: Sequence[Tuple[str, str]] = (),
) -> PreparedStructure:
    """Write the selected structure and a visible CBLC/SCBLC header."""

    source = source.expanduser().resolve()
    text = source.read_text(encoding="utf-8", errors="replace")
    atom_lines = _selected_atom_lines(text, model, chains)
    mapping_rows: List[Dict[str, str]] = []
    if chemistry == "nucleic_acid":
        atom_lines, mapping_rows = _rename_nucleic_atoms(atom_lines, topology_profile)
    elif chemistry == "protein" and disulfide_pairs:
        selected_cyx = {selector for pair in disulfide_pairs for selector in pair}
        renamed = []
        for raw in atom_lines:
            line = raw.ljust(80)
            selector = ResidueIdentifier(
                line[21].strip(), line[22:26].strip(), line[26].strip(), line[17:20].strip()
            ).selector
            if selector in selected_cyx:
                if line[12:16].strip().upper() in ("HG", "HG1"):
                    continue
                line = line[:17] + "CYX" + line[20:]
            renamed.append(line.rstrip())
        atom_lines = renamed

    actual_chains: List[str] = []
    residue_keys = set()
    for raw in atom_lines:
        line = raw.ljust(80)
        chain = line[21].strip() or "A"
        if chain not in actual_chains:
            actual_chains.append(chain)
        residue_keys.add((chain, line[22:26].strip(), line[26].strip()))

    if header_mode == "regular":
        header = "CBLC >" + "".join(actual_chains)
    elif header_mode == "successive":
        header = "CBLC ~" + "".join(actual_chains)
    elif header_mode == "none":
        header = ""
    else:
        raise ValueError("unsupported MOSAICS header mode: {}".format(header_mode))

    destination = destination.expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    lines = ([header] if header else []) + atom_lines + ["END"]
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")
    mapping_path = None
    if mapping_rows:
        mapping_path = destination.with_name(destination.stem + ".mapping.tsv")
        _write_mapping(mapping_path, mapping_rows)
    return PreparedStructure(destination, mapping_path, header, len(atom_lines), len(residue_keys))
