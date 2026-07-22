"""Pre-run PDB-to-RTF atom validation with residue-level diagnostics."""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple


@dataclass(frozen=True)
class TopologyIssue:
    selector: str
    pdb_residue: str
    topology_residue: str
    missing_atoms: Tuple[str, ...]
    extra_atoms: Tuple[str, ...]
    duplicate_atoms: Tuple[str, ...]

    def message(self) -> str:
        details = []
        if self.missing_atoms:
            details.append("missing " + ", ".join(self.missing_atoms))
        if self.extra_atoms:
            details.append("unexpected " + ", ".join(self.extra_atoms))
        if self.duplicate_atoms:
            details.append("duplicate " + ", ".join(self.duplicate_atoms))
        return "{} {} → {}: {}".format(
            self.selector, self.pdb_residue, self.topology_residue, "; ".join(details)
        )


@dataclass(frozen=True)
class RtfChiDefinition:
    chigroups: Tuple[Tuple[str, ...], ...]
    chiprims: Tuple[Tuple[str, ...], ...]


@dataclass(frozen=True)
class RtfChiIssue:
    residue: str
    message: str


def read_rtf_chi_definitions(path: Path) -> Dict[str, RtfChiDefinition]:
    """Read CHIGROUP/CHIPRIM records without altering legacy RTF syntax."""

    records: Dict[str, Dict[str, List[Tuple[str, ...]]]] = {}
    current = None
    for raw in path.expanduser().read_text(encoding="utf-8", errors="replace").splitlines():
        fields = raw.split()
        if not fields or fields[0].startswith(("!", "*")):
            continue
        keyword = fields[0].upper()
        if keyword == "RESI" and len(fields) >= 2:
            current = fields[1].upper()
            records[current] = {"CHIGROUP": [], "CHIPRIM": []}
        elif keyword in ("PRES", "END"):
            current = None
        elif current and keyword in ("CHIGROUP", "CHIPRIM"):
            records[current][keyword].append(tuple(field.upper() for field in fields[1:]))
    return {
        residue: RtfChiDefinition(
            tuple(values["CHIGROUP"]), tuple(values["CHIPRIM"])
        )
        for residue, values in records.items()
    }


def validate_nucleic_chi_definitions(path: Path) -> Tuple[RtfChiIssue, ...]:
    """Require an explicit glycosidic chi move for every DNA/RNA residue.

    MOSAICS accepts an RTF with zero CHIGROUP and zero CHIPRIM records, but then
    constructs backbone transformations only.  Atom-level topology validation
    cannot detect that scientifically incomplete move set.
    """

    atoms = read_rtf_atom_templates(path)
    definitions = read_rtf_chi_definitions(path)
    issues = []
    for residue, atom_names in atoms.items():
        atom_set = set(atom_names)
        if not {"C1'", "C2'"}.issubset(atom_set):
            continue
        if {"N9", "C8"}.issubset(atom_set):
            primitive = ("C2'", "C1'", "N9", "C8")
        elif {"N1", "C2"}.issubset(atom_set):
            primitive = ("C2'", "C1'", "N1", "C2")
        else:
            continue
        definition = definitions.get(residue, RtfChiDefinition((), ()))
        reasons = []
        if len(definition.chigroups) != len(definition.chiprims):
            reasons.append(
                "{} CHIGROUP record(s) but {} CHIPRIM record(s)".format(
                    len(definition.chigroups), len(definition.chiprims)
                )
            )
        if primitive not in definition.chiprims:
            reasons.append("missing glycosidic CHIPRIM {}".format(" ".join(primitive)))
        if not any(primitive[3] in group for group in definition.chigroups):
            reasons.append("missing glycosidic CHIGROUP containing {}".format(primitive[3]))
        if reasons:
            issues.append(RtfChiIssue(residue, "{}: {}".format(residue, "; ".join(reasons))))
    return tuple(issues)


def read_rtf_atom_templates(path: Path) -> Dict[str, Tuple[str, ...]]:
    templates: Dict[str, List[str]] = {}
    current = None
    for raw in path.expanduser().read_text(encoding="utf-8", errors="replace").splitlines():
        fields = raw.split()
        if not fields or fields[0].startswith(("!", "*")):
            continue
        keyword = fields[0].upper()
        if keyword == "RESI" and len(fields) >= 2:
            current = fields[1].upper()
            templates[current] = []
        elif keyword in ("PRES", "END"):
            current = None
        elif keyword == "ATOM" and current and len(fields) >= 2:
            templates[current].append(fields[1].upper())
    return {name: tuple(atoms) for name, atoms in templates.items()}


def _pdb_residues(path: Path):
    residues = []
    positions = {}
    for raw in path.expanduser().read_text(encoding="utf-8", errors="replace").splitlines():
        if raw[:6].strip().upper() not in ("ATOM", "HETATM"):
            continue
        line = raw.ljust(80)
        if line[16].strip() not in ("", "A"):
            continue
        key = (line[21].strip(), line[22:26].strip(), line[26].strip())
        if key not in positions:
            positions[key] = len(residues)
            residues.append([key, line[17:20].strip().upper(), []])
        residues[positions[key]][2].append(line[12:16].strip().upper())
    return residues


def validate_pdb_against_rtf(path: Path, rtf_path: Path, chemistry: str) -> Tuple[TopologyIssue, ...]:
    templates = read_rtf_atom_templates(rtf_path)
    residues = _pdb_residues(path)
    chain_indices: Dict[str, List[int]] = {}
    for index, (key, _, _) in enumerate(residues):
        chain_indices.setdefault(key[0], []).append(index)

    issues = []
    for index, (key, residue_name, atom_names) in enumerate(residues):
        topology_name = residue_name
        if chemistry == "protein":
            chain_order = chain_indices[key[0]]
            if index == chain_order[0] and "N" + residue_name in templates:
                topology_name = "N" + residue_name
            elif index == chain_order[-1] and "C" + residue_name in templates:
                topology_name = "C" + residue_name
        expected = templates.get(topology_name)
        actual = set(atom_names)
        duplicates = tuple(sorted({name for name in atom_names if atom_names.count(name) > 1}))
        if expected is None:
            issues.append(
                TopologyIssue(
                    "{}:{}{}".format(key[0] or "_", key[1], key[2]),
                    residue_name,
                    topology_name,
                    (),
                    tuple(sorted(actual)),
                    duplicates,
                )
            )
            continue
        expected_set = set(expected)
        missing = tuple(sorted(expected_set - actual))
        extra = tuple(sorted(actual - expected_set))
        if missing or extra or duplicates or len(atom_names) != len(expected):
            issues.append(
                TopologyIssue(
                    "{}:{}{}".format(key[0] or "_", key[1], key[2]),
                    residue_name,
                    topology_name,
                    missing,
                    extra,
                    duplicates,
                )
            )
    return tuple(issues)


def format_topology_issues(issues: Tuple[TopologyIssue, ...], maximum: int = 12) -> str:
    lines = [issue.message() for issue in issues[:maximum]]
    if len(issues) > maximum:
        lines.append("… and {} more residue mismatch(es)".format(len(issues) - maximum))
    return "\n".join(lines)
