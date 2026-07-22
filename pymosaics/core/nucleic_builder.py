"""Validated plans for constructing nucleic acids through a PyMOL adapter."""

from dataclasses import dataclass
from typing import Optional, Tuple


BUILD_KINDS = {
    "single-dna": (("DNA", "A"),),
    "single-rna": (("RNA", "A"),),
    "dna-duplex": (("DNA", "A"), ("DNA", "B")),
    "rna-duplex": (("RNA", "A"), ("RNA", "B")),
    "dna-rna-hybrid": (("DNA", "A"), ("RNA", "B")),
}
VALID_BASES = set("ACGTU")
SUGAR_HYDROGENS = {
    "C5'": ("H5'", "H5''"),
    "O5'": ("HO5'",),
    "C4'": ("H4'",),
    "C3'": ("H3'",),
    "C1'": ("H1'",),
    "O3'": ("HO3'",),
}
BASE_HYDROGENS = {
    "A": {"C2": ("H2",), "N6": ("H61", "H62"), "C8": ("H8",)},
    "C": {"N4": ("H41", "H42"), "C5": ("H5",), "C6": ("H6",)},
    "G": {"N1": ("H1",), "N2": ("H21", "H22"), "C8": ("H8",)},
    "T": {
        "N3": ("H3",),
        "C5M": ("H51", "H52", "H53"),
        "C7": ("H51", "H52", "H53"),
        "C6": ("H6",),
    },
    "U": {"N3": ("H3",), "C5": ("H5",), "C6": ("H6",)},
}


@dataclass(frozen=True)
class NucleicStrandPlan:
    sequence: str
    polymer: str
    form: str
    chain: str
    expected_pucker: str


@dataclass(frozen=True)
class NucleicBuildPlan:
    kind: str
    strands: Tuple[NucleicStrandPlan, ...]
    template_form: str
    warnings: Tuple[str, ...]


def _normalize_sequence(sequence: str) -> str:
    normalized = "".join(sequence.upper().split())
    if not normalized:
        raise ValueError("a nucleotide sequence is required")
    unsupported = sorted(set(normalized) - VALID_BASES)
    if unsupported:
        raise ValueError("unsupported nucleotide code(s): {}".format(", ".join(unsupported)))
    return normalized


def _for_polymer(sequence: str, polymer: str) -> str:
    return sequence.replace("U", "T") if polymer == "DNA" else sequence.replace("T", "U")


def reverse_complement(sequence: str, target_polymer: str) -> str:
    target_polymer = target_polymer.upper()
    if target_polymer not in ("DNA", "RNA"):
        raise ValueError("target polymer must be DNA or RNA")
    normalized = _normalize_sequence(sequence).replace("U", "T")
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    result = "".join(complement[base] for base in reversed(normalized))
    return _for_polymer(result, target_polymer)


def hydrogen_names_for_parent(
    residue_name: str, parent_atom: str, terminal_hydrogens: bool = True
) -> Tuple[str, ...]:
    """Return force-field-facing hydrogen names for one PyMOL heavy atom."""

    residue = residue_name.strip().upper()
    if residue.startswith("D") and len(residue) >= 2:
        polymer, base = "DNA", residue[1]
    elif residue.startswith("R") and len(residue) >= 2:
        polymer, base = "RNA", residue[1]
    else:
        polymer, base = ("DNA", "T") if residue == "T" else ("RNA", residue[:1])
    atom = parent_atom.strip().upper().replace("*", "'")
    if atom in ("O1P", "O2P", "OP1", "OP2"):
        return ()
    if atom == "C2'":
        # These are source-PDB names. prepare_structure maps RNA H2'/HO2'
        # onto the MOSAICS H2''/H2' convention without ambiguity.
        return ("H2'", "H2''") if polymer == "DNA" else ("H2'",)
    if atom == "O2'":
        return ("HO2'",) if polymer == "RNA" else ()
    if atom in ("O5'", "O3'") and not terminal_hydrogens:
        return ()
    if atom in SUGAR_HYDROGENS:
        return SUGAR_HYDROGENS[atom]
    return BASE_HYDROGENS.get(base, {}).get(atom, ())


def plan_nucleic_acid_build(
    kind: str,
    sequence1: str,
    strand1_form: str,
    strand2_form: Optional[str] = None,
    sequence2: Optional[str] = None,
) -> NucleicBuildPlan:
    if kind not in BUILD_KINDS:
        raise ValueError("unsupported nucleic-acid build kind: {}".format(kind))
    forms = (strand1_form.upper(), (strand2_form or strand1_form).upper())
    if any(form not in ("A", "B") for form in forms):
        raise ValueError("strand form must be A or B")

    definitions = BUILD_KINDS[kind]
    first_polymer = definitions[0][0]
    first_sequence = _for_polymer(_normalize_sequence(sequence1), first_polymer)
    strands = [
        NucleicStrandPlan(
            first_sequence,
            first_polymer,
            forms[0],
            "A",
            "A-like / C3'-endo" if forms[0] == "A" else "B-like / C2'-endo",
        )
    ]
    warnings = []
    if first_polymer == "RNA" and forms[0] == "B":
        warnings.append("B-form RNA is a noncanonical structural hypothesis; A-form is recommended.")

    if len(definitions) == 2:
        second_polymer = definitions[1][0]
        expected = reverse_complement(first_sequence, second_polymer)
        if sequence2 is not None and sequence2.strip():
            supplied = _for_polymer(_normalize_sequence(sequence2), second_polymer)
            if supplied != expected:
                raise ValueError(
                    "strand 2 must be the antiparallel Watson-Crick reverse complement ({})".format(expected)
                )
        second_sequence = expected
        strands.append(
            NucleicStrandPlan(
                second_sequence,
                second_polymer,
                forms[1],
                "B",
                "A-like / C3'-endo" if forms[1] == "A" else "B-like / C2'-endo",
            )
        )
        if second_polymer == "RNA" and forms[1] == "B":
            warnings.append("B-form RNA is a noncanonical structural hypothesis; A-form is recommended.")
        if forms[0] != forms[1]:
            warnings.append(
                "Mixed A/B strand forms are fitted to one duplex template and should be relaxed before production use."
            )

    template_form = "A" if any(strand.polymer == "RNA" for strand in strands) else forms[0]
    return NucleicBuildPlan(kind, tuple(strands), template_form, tuple(dict.fromkeys(warnings)))
