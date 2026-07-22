"""Transparent AMBER-compatible all-atom preparation through PDB2PQR."""

import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence, Tuple

from .structures import _selected_atom_lines


@dataclass(frozen=True)
class ProteinPreparationResult:
    pdb_path: Path
    pqr_path: Path
    selected_input_path: Path
    log_path: Path
    command: Tuple[str, ...]
    ph: float


def discover_pdb2pqr(explicit: Optional[Path] = None) -> Optional[Path]:
    """Find a local PDB2PQR command without modifying the host environment."""

    candidates = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    environment = os.environ.get("PYMOSAICS_PDB2PQR")
    if environment:
        candidates.append(Path(environment).expanduser())
    executable_directory = Path(sys.executable).resolve().parent
    for name in ("pdb2pqr", "pdb2pqr30", "pdb2pqr.exe", "pdb2pqr30.exe"):
        candidates.append(executable_directory / name)
        discovered = shutil.which(name)
        if discovered:
            candidates.append(Path(discovered))

    visited = set()
    for candidate in candidates:
        try:
            resolved = candidate.resolve()
        except OSError:
            continue
        marker = os.path.normcase(str(resolved))
        if marker in visited:
            continue
        visited.add(marker)
        if resolved.is_file():
            return resolved
    return None


def _write_selected_protein(source: Path, destination: Path, model: str, chains: Sequence[str]) -> None:
    text = source.expanduser().resolve().read_text(encoding="utf-8", errors="replace")
    selected = [
        line
        for line in _selected_atom_lines(text, model, chains)
        if line[:6].strip().upper() == "ATOM"
    ]
    if not selected:
        raise ValueError("the selected model and chains contain no protein ATOM records")

    output = []
    previous_chain = None
    for raw in selected:
        line = raw.ljust(80)
        chain = line[21].strip()
        if previous_chain is not None and chain != previous_chain:
            output.append("TER")
        output.append(line.rstrip())
        previous_chain = chain
    output.extend(("TER", "END"))
    destination.write_text("\n".join(output) + "\n", encoding="utf-8")


def _pqr_atom_to_pdb(raw: str) -> str:
    line = raw.rstrip("\r\n").ljust(70)
    try:
        serial = int(line[6:11])
        atom = line[12:16]
        residue = line[17:20].strip()
        chain = line[21].strip()
        residue_number = int(line[22:26])
        insertion_code = line[26].strip()
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
    except (ValueError, IndexError):
        fields = raw.split()
        if len(fields) < 10:
            raise ValueError("PDB2PQR produced an unreadable atom record: {}".format(raw))
        try:
            serial = int(fields[1])
            atom = fields[2]
            residue = fields[3]
            has_chain = not fields[4].lstrip("+-").isdigit()
            chain = fields[4] if has_chain else ""
            residue_number = int(fields[5] if has_chain else fields[4])
            coordinate_start = 6 if has_chain else 5
            x, y, z = map(float, fields[coordinate_start : coordinate_start + 3])
            insertion_code = ""
        except (ValueError, IndexError) as exc:
            raise ValueError("PDB2PQR produced an unreadable atom record: {}".format(raw)) from exc

    if not 0 <= serial <= 99999:
        raise ValueError("PDB2PQR atom serial exceeds the PDB format limit")
    if not -999 <= residue_number <= 9999:
        raise ValueError("PDB2PQR residue number exceeds the PDB format limit")
    if len(chain) > 1:
        raise ValueError("PDB2PQR produced a multi-character chain identifier")
    atom_name = atom.strip()
    if not atom_name:
        raise ValueError("PDB2PQR produced an atom without a name")
    formatted_atom = (
        "{:<4}" if len(atom_name) >= 4 or atom_name[0].isdigit() else " {:<3}"
    ).format(atom_name)[:4]
    element = next((character.upper() for character in atom_name if character.isalpha()), "")
    return (
        "ATOM  {serial:5d} {atom} {residue:>3s} {chain:1s}{residue_number:4d}{icode:1s}   "
        "{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}  "
    ).format(
        serial=serial,
        atom=formatted_atom,
        residue=residue,
        chain=chain,
        residue_number=residue_number,
        icode=insertion_code,
        x=x,
        y=y,
        z=z,
        element=element,
    )


def pqr_to_pdb(source: Path, destination: Path) -> Path:
    """Convert PDB2PQR coordinates to a strict PDB without hiding its PQR output."""

    output = []
    for raw in source.read_text(encoding="utf-8", errors="replace").splitlines():
        record = raw[:6].strip().upper()
        if record in ("ATOM", "HETATM"):
            output.append(_pqr_atom_to_pdb(raw))
        elif record == "TER" and output and output[-1] != "TER":
            output.append("TER")
    if not any(line.startswith("ATOM") for line in output):
        raise ValueError("PDB2PQR output contains no protein atoms")
    if output[-1] != "TER":
        output.append("TER")
    output.append("END")
    destination.write_text("\n".join(output) + "\n", encoding="utf-8")
    return destination


def _pdb2pqr_options(help_text: str, ph: float) -> Tuple[str, ...]:
    if "--keep-chain" in help_text:
        chain_option = "--keep-chain"
    elif "--chain" in help_text:
        chain_option = "--chain"
    else:
        raise ValueError("this PDB2PQR version cannot preserve chain identifiers")

    if "--titration-state-method" in help_text:
        ph_method = "--titration-state-method=propka"
    elif "--ph-calc-method" in help_text:
        ph_method = "--ph-calc-method=propka"
    else:
        raise ValueError("this PDB2PQR version does not expose PROPKA titration-state assignment")
    if "--with-ph" not in help_text:
        raise ValueError("this PDB2PQR version does not expose an explicit pH setting")
    return (
        "--ff=amber",
        "--ffout=amber",
        chain_option,
        "--drop-water",
        ph_method,
        "--with-ph={:.2f}".format(ph),
    )


def prepare_protein_with_pdb2pqr(
    source: Path,
    output_directory: Path,
    executable: Path,
    model: str,
    chains: Sequence[str],
    ph: float = 7.0,
    timeout: int = 300,
) -> ProteinPreparationResult:
    """Repair/protonate a selected protein with visible AMBER/PROPKA provenance."""

    if not 0.0 <= ph <= 14.0:
        raise ValueError("protein preparation pH must be between 0 and 14")
    command_path = discover_pdb2pqr(executable)
    if command_path is None or command_path != Path(executable).expanduser().resolve():
        raise ValueError("the selected PDB2PQR executable does not exist")

    output_directory = output_directory.expanduser().resolve()
    output_directory.mkdir(parents=True, exist_ok=True)
    selected_input = output_directory / "selected-protein.pdb"
    pqr_path = output_directory / "pdb2pqr-amber.pqr"
    pdb_path = output_directory / "pdb2pqr-amber.pdb"
    log_path = output_directory / "pdb2pqr.log"
    _write_selected_protein(source, selected_input, model, chains)

    try:
        help_process = subprocess.run(
            (str(command_path), "--help"),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=30,
            check=False,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        raise ValueError("PDB2PQR could not be started: {}".format(exc)) from exc
    help_text = help_process.stdout or ""
    options = _pdb2pqr_options(help_text, ph)
    command = options + (str(selected_input), str(pqr_path))
    full_command = (str(command_path),) + command
    log_header = "PymoSAICS protein preparation\npH: {:.2f}\nArguments:\n{}\n\n".format(
        ph, "\n".join(full_command)
    )
    try:
        process = subprocess.run(
            full_command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        log_path.write_text(log_header + "PDB2PQR timed out.\n", encoding="utf-8")
        raise ValueError("PDB2PQR timed out after {} seconds; inspect {}".format(timeout, log_path)) from exc
    except OSError as exc:
        log_path.write_text(log_header + str(exc) + "\n", encoding="utf-8")
        raise ValueError("PDB2PQR could not be started; inspect {}".format(log_path)) from exc

    log_path.write_text(log_header + (process.stdout or ""), encoding="utf-8")
    if process.returncode != 0 or not pqr_path.is_file():
        raise ValueError(
            "PDB2PQR could not prepare this protein (exit code {}); inspect {}".format(
                process.returncode, log_path
            )
        )
    pqr_to_pdb(pqr_path, pdb_path)
    return ProteinPreparationResult(
        pdb_path=pdb_path,
        pqr_path=pqr_path,
        selected_input_path=selected_input,
        log_path=log_path,
        command=full_command,
        ph=ph,
    )
