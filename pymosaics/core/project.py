"""MOSAICS input preparation and project validation."""

import hashlib
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

from .models import Diagnostic, PreparedRun, RuntimeConfig
from .runtime import build_command, has_errors, validate_runtime


PLACEHOLDERS = ("${PYMOSAICS_FORCEFIELD_DIR}", "${PROJECT_DIR}")
ENTRY_PATTERN = re.compile(r"\\([A-Za-z][A-Za-z0-9_]*)\{([^{}\r\n]*)\}")
INPUT_FILE_KEYS = {
    "mol_parm_file",
    "bond_database_file",
    "bend_database_file",
    "tors_database_file",
    "onfo_database_file",
    "inter_database_file",
    "region_database_file",
    "pos_init_file",
}
OUTPUT_PATTERNS = ("simulation.pdb", "simulation_result.pdb", "*.pos_out.pdb")


class PreparationError(RuntimeError):
    """Raised when an invalid project is prepared for execution."""


def _portable(path: Path) -> str:
    return str(path.expanduser().resolve()).replace("\\", "/")


def resolve_placeholders(text: str, project_directory: Path, config: RuntimeConfig) -> str:
    replacements: Dict[str, str] = {
        "${PYMOSAICS_FORCEFIELD_DIR}": _portable(config.forcefield_directory),
        "${PROJECT_DIR}": _portable(project_directory),
    }
    for source, destination in replacements.items():
        text = text.replace(source, destination)
    return text


def planned_parameter_input(parameter_input: Path, config: RuntimeConfig) -> Path:
    """Return the exact input path that prepare_run will pass to MOSAICS."""

    parameter_input = parameter_input.expanduser().resolve()
    source = parameter_input.read_text(encoding="utf-8")
    resolved = resolve_placeholders(source, parameter_input.parent, config)
    if resolved == source:
        return parameter_input
    digest = hashlib.sha256(resolved.encode("utf-8")).hexdigest()[:12]
    return parameter_input.parent / ".pymosaics" / "resolved" / (
        "{}-{}.input".format(parameter_input.stem, digest)
    )


def referenced_input_files(text: str, project_directory: Path) -> Tuple[Path, ...]:
    references: List[Path] = []
    for key, value in ENTRY_PATTERN.findall(text):
        if key not in INPUT_FILE_KEYS:
            continue
        value = value.strip()
        if not value:
            continue
        path = Path(value).expanduser()
        if not path.is_absolute():
            path = project_directory / path
        references.append(path)
    return tuple(references)


def validate_project(parameter_input: Path, config: RuntimeConfig) -> Tuple[Diagnostic, ...]:
    diagnostics: List[Diagnostic] = list(validate_runtime(config))
    parameter_input = parameter_input.expanduser()

    if not parameter_input.is_file():
        diagnostics.append(Diagnostic("error", "Parameter input does not exist: {}".format(parameter_input)))
        return tuple(diagnostics)

    try:
        source = parameter_input.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as exc:
        diagnostics.append(Diagnostic("error", "Cannot read parameter input: {}".format(exc)))
        return tuple(diagnostics)

    resolved = resolve_placeholders(source, parameter_input.parent, config)
    unresolved = sorted(set(re.findall(r"\$\{[A-Za-z_][A-Za-z0-9_]*\}", resolved)))
    if unresolved:
        diagnostics.append(
            Diagnostic("error", "Unresolved input placeholder(s): {}".format(", ".join(unresolved)))
        )

    references = referenced_input_files(resolved, parameter_input.parent)
    for path in references:
        if not path.is_file():
            diagnostics.append(Diagnostic("error", "Referenced input file does not exist: {}".format(path)))
    if not references:
        diagnostics.append(
            Diagnostic("warning", "No recognized input-file references were found; review the parameter file manually")
        )

    if not any(item.severity == "error" for item in diagnostics):
        diagnostics.append(Diagnostic("info", "Project input is ready to run"))
    return tuple(diagnostics)


def prepare_run(parameter_input: Path, config: RuntimeConfig) -> PreparedRun:
    parameter_input = parameter_input.expanduser().resolve()
    diagnostics = validate_project(parameter_input, config)
    if has_errors(diagnostics):
        messages = "; ".join(item.message for item in diagnostics if item.severity == "error")
        raise PreparationError(messages)

    source = parameter_input.read_text(encoding="utf-8")
    resolved = resolve_placeholders(source, parameter_input.parent, config)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    state_directory = parameter_input.parent / ".pymosaics"
    log_directory = state_directory / "logs"
    log_directory.mkdir(parents=True, exist_ok=True)

    if resolved != source:
        resolved_input = planned_parameter_input(parameter_input, config)
        resolved_input.parent.mkdir(parents=True, exist_ok=True)
        resolved_input.write_text(resolved, encoding="utf-8")
    else:
        resolved_input = parameter_input

    log_file = log_directory / ("run-{}.log".format(stamp))
    return PreparedRun(
        command=build_command(config.executable, resolved_input),
        working_directory=parameter_input.parent,
        source_input=parameter_input,
        resolved_input=resolved_input,
        log_file=log_file,
    )


def discover_outputs(project_directory: Path) -> Tuple[Path, ...]:
    candidates = []
    seen = set()
    for pattern in OUTPUT_PATTERNS:
        for path in project_directory.glob(pattern):
            resolved = path.resolve()
            if resolved.is_file() and resolved not in seen:
                seen.add(resolved)
                candidates.append(resolved)
    candidates.sort(key=lambda path: path.stat().st_mtime, reverse=True)
    return tuple(candidates)


def format_diagnostics(diagnostics: Iterable[Diagnostic]) -> str:
    icons = {"error": "ERROR", "warning": "WARNING", "info": "OK"}
    return "\n".join("[{}] {}".format(icons.get(item.severity, item.severity.upper()), item.message) for item in diagnostics)
