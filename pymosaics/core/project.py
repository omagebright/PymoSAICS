"""MOSAICS input preparation and project validation."""

import hashlib
import re
import os
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

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
    "cryo_em_database_file",
    "constraint_database_file",
}
OUTPUT_FILE_KEYS = {
    "pos_out_file",
    "atom_pos_file",
    "param_out_file",
    "epot_file",
    "einter_file",
    "tors_pos_file",
    "hessian_file",
    "eighess_file",
}
UNSUPPORTED_PORTABLE_KEYS = {
    # MOSAICS 3.9.1 and the bundled experimental runtime always write this
    # record as ./sim_param.out; neither parser accepts a configurable path.
    "param_out_file": "MOSAICS writes sim_param.out in the working directory",
}
OUTPUT_PATTERNS = ("simulation.pdb", "simulation_result.pdb", "*.pos_out.pdb")


class PreparationError(RuntimeError):
    """Raised when an invalid project is prepared for execution."""


@dataclass(frozen=True)
class PathRewrite:
    key: str
    original: str
    replacement: str


@dataclass(frozen=True)
class PathIssue:
    key: str
    value: str
    reason: str


@dataclass(frozen=True)
class PortableInput:
    content: str
    rewrites: Tuple[PathRewrite, ...]
    unresolved: Tuple[PathIssue, ...]


@dataclass(frozen=True)
class ImportedInput:
    input_path: Path
    source_path: Path
    rewrites: Tuple[PathRewrite, ...]


class MosaicsInputDocument:
    """Lossless view of a MOSAICS deck with surgical directive updates.

    MOSAICS accepts repeated directives such as ``energy_term`` and historical
    decks often carry options unknown to this PymoSAICS release.  This class
    deliberately avoids reformatting or regenerating the file: callers can read
    any directive and replace selected scalar values while every other byte is
    retained.
    """

    def __init__(self, text: str):
        self.text = text

    def values(self, key: str) -> Tuple[str, ...]:
        return tuple(
            match.group(2).strip()
            for match in ENTRY_PATTERN.finditer(self.text)
            if match.group(1) == key
        )

    def value(self, key: str, default: Optional[str] = None) -> Optional[str]:
        values = self.values(key)
        return values[0] if values else default

    def updated(self, changes: Mapping[str, object]) -> str:
        text = self.text
        for key, value in changes.items():
            if not re.match(r"^[A-Za-z][A-Za-z0-9_]*$", key):
                raise ValueError("invalid MOSAICS directive name: {}".format(key))
            replacement = str(value)
            if any(character in replacement for character in "{}\r\n"):
                raise ValueError("invalid value for MOSAICS directive {}".format(key))
            pattern = re.compile(r"(\\{}\{{)[^{{}}\r\n]*(\}})".format(re.escape(key)))
            text, count = pattern.subn(
                lambda match: match.group(1) + replacement + match.group(2),
                text,
                count=1,
            )
            if count == 0:
                raise KeyError("MOSAICS directive is not present: {}".format(key))
        return text


def _project_file_index(project_directory: Path) -> Dict[str, Tuple[Path, ...]]:
    grouped: Dict[str, List[Path]] = {}
    for path in project_directory.rglob("*"):
        if not path.is_file() or ".pymosaics" in path.parts or "output" in path.parts:
            continue
        grouped.setdefault(path.name, []).append(path.resolve())
    return {
        name: tuple(sorted(paths, key=lambda item: item.as_posix()))
        for name, paths in grouped.items()
    }


def _project_placeholder(path: Path, project_directory: Path) -> str:
    relative = path.resolve().relative_to(project_directory.resolve())
    return "${PROJECT_DIR}/" + relative.as_posix()


def make_portable_input(
    text: str, project_directory: Path, source_directory: Optional[Path] = None
) -> PortableInput:
    """Remap foreign absolute paths without guessing between local candidates.

    Input files are matched by unique basename within the project. Output files
    are moved beneath a visible ``output`` directory. The returned text is safe
    to move because all rewritten paths use ``${PROJECT_DIR}``.
    """

    project_directory = project_directory.expanduser().resolve()
    source_directory = (
        source_directory.expanduser().resolve()
        if source_directory is not None
        else project_directory
    )
    index = _project_file_index(project_directory)
    rewrites: List[PathRewrite] = []
    unresolved: List[PathIssue] = []

    def replace(match: re.Match) -> str:
        key, raw_value = match.group(1), match.group(2)
        value = raw_value.strip()
        if not value or key not in INPUT_FILE_KEYS | OUTPUT_FILE_KEYS:
            return match.group(0)
        if key in UNSUPPORTED_PORTABLE_KEYS:
            explanation = UNSUPPORTED_PORTABLE_KEYS[key]
            rewrites.append(PathRewrite(key, value, "removed: " + explanation))
            return ""
        if "${PROJECT_DIR}" in value or "${PYMOSAICS_FORCEFIELD_DIR}" in value:
            return match.group(0)

        path = Path(value).expanduser()
        replacement = None
        if key in OUTPUT_FILE_KEYS:
            name = path.name
            if not name:
                unresolved.append(PathIssue(key, value, "output path has no filename"))
                return match.group(0)
            replacement = "${PROJECT_DIR}/output/" + name
        else:
            candidate = path if path.is_absolute() else source_directory / path
            try:
                if candidate.is_file() and candidate.resolve().is_relative_to(project_directory):
                    replacement = _project_placeholder(candidate, project_directory)
            except (OSError, ValueError):
                replacement = None
            if replacement is None:
                candidates = index.get(path.name, ())
                if len(candidates) == 1:
                    replacement = _project_placeholder(candidates[0], project_directory)
                elif len(candidates) > 1:
                    unresolved.append(
                        PathIssue(key, value, "ambiguous basename matches {} local files".format(len(candidates)))
                    )
                else:
                    unresolved.append(PathIssue(key, value, "no local file has this basename"))
        if replacement is None or replacement == value:
            return match.group(0)
        rewrites.append(PathRewrite(key, value, replacement))
        return "\\{}{{{}}}".format(key, replacement)

    content = ENTRY_PATTERN.sub(replace, text)
    return PortableInput(content, tuple(rewrites), tuple(unresolved))


def import_project_input(
    source_input: Path,
    destination_input: Optional[Path] = None,
    project_directory: Optional[Path] = None,
) -> ImportedInput:
    """Create a managed portable copy of an existing MOSAICS input deck."""

    source_input = source_input.expanduser().resolve()
    if not source_input.is_file():
        raise PreparationError("MOSAICS input does not exist: {}".format(source_input))
    destination_input = (
        destination_input.expanduser().resolve()
        if destination_input is not None
        else source_input.parent / "mcmc.input"
    )
    project_directory = (
        project_directory.expanduser().resolve()
        if project_directory is not None
        else destination_input.parent
    )
    result = make_portable_input(
        source_input.read_text(encoding="utf-8"),
        project_directory,
        source_directory=source_input.parent,
    )
    if result.unresolved:
        details = "; ".join(
            "{}={}: {}".format(issue.key, issue.value, issue.reason)
            for issue in result.unresolved
        )
        raise PreparationError("Cannot import project portably: " + details)

    destination_input.parent.mkdir(parents=True, exist_ok=True)
    (destination_input.parent / "output").mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=str(destination_input.parent),
            prefix=".pymosaics-import-",
            suffix=".tmp",
            delete=False,
        ) as handle:
            handle.write(result.content)
            temporary = Path(handle.name)
        os.replace(str(temporary), str(destination_input))
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()
    return ImportedInput(destination_input, source_input, result.rewrites)


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
    # A run log is a primary scientific record, not internal application state.
    # Keep it in a visible project directory so it is easy to find, archive,
    # compare, and share outside PymoSAICS.
    log_directory = parameter_input.parent / "logs"
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
