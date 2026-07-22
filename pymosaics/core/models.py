"""Data models shared by the PymoSAICS core and GUI."""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple


@dataclass(frozen=True)
class RuntimeConfig:
    """Locations supplied by the user for an external MOSAICS installation."""

    executable: Path
    forcefield_directory: Path
    default_workspace: Optional[Path] = None


@dataclass(frozen=True)
class Diagnostic:
    """A validation result suitable for both tests and user interfaces."""

    severity: str
    message: str


@dataclass(frozen=True)
class PreparedRun:
    """All paths and arguments needed to start one MOSAICS process."""

    command: Tuple[str, ...]
    working_directory: Path
    source_input: Path
    resolved_input: Path
    log_file: Path
