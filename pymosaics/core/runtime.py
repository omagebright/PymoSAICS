"""Validation and safe process arguments for an external MOSAICS runtime."""

import os
import sys
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

from .models import Diagnostic, RuntimeConfig


def validate_runtime(
    config: RuntimeConfig, platform_name: Optional[str] = None
) -> Tuple[Diagnostic, ...]:
    diagnostics: List[Diagnostic] = []
    platform_name = platform_name or sys.platform
    executable = config.executable.expanduser()
    forcefields = config.forcefield_directory.expanduser()

    if not executable.is_file():
        diagnostics.append(Diagnostic("error", "MOSAICS executable does not exist: {}".format(executable)))
    elif platform_name.startswith("win"):
        if executable.suffix.lower() not in (".exe", ".com"):
            diagnostics.append(Diagnostic("error", "Select a Windows .exe or .com MOSAICS executable"))
    elif not os.access(str(executable), os.X_OK):
        diagnostics.append(Diagnostic("error", "MOSAICS file is not executable: {}".format(executable)))

    if not forcefields.is_dir():
        diagnostics.append(Diagnostic("error", "Force-field directory does not exist: {}".format(forcefields)))
    else:
        try:
            from .catalog import force_field_profile, runtime_supports_force_field

            profile = force_field_profile(config.force_field_id)
        except (ImportError, StopIteration):
            diagnostics.append(
                Diagnostic(
                    "error",
                    "Unknown force-field profile: {}".format(config.force_field_id),
                )
            )
        else:
            if not runtime_supports_force_field(config.runtime_id, config.force_field_id):
                diagnostics.append(
                    Diagnostic(
                        "error",
                        "{} is not compatible with force-field profile {}".format(
                            config.runtime_id, profile.label
                        ),
                    )
                )
            missing = [path for path in profile.all_paths() if not path.is_file()]
            if missing:
                diagnostics.append(
                    Diagnostic(
                        "error",
                        "Selected force-field profile is incomplete: {}".format(
                            ", ".join(path.name for path in missing)
                        ),
                    )
                )
            else:
                diagnostics.append(
                    Diagnostic("info", "Force-field profile is complete: {}".format(profile.label))
                )

    if config.default_workspace is not None and not config.default_workspace.expanduser().is_dir():
        diagnostics.append(
            Diagnostic("warning", "Default workspace does not exist: {}".format(config.default_workspace))
        )

    if not any(item.severity == "error" for item in diagnostics):
        diagnostics.append(Diagnostic("info", "External MOSAICS runtime is configured"))
    return tuple(diagnostics)


def build_command(executable: Path, parameter_input: Path) -> Tuple[str, ...]:
    """Return an argument vector; callers must never pass it through a shell."""

    return (str(executable.expanduser().resolve()), str(parameter_input.expanduser().resolve()))


def has_errors(diagnostics: Sequence[Diagnostic]) -> bool:
    return any(item.severity == "error" for item in diagnostics)
