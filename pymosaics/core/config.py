"""Cross-platform user configuration for PymoSAICS."""

import json
import os
import sys
import tempfile
from pathlib import Path
from typing import Mapping, Optional

from .models import RuntimeConfig


class ConfigError(RuntimeError):
    """Raised when the saved configuration cannot be read safely."""


def config_directory(
    platform_name: Optional[str] = None,
    environment: Optional[Mapping[str, str]] = None,
    home: Optional[Path] = None,
) -> Path:
    """Return the native per-user configuration directory."""

    platform_name = platform_name or sys.platform
    environment = environment if environment is not None else os.environ
    home = Path(home) if home is not None else Path.home()

    if platform_name.startswith("win"):
        appdata = environment.get("APPDATA")
        base = Path(appdata) if appdata else home / "AppData" / "Roaming"
        return base / "PymoSAICS"
    if platform_name == "darwin":
        return home / "Library" / "Application Support" / "PymoSAICS"

    xdg_config_home = environment.get("XDG_CONFIG_HOME")
    base = Path(xdg_config_home) if xdg_config_home else home / ".config"
    return base / "pymosaics"


class ConfigStore:
    """Read and atomically write the plugin's JSON configuration."""

    SCHEMA_VERSION = 2

    def __init__(self, path: Optional[Path] = None) -> None:
        override = os.environ.get("PYMOSAICS_CONFIG_FILE")
        self.path = Path(path or override or (config_directory() / "config.json"))

    def load(self) -> Optional[RuntimeConfig]:
        if not self.path.exists():
            return None
        try:
            payload = json.loads(self.path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError) as exc:
            raise ConfigError("Cannot read PymoSAICS configuration: {}".format(exc)) from exc

        schema_version = payload.get("schema_version")
        if schema_version not in (1, self.SCHEMA_VERSION):
            raise ConfigError("Unsupported PymoSAICS configuration version")

        try:
            workspace = payload.get("default_workspace") or None
            return RuntimeConfig(
                executable=Path(payload["executable"]).expanduser(),
                forcefield_directory=Path(payload["forcefield_directory"]).expanduser(),
                default_workspace=Path(workspace).expanduser() if workspace else None,
                runtime_id=payload.get("runtime_id", "custom"),
                force_field_id=payload.get("force_field_id", "ol24-ol3-standard"),
            )
        except (KeyError, TypeError) as exc:
            raise ConfigError("PymoSAICS configuration is incomplete") from exc

    def save(self, config: RuntimeConfig) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "schema_version": self.SCHEMA_VERSION,
            "executable": str(config.executable.expanduser().resolve()),
            "forcefield_directory": str(config.forcefield_directory.expanduser().resolve()),
            "default_workspace": (
                str(config.default_workspace.expanduser().resolve())
                if config.default_workspace is not None
                else ""
            ),
            "runtime_id": config.runtime_id,
            "force_field_id": config.force_field_id,
        }

        temporary_name = None
        try:
            with tempfile.NamedTemporaryFile(
                mode="w",
                encoding="utf-8",
                dir=str(self.path.parent),
                prefix="config-",
                suffix=".tmp",
                delete=False,
            ) as handle:
                json.dump(payload, handle, indent=2, sort_keys=True)
                handle.write("\n")
                temporary_name = handle.name
            os.replace(temporary_name, self.path)
            if os.name != "nt":
                self.path.chmod(0o600)
        except OSError as exc:
            if temporary_name:
                try:
                    Path(temporary_name).unlink()
                except OSError:
                    pass
            raise ConfigError("Cannot save PymoSAICS configuration: {}".format(exc)) from exc
