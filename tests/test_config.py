import json
import tempfile
import unittest
from pathlib import Path

from pymosaics.core.config import ConfigError, ConfigStore, config_directory
from pymosaics.core.models import RuntimeConfig


class ConfigDirectoryTests(unittest.TestCase):
    def test_windows_uses_appdata(self):
        result = config_directory("win32", {"APPDATA": "C:/Users/Test/AppData/Roaming"}, Path("C:/Users/Test"))
        self.assertEqual(str(result).replace("\\", "/"), "C:/Users/Test/AppData/Roaming/PymoSAICS")

    def test_macos_uses_application_support(self):
        result = config_directory("darwin", {}, Path("/Users/test"))
        self.assertEqual(result, Path("/Users/test/Library/Application Support/PymoSAICS"))

    def test_linux_uses_xdg_when_set(self):
        result = config_directory("linux", {"XDG_CONFIG_HOME": "/tmp/config"}, Path("/home/test"))
        self.assertEqual(result, Path("/tmp/config/pymosaics"))


class ConfigStoreTests(unittest.TestCase):
    def test_round_trip(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "config.json"
            store = ConfigStore(path)
            expected = RuntimeConfig(
                executable=Path(temporary) / "MOSAICS executable",
                forcefield_directory=Path(temporary) / "force fields",
                default_workspace=Path(temporary) / "workspace",
            )
            store.save(expected)
            self.assertEqual(
                store.load(),
                RuntimeConfig(
                    expected.executable.resolve(),
                    expected.forcefield_directory.resolve(),
                    expected.default_workspace.resolve(),
                ),
            )
            self.assertEqual(json.loads(path.read_text())["schema_version"], 2)

    def test_schema_one_is_migrated_with_safe_defaults(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "config.json"
            path.write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "executable": "/tmp/mosaics",
                        "forcefield_directory": "/tmp/forcefields",
                        "default_workspace": "",
                    }
                ),
                encoding="utf-8",
            )
            loaded = ConfigStore(path).load()
            self.assertEqual(loaded.runtime_id, "custom")
            self.assertEqual(loaded.force_field_id, "ol24-ol3-standard")

    def test_invalid_json_is_reported(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "config.json"
            path.write_text("not json", encoding="utf-8")
            with self.assertRaises(ConfigError):
                ConfigStore(path).load()


if __name__ == "__main__":
    unittest.main()
