import os
import tempfile
import unittest
from pathlib import Path

from pymosaics.core.catalog import FORCEFIELD_ROOT
from pymosaics.core.models import RuntimeConfig
from pymosaics.core.runtime import build_command, has_errors, validate_runtime


class RuntimeTests(unittest.TestCase):
    def _runtime_tree(self, root: Path, executable_name="mosaics.x"):
        executable = root / executable_name
        executable.write_text("runtime", encoding="utf-8")
        executable.chmod(0o755)
        forcefields = root / "force fields"
        (forcefields / "top_database").mkdir(parents=True)
        (forcefields / "pot_database").mkdir()
        return RuntimeConfig(executable, forcefields)

    def test_valid_posix_runtime(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = self._runtime_tree(Path(temporary))
            self.assertFalse(has_errors(validate_runtime(config, "darwin")))

    def test_windows_requires_executable_extension(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = self._runtime_tree(Path(temporary))
            self.assertTrue(has_errors(validate_runtime(config, "win32")))

    def test_windows_exe_is_accepted(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = self._runtime_tree(Path(temporary), "mosaics.exe")
            self.assertFalse(has_errors(validate_runtime(config, "win32")))

    def test_command_is_an_argument_tuple_not_a_shell_string(self):
        with tempfile.TemporaryDirectory(prefix="project with spaces; ") as temporary:
            root = Path(temporary)
            executable = root / "mosaics executable"
            parameter_input = root / "parameters input.input"
            executable.write_text("", encoding="utf-8")
            parameter_input.write_text("", encoding="utf-8")
            command = build_command(executable, parameter_input)
            self.assertIsInstance(command, tuple)
            self.assertEqual(len(command), 2)
            self.assertIn("spaces", command[0])
            self.assertNotIn(" > ", command)

    def test_stable_runtime_rejects_newer_force_field_profile(self):
        with tempfile.TemporaryDirectory() as temporary:
            executable = Path(temporary) / "mosaics"
            executable.write_text("runtime", encoding="utf-8")
            executable.chmod(0o755)
            config = RuntimeConfig(
                executable,
                FORCEFIELD_ROOT,
                runtime_id="mosaics-3.9.1",
                force_field_id="ol24-ol3-standard",
            )
            diagnostics = validate_runtime(config, platform_name="linux")
            self.assertTrue(any("not compatible" in item.message for item in diagnostics))


if __name__ == "__main__":
    unittest.main()
