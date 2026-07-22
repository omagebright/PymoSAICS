import os
import tempfile
import time
import unittest
from pathlib import Path

from pymosaics.core.models import RuntimeConfig
from pymosaics.core.project import (
    PreparationError,
    discover_outputs,
    planned_parameter_input,
    prepare_run,
    resolve_placeholders,
    validate_project,
)
from pymosaics.core.runtime import has_errors


class ProjectTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="PymoSAICS project ")
        self.root = Path(self.temporary.name)
        executable_name = "mosaics.exe" if os.name == "nt" else "mosaics.x"
        self.executable = self.root / "bin" / executable_name
        self.executable.parent.mkdir()
        self.executable.write_text("runtime", encoding="utf-8")
        self.executable.chmod(0o755)
        self.forcefields = self.root / "force fields"
        (self.forcefields / "top_database").mkdir(parents=True)
        (self.forcefields / "pot_database").mkdir()
        (self.forcefields / "top_database" / "test.rtf").write_text("topology", encoding="utf-8")
        self.project = self.root / "project with spaces"
        self.project.mkdir()
        (self.project / "start.pdb").write_text("END\n", encoding="utf-8")
        self.config = RuntimeConfig(self.executable, self.forcefields)

    def tearDown(self):
        self.temporary.cleanup()

    def _input(self, extra=""):
        path = self.project / "parameters.input"
        path.write_text(
            "~sim_mol_def[\n"
            "\\mol_parm_file{${PYMOSAICS_FORCEFIELD_DIR}/top_database/test.rtf}\n"
            "\\pos_init_file{${PROJECT_DIR}/start.pdb}\n"
            + extra
            + "\n]",
            encoding="utf-8",
        )
        return path

    def test_placeholders_resolve_to_portable_absolute_paths(self):
        resolved = resolve_placeholders("${PYMOSAICS_FORCEFIELD_DIR}|${PROJECT_DIR}", self.project, self.config)
        self.assertNotIn("${", resolved)
        self.assertIn(str(self.forcefields.resolve()).replace("\\", "/"), resolved)

    def test_prepare_preserves_source_and_writes_resolved_copy(self):
        parameter_input = self._input()
        original = parameter_input.read_text(encoding="utf-8")
        prepared = prepare_run(parameter_input, self.config)
        self.assertEqual(parameter_input.read_text(encoding="utf-8"), original)
        self.assertNotEqual(prepared.resolved_input, parameter_input)
        self.assertTrue(prepared.resolved_input.is_file())
        self.assertNotIn("${", prepared.resolved_input.read_text(encoding="utf-8"))
        self.assertEqual(prepared.command[0], str(self.executable.resolve()))
        self.assertEqual(prepared.working_directory, self.project.resolve())
        self.assertEqual(prepared.resolved_input, planned_parameter_input(parameter_input, self.config))

    def test_resolved_input_name_is_content_addressed(self):
        parameter_input = self._input()
        first = planned_parameter_input(parameter_input, self.config)
        second = planned_parameter_input(parameter_input, self.config)
        self.assertEqual(first, second)
        parameter_input.write_text(parameter_input.read_text() + "\n# changed\n", encoding="utf-8")
        self.assertNotEqual(first, planned_parameter_input(parameter_input, self.config))

    def test_missing_referenced_input_is_an_error(self):
        parameter_input = self._input("\\region_database_file{${PROJECT_DIR}/missing.region}")
        diagnostics = validate_project(parameter_input, self.config)
        self.assertTrue(has_errors(diagnostics))
        self.assertTrue(any("missing.region" in item.message for item in diagnostics))

    def test_unresolved_placeholder_is_an_error(self):
        parameter_input = self._input("\\region_database_file{${UNKNOWN}/region.txt}")
        diagnostics = validate_project(parameter_input, self.config)
        self.assertTrue(any("Unresolved" in item.message for item in diagnostics))

    def test_invalid_project_cannot_be_prepared(self):
        with self.assertRaises(PreparationError):
            prepare_run(self.project / "missing.input", self.config)

    def test_output_discovery_prefers_newest(self):
        older = self.project / "simulation.pdb"
        newer = self.project / "simulation_result.pdb"
        older.write_text("END\n", encoding="utf-8")
        time.sleep(0.01)
        newer.write_text("END\n", encoding="utf-8")
        self.assertEqual(discover_outputs(self.project)[0], newer.resolve())


if __name__ == "__main__":
    unittest.main()
