import os
import tempfile
import time
import unittest
from pathlib import Path

from pymosaics.core.models import RuntimeConfig
from pymosaics.core.project import (
    MosaicsInputDocument,
    PreparationError,
    discover_outputs,
    import_project_input,
    make_portable_input,
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
        self.assertEqual(prepared.log_file.parent, self.project.resolve() / "logs")
        self.assertNotIn(".pymosaics", prepared.log_file.parts)

    def test_prepare_archives_existing_outputs_before_a_new_run(self):
        trajectory = self.project / "simulation.trajectory.pdb"
        energy = self.project / "simulation.potential_energy.dat"
        fixed = self.project / "sim_param.out"
        trajectory.write_text("old trajectory\n", encoding="utf-8")
        energy.write_text("old energy\n", encoding="utf-8")
        fixed.write_text("old parameters\n", encoding="utf-8")
        parameter_input = self._input(
            "\\atom_pos_file{simulation.trajectory.pdb}\n"
            "\\epot_file{simulation.potential_energy.dat}\n"
        )

        prepared = prepare_run(parameter_input, self.config)

        self.assertFalse(trajectory.exists())
        self.assertFalse(energy.exists())
        self.assertFalse(fixed.exists())
        self.assertEqual(len(prepared.archived_outputs), 3)
        for archived in prepared.archived_outputs:
            self.assertTrue(archived.is_file())
            self.assertIn("run_history", archived.parts)

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

    def test_input_document_round_trips_unknown_and_repeated_directives(self):
        source = (
            "# keep this scientific note\n"
            "~sim_gen_def[\n"
            "  \\temperature{310.15}\n"
            "  \\total_step_mc{100000}\n"
            "]\n"
            "~sim_mol_def[\n"
            "  \\energy_term{inter}\n"
            "  \\energy_term{cryo_em}\n"
            "]\n"
        )
        document = MosaicsInputDocument(source)
        self.assertEqual(document.value("temperature"), "310.15")
        self.assertEqual(document.values("energy_term"), ("inter", "cryo_em"))
        updated = document.updated({"temperature": "300", "total_step_mc": "80000"})
        self.assertIn("# keep this scientific note", updated)
        self.assertIn("\\temperature{300}", updated)
        self.assertIn("\\total_step_mc{80000}", updated)
        self.assertEqual(updated.count("\\energy_term{"), 2)

    def test_tom_style_absolute_paths_are_remapped_by_unique_local_basename(self):
        for name in ("top_3pt_prot_na.rtf", "par_3pt_prot_na.prm", "3pt_strides.pdb", "regions.data"):
            (self.project / name).write_text(name, encoding="utf-8")
        source = self.project / "mosaics.input"
        source.write_text(
            "~sim_mol_def[\n"
            "  \\mol_parm_file{/Users/tom/mo_temp/input/top_3pt_prot_na.rtf}\n"
            "  \\bond_database_file{/Users/tom/mo_temp/input/par_3pt_prot_na.prm}\n"
            "  \\pos_init_file{/Users/tom/mo_temp/input/3pt_strides.pdb}\n"
            "  \\region_database_file{/Users/tom/mo_temp/input/regions.data}\n"
            "  \\pos_out_file{/Users/tom/mo_temp/output/sampled.pdb}\n"
            "  \\atom_pos_file{/Users/tom/mo_temp/output/sampled.pos}\n"
            "  \\param_out_file{/Users/tom/mo_temp/output/sim_param.out}\n"
            "]\n",
            encoding="utf-8",
        )
        portable = make_portable_input(source.read_text(), self.project)
        self.assertEqual(portable.unresolved, ())
        self.assertIn("\\pos_init_file{${PROJECT_DIR}/3pt_strides.pdb}", portable.content)
        self.assertIn("\\pos_out_file{${PROJECT_DIR}/output/sampled.pdb}", portable.content)
        self.assertNotIn("param_out_file", portable.content)

        managed = import_project_input(source, self.project / "mcmc.input")
        self.assertEqual(managed.input_path, (self.project / "mcmc.input").resolve())
        self.assertTrue((self.project / "output").is_dir())
        self.assertIn("/Users/tom/", source.read_text())
        self.assertNotIn("/Users/tom/", managed.input_path.read_text())

    def test_ambiguous_foreign_basename_is_reported_not_guessed(self):
        (self.project / "a").mkdir()
        (self.project / "b").mkdir()
        for parent in (self.project / "a", self.project / "b"):
            (parent / "structure.pdb").write_text("END\n", encoding="utf-8")
        source = "\\pos_init_file{/another/computer/structure.pdb}\n"
        result = make_portable_input(source, self.project)
        self.assertEqual(result.content, source)
        self.assertEqual(len(result.unresolved), 1)
        self.assertIn("ambiguous", result.unresolved[0].reason)

    def test_nested_historical_deck_resolves_paths_from_its_own_directory(self):
        level = self.project / "level1"
        database = self.project / "top_database"
        level.mkdir()
        database.mkdir()
        topology = database / "top_3pt_prot_na.rtf"
        topology.write_text("topology", encoding="utf-8")
        source = level / "refine.input"
        source.write_text(
            "\\mol_parm_file{../top_database/top_3pt_prot_na.rtf}\n",
            encoding="utf-8",
        )
        imported = import_project_input(
            source,
            self.project / "mcmc.input",
            project_directory=self.project,
        )
        self.assertIn(
            "\\mol_parm_file{${PROJECT_DIR}/top_database/top_3pt_prot_na.rtf}",
            imported.input_path.read_text(),
        )


if __name__ == "__main__":
    unittest.main()
