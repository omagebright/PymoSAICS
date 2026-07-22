import tempfile
import unittest
from pathlib import Path

from pymosaics.core.analysis import (
    discover_energy_series,
    discover_pdb_outputs,
    discover_project_files,
    parse_acceptance_log,
    read_text_file,
)


class AnalysisTests(unittest.TestCase):
    def test_energy_and_acceptance_parsing(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "simulation.potential_energy.dat").write_text("header\n-10.5 0\n-8.0 1\n")
            log = root / "run.log"
            log.write_text("Chain 0 25 100 0.25\nnoise\n")
            series = discover_energy_series(root)
            self.assertEqual(series[0].values, (-10.5, -8.0))
            self.assertEqual(parse_acceptance_log(log)[0].ratio, 0.25)

    def test_project_files_include_logs_and_identify_pdbs(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "result.pdb").write_text("END\n")
            logs = root / ".pymosaics" / "logs"
            logs.mkdir(parents=True)
            (logs / "run.log").write_text("ok\n")
            hidden = root / ".pymosaics" / "resolved"
            hidden.mkdir()
            (hidden / "copy.input").write_text("internal\n")
            files = discover_project_files(root)
            self.assertEqual({item.path.name for item in files}, {"result.pdb", "run.log"})
            pdb = next(item for item in files if item.path.suffix == ".pdb")
            self.assertTrue(pdb.text_readable)
            self.assertTrue(pdb.loadable_in_pymol)

    def test_text_reader_rejects_binary(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "binary.dat"
            path.write_bytes(b"a\x00b")
            with self.assertRaisesRegex(ValueError, "binary"):
                read_text_file(path)

    def test_generated_input_structure_is_not_reported_as_an_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "structure.pdb").write_text("END\n")
            result = root / "simulation.pos_out.pdb"
            result.write_text("END\n")
            outputs = discover_pdb_outputs(root)
        self.assertEqual(outputs, (result.resolve(),))


if __name__ == "__main__":
    unittest.main()
