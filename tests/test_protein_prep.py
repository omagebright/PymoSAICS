import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from pymosaics.core.protein_prep import (
    _pdb2pqr_options,
    discover_pdb2pqr,
    pqr_to_pdb,
    prepare_protein_with_pdb2pqr,
)


PDB_ATOM = "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N"
PQR_TEXT = """ATOM      1  N   ALA A   1       0.000   0.000   0.000 -0.4157 1.8240
ATOM      2  H1  ALA A   1       0.900   0.000   0.000  0.2719 0.6000
TER
"""


class ProteinPreparationTests(unittest.TestCase):
    def test_discovers_explicit_executable(self):
        with tempfile.TemporaryDirectory() as temporary:
            executable = Path(temporary) / "pdb2pqr"
            executable.write_text("tool", encoding="utf-8")
            self.assertEqual(discover_pdb2pqr(executable), executable.resolve())

    def test_supports_old_and_current_pdb2pqr_flags(self):
        old = _pdb2pqr_options("--chain --ph-calc-method --with-ph", 6.5)
        current = _pdb2pqr_options(
            "--keep-chain --titration-state-method --with-ph", 7.4
        )
        self.assertIn("--chain", old)
        self.assertIn("--ph-calc-method=propka", old)
        self.assertIn("--with-ph=6.50", old)
        self.assertIn("--keep-chain", current)
        self.assertIn("--titration-state-method=propka", current)
        self.assertIn("--with-ph=7.40", current)

    def test_pqr_conversion_preserves_amber_names_and_pdb_columns(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "prepared.pqr"
            destination = root / "prepared.pdb"
            source.write_text(PQR_TEXT, encoding="utf-8")
            pqr_to_pdb(source, destination)
            lines = destination.read_text(encoding="utf-8").splitlines()
            self.assertEqual(lines[0][12:16].strip(), "N")
            self.assertEqual(lines[1][12:16].strip(), "H1")
            self.assertEqual(lines[0][17:20], "ALA")
            self.assertEqual(lines[-1], "END")

    def test_preparation_is_argument_based_and_keeps_a_visible_log(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source.pdb"
            source.write_text(
                PDB_ATOM + "\n"
                + "HETATM    2  O   HOH A   2       1.000   0.000   0.000  1.00  0.00           O\nEND\n",
                encoding="utf-8",
            )
            executable = root / "pdb2pqr"
            executable.write_text("tool", encoding="utf-8")

            def fake_run(arguments, **kwargs):
                self.assertIsInstance(arguments, tuple)
                self.assertNotIn("shell", kwargs)
                if arguments[-1] == "--help":
                    return SimpleNamespace(
                        stdout="--keep-chain --titration-state-method --with-ph",
                        returncode=0,
                    )
                Path(arguments[-1]).write_text(PQR_TEXT, encoding="utf-8")
                return SimpleNamespace(stdout="prepared\n", returncode=0)

            with patch("pymosaics.core.protein_prep.subprocess.run", side_effect=fake_run):
                result = prepare_protein_with_pdb2pqr(
                    source,
                    root / "output",
                    executable,
                    "1",
                    ("A",),
                    ph=7.2,
                )
            self.assertTrue(result.pdb_path.is_file())
            self.assertTrue(result.pqr_path.is_file())
            self.assertIn("pH: 7.20", result.log_path.read_text(encoding="utf-8"))
            selected = result.selected_input_path.read_text(encoding="utf-8")
            self.assertIn("ATOM", selected)
            self.assertNotIn("HETATM", selected)
            self.assertIn("--with-ph=7.20", result.command)


if __name__ == "__main__":
    unittest.main()
