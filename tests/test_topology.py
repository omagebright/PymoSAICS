import tempfile
import unittest
from pathlib import Path

from pymosaics.core.topology import (
    read_rtf_atom_templates,
    validate_nucleic_chi_definitions,
    validate_pdb_against_rtf,
)


def atom_line(serial, atom, residue, chain, number):
    return (
        "ATOM  {serial:5d} {atom:>4s} {residue:>3s} {chain:1s}{number:4d}    "
        "{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00          {element:>2s}"
    ).format(serial=serial, atom=atom, residue=residue, chain=chain, number=number, x=serial * 1.0, y=0.0, z=0.0, element=atom[0])


class TopologyTests(unittest.TestCase):
    def test_nucleic_residues_require_glycosidic_chi_groups_and_primitives(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            missing = root / "missing.rtf"
            missing.write_text(
                "RESI GUD -1\n"
                "ATOM C2' CT 0\nATOM C1' CT 0\nATOM N9 NS 0\n"
                "ATOM C8 CK 0\nATOM C4 CB 0\nEND\n",
                encoding="utf-8",
            )
            issues = validate_nucleic_chi_definitions(missing)
            self.assertEqual(len(issues), 1)
            self.assertIn("glycosidic", issues[0].message)

            complete = root / "complete.rtf"
            complete.write_text(
                missing.read_text(encoding="utf-8").replace(
                    "END\n",
                    "CHIGROUP C8 C4\nCHIPRIM C2' C1' N9 C8\nEND\n",
                ),
                encoding="utf-8",
            )
            self.assertEqual(validate_nucleic_chi_definitions(complete), ())

    def test_reads_rtf_templates_and_reports_missing_atoms(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            rtf = root / "test.rtf"
            rtf.write_text("RESI ALA 0.0\nATOM N X 0\nATOM CA X 0\nATOM H X 0\nEND\n")
            pdb = root / "test.pdb"
            pdb.write_text(atom_line(1, "N", "ALA", "A", 1) + "\n" + atom_line(2, "CA", "ALA", "A", 1) + "\nEND\n")
            self.assertEqual(read_rtf_atom_templates(rtf)["ALA"], ("N", "CA", "H"))
            issues = validate_pdb_against_rtf(pdb, rtf, "nucleic_acid")
            self.assertEqual(issues[0].missing_atoms, ("H",))

    def test_protein_chain_ends_use_n_and_c_templates(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            rtf = root / "test.rtf"
            rtf.write_text(
                "RESI NGLY 0\nATOM N X 0\nATOM H1 X 0\n"
                "RESI CVAL 0\nATOM C X 0\nATOM OXT X 0\nEND\n"
            )
            pdb = root / "test.pdb"
            pdb.write_text(
                atom_line(1, "N", "GLY", "A", 1) + "\n"
                + atom_line(2, "H1", "GLY", "A", 1) + "\n"
                + atom_line(3, "C", "VAL", "A", 2) + "\n"
                + atom_line(4, "OXT", "VAL", "A", 2) + "\nEND\n"
            )
            self.assertEqual(validate_pdb_against_rtf(pdb, rtf, "protein"), ())


if __name__ == "__main__":
    unittest.main()
