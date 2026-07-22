import tempfile
import unittest
from pathlib import Path

from pymosaics.core.regions import RegionSettings, generate_region_file
from pymosaics.core.structures import (
    detect_disulfides_text,
    fetch_rcsb_pdb,
    inspect_pdb_text,
    prepare_structure,
    unambiguous_disulfide_keys,
)


def atom_line(serial, atom, residue, chain, number, x, y, z, element):
    return (
        "ATOM  {serial:5d} {atom:>4s} {residue:>3s} {chain:1s}{number:4d}    "
        "{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00          {element:>2s}"
    ).format(**locals())


class StructureTests(unittest.TestCase):
    def test_rcsb_identifier_is_validated_before_network_access(self):
        with self.assertRaisesRegex(ValueError, "four"):
            fetch_rcsb_pdb("1BNA;rm", Path("unused.pdb"))

    def test_inspects_models_chains_and_residues(self):
        text = "\n".join(
            (
                "MODEL        1",
                atom_line(1, "P", "DA", "A", 1, 0, 0, 0, "P"),
                atom_line(2, "P", "U", "B", 1, 1, 0, 0, "P"),
                "ENDMDL",
                "MODEL        2",
                atom_line(3, "P", "DA", "A", 1, 0, 1, 0, "P"),
                "ENDMDL",
            )
        )
        metadata = inspect_pdb_text(text)
        self.assertEqual(metadata.models, ("1", "2"))
        self.assertEqual(metadata.chains_by_model["1"], ("A", "B"))

    def test_nucleic_acid_standard_and_terminal_naming(self):
        text = "\n".join(
            (
                atom_line(1, "P", "DA", "A", 1, 0, 0, 0, "P"),
                atom_line(2, "OP1", "DA", "A", 1, 1, 0, 0, "O"),
                atom_line(3, "P", "DC", "A", 2, 2, 0, 0, "P"),
                "END",
            )
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source.pdb"
            source.write_text(text)
            standard = prepare_structure(source, root / "standard.pdb", "1", ("A",), "nucleic_acid", "standard", "regular")
            standard_text = standard.pdb_path.read_text()
            self.assertTrue(standard_text.startswith("CBLC >A\n"))
            self.assertIn("ADD", standard_text)
            self.assertIn("CYD", standard_text)
            self.assertIn("O1P", standard_text)
            terminal = prepare_structure(source, root / "terminal.pdb", "1", ("A",), "nucleic_acid", "terminal", "successive")
            terminal_text = terminal.pdb_path.read_text()
            self.assertTrue(terminal_text.startswith("CBLC ~A\n"))
            self.assertIn("AD5", terminal_text)
            self.assertIn("CD3", terminal_text)

    def test_existing_mosaics_atom_names_are_not_converted_twice(self):
        text = "\n".join(
            (
                atom_line(1, "H2'", "URA", "A", 1, 0, 0, 0, "H"),
                atom_line(2, "H2''", "URA", "A", 1, 1, 0, 0, "H"),
                atom_line(3, "P", "ADE", "A", 2, 2, 0, 0, "P"),
                "END",
            )
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source.pdb"
            source.write_text(text)
            prepared = prepare_structure(source, root / "prepared.pdb", "1", ("A",), "nucleic_acid", "standard", "none")
            atom_names = [line[12:16].strip() for line in prepared.pdb_path.read_text().splitlines() if line.startswith("ATOM")]
            self.assertEqual(atom_names[:2], ["H2'", "H2''"])

    def test_disulfide_detection_and_amber_processing(self):
        text = "\n".join(
            (
                atom_line(1, "SG", "CYS", "A", 4, 0, 0, 0, "S"),
                atom_line(2, "HG", "CYS", "A", 4, 0, 0, 1, "H"),
                atom_line(3, "SG", "CYS", "A", 9, 2.03, 0, 0, "S"),
                atom_line(4, "HG", "CYS", "A", 9, 2.03, 0, 1, "H"),
                "END",
            )
        )
        candidates = detect_disulfides_text(text)
        self.assertEqual(len(candidates), 1)
        self.assertAlmostEqual(candidates[0].distance_angstrom, 2.03)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source.pdb"
            source.write_text(text)
            prepared = prepare_structure(
                source,
                root / "prepared.pdb",
                "1",
                ("A",),
                "protein",
                "protein",
                "none",
                (candidates[0].key,),
            )
            result = prepared.pdb_path.read_text()
            self.assertEqual(result.count("CYX"), 2)
            self.assertNotIn(" HG ", result)

    def test_only_isolated_disulfide_pairs_are_recommended(self):
        text = "\n".join(
            (
                atom_line(1, "SG", "CYS", "A", 1, 0, 0, 0, "S"),
                atom_line(2, "SG", "CYS", "A", 2, 2.0, 0, 0, "S"),
                atom_line(3, "SG", "CYS", "A", 3, 0, 2.0, 0, "S"),
                atom_line(4, "SG", "CYS", "A", 4, 10, 0, 0, "S"),
                atom_line(5, "SG", "CYS", "A", 5, 12.0, 0, 0, "S"),
                "END",
            )
        )
        candidates = detect_disulfides_text(text)
        self.assertEqual(unambiguous_disulfide_keys(candidates), (("A:4", "A:5"),))


class RegionTests(unittest.TestCase):
    def test_graphical_settings_generate_parser_keywords(self):
        text = generate_region_file(
            RegionSettings(
                residues=("A:1", "A:2", "B:3"),
                centers=("A:1",),
                residue_pairs=(("A:2", "B:3"),),
            )
        )
        self.assertIn("\\nres{3}", text)
        self.assertIn("\\centers{A:1}", text)
        self.assertIn("\\residue_pair{A:2,B:3}", text)

    def test_region_rejects_unknown_center(self):
        with self.assertRaisesRegex(ValueError, "belong"):
            generate_region_file(RegionSettings(("A:1",), ("A:2",), ()))

    def test_region_rejects_invalid_widths_and_duplicate_pairs(self):
        with self.assertRaisesRegex(ValueError, "numeric"):
            generate_region_file(
                RegionSettings(("A:1",), (), (), translation_sigma="not-a-number")
            )
        with self.assertRaisesRegex(ValueError, "only once"):
            generate_region_file(
                RegionSettings(
                    ("A:1", "A:2"),
                    (),
                    (("A:1", "A:2"), ("A:2", "A:1")),
                )
            )


if __name__ == "__main__":
    unittest.main()
