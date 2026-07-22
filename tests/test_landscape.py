import math
import tempfile
import unittest
from pathlib import Path

try:
    import numpy as np
except ImportError:
    np = None

from pymosaics.core.landscape import (
    build_landscape,
    kabsch_rmsd,
    read_coordinate_frames,
    write_landscape_table,
)
from pymosaics.core.trajectory import analyze_trajectory, classify_pucker, pseudorotation_phase


def atom_line(serial, atom, residue, chain, number, xyz, element="C"):
    x, y, z = xyz
    return (
        "ATOM  {serial:5d} {atom:>4s} {residue:>3s} {chain:1s}{number:4d}    "
        "{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00          {element:>2s}"
    ).format(**locals())


@unittest.skipIf(np is None, "NumPy is optional")
class LandscapeTests(unittest.TestCase):
    def test_pseudorotation_classification_distinguishes_a_and_b_like_sugars(self):
        self.assertEqual(classify_pucker(18.0), "A-like / C3'-endo")
        self.assertEqual(classify_pucker(162.0), "B-like / C2'-endo")

    def test_pseudorotation_matches_canonical_a_rna_and_b_dna_models(self):
        # Coordinates are from PyMOL fnab canonical A-RNA and B-DNA models.
        a_rna = {
            "O4'": (6.486, 6.363, -1.919),
            "C1'": (6.874, 5.000, -1.877),
            "C2'": (7.566, 4.726, -3.208),
            "C3'": (6.776, 5.641, -4.140),
            "C4'": (6.718, 6.889, -3.258),
        }
        b_dna = {
            "O4'": (1.329, 5.973, 1.328),
            "C1'": (1.256, 5.500, -0.008),
            "C2'": (1.191, 6.754, -0.874),
            "C3'": (2.295, 7.526, -0.151),
            "C4'": (2.069, 7.228, 1.331),
        }
        a_phase = pseudorotation_phase(a_rna)
        b_phase = pseudorotation_phase(b_dna)
        self.assertGreaterEqual(a_phase, 0.0)
        self.assertLess(a_phase, 36.0)
        self.assertGreaterEqual(b_phase, 144.0)
        self.assertLess(b_phase, 180.0)
        self.assertEqual(classify_pucker(a_phase), "A-like / C3'-endo")
        self.assertEqual(classify_pucker(b_phase), "B-like / C2'-endo")

    def test_kabsch_removes_rigid_rotation_and_translation(self):
        first = np.asarray(((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)), dtype=float)
        rotation = np.asarray(((0, -1, 0), (1, 0, 0), (0, 0, 1)), dtype=float)
        second = first @ rotation + np.asarray((5, -2, 3))
        self.assertAlmostEqual(kabsch_rmsd(first, second), 0.0, places=12)

    def test_trajectory_analysis_reports_puckers_rmsd_and_terminal_mobility(self):
        ring = {
            "O4'": (6.486, 6.363, -1.919),
            "C1'": (6.874, 5.000, -1.877),
            "C2'": (7.566, 4.726, -3.208),
            "C3'": (6.776, 5.641, -4.140),
            "C4'": (6.718, 6.889, -3.258),
        }
        lines = []
        for frame in (1, 2):
            lines.append("MODEL{:9d}".format(frame))
            for serial, (name, xyz) in enumerate(ring.items(), start=1):
                shifted = np.asarray(xyz) + np.asarray((10.0, -3.0, 2.0)) * (frame - 1)
                if frame == 2 and name == "C3'":
                    shifted += np.asarray((0.35, 0.0, 0.0))
                lines.append(atom_line(serial, name, "RA", "A", 1, shifted))
            lines.append("ENDMDL")
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "rna_trajectory.pdb"
            path.write_text("\n".join(lines) + "\n")
            result = analyze_trajectory(path)
        self.assertEqual(len(result.rmsd_to_first), 2)
        self.assertGreater(result.maximum_rmsd, 0.0)
        self.assertEqual(len(result.puckers), 2)
        self.assertEqual({item.end for item in result.terminal_mobility}, {"5'", "3'"})

    def test_builds_two_dimensional_map_and_representatives(self):
        base = ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
        frames = []
        for frame, stretch in enumerate((1.0, 1.2, 1.7), start=1):
            frames.append("MODEL{:9d}".format(frame))
            for serial, point in enumerate(base, start=1):
                xyz = (point[0] * stretch, point[1], point[2])
                frames.append(atom_line(serial, "CA", "ALA", "A", serial, xyz))
            frames.append("ENDMDL")
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "trajectory.pdb"
            path.write_text("\n".join(frames) + "\n")
            result = build_landscape(path, representatives=2)
            table = write_landscape_table(Path(temporary) / "landscape.tsv", result)
            table_text = table.read_text()
        self.assertEqual(len(result.coordinates), 3)
        self.assertEqual(result.frame_numbers, (1, 2, 3))
        self.assertEqual(len(result.representative_frames), 2)
        self.assertEqual(result.atom_count, 4)
        self.assertGreater(result.rmsd_matrix[0][2], result.rmsd_matrix[0][1])
        self.assertIn("frame\tcoordinate_1\tcoordinate_2\trepresentative\tenergy", table_text)

    def test_streams_and_uniformly_samples_trajectory_frames(self):
        frames = []
        for frame in range(1, 7):
            frames.append("MODEL{:9d}".format(frame))
            for serial, point in enumerate(((0, 0, 0), (1, 0, 0), (0, 1, frame)), start=1):
                frames.append(atom_line(serial, "CA", "ALA", "A", serial, point))
            frames.append("ENDMDL")
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "trajectory.pdb"
            path.write_text("\n".join(frames) + "\n")
            coordinates, frame_numbers = read_coordinate_frames(path, maximum_frames=3)
        self.assertEqual(coordinates.shape, (3, 3, 3))
        self.assertEqual(frame_numbers, (1, 3, 6))


if __name__ == "__main__":
    unittest.main()
