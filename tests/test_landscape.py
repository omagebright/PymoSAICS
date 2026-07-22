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


def atom_line(serial, atom, residue, chain, number, xyz, element="C"):
    x, y, z = xyz
    return (
        "ATOM  {serial:5d} {atom:>4s} {residue:>3s} {chain:1s}{number:4d}    "
        "{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00          {element:>2s}"
    ).format(**locals())


@unittest.skipIf(np is None, "NumPy is optional")
class LandscapeTests(unittest.TestCase):
    def test_kabsch_removes_rigid_rotation_and_translation(self):
        first = np.asarray(((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)), dtype=float)
        rotation = np.asarray(((0, -1, 0), (1, 0, 0), (0, 0, 1)), dtype=float)
        second = first @ rotation + np.asarray((5, -2, 3))
        self.assertAlmostEqual(kabsch_rmsd(first, second), 0.0, places=12)

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
