import tempfile
import unittest
import zipfile
from pathlib import Path

from scripts.build_release import VERSION, build


class ReleaseTests(unittest.TestCase):
    def test_plugin_manager_archive_contains_assets_and_executable_modes(self):
        with tempfile.TemporaryDirectory() as temporary:
            archive = build(Path(temporary))
            self.assertEqual(archive.name, "PymoSAICS-{}.zip".format(VERSION))
            with zipfile.ZipFile(archive) as bundle:
                names = set(bundle.namelist())
                stable = "pymosaics/assets/runtimes/darwin-arm64/mosaics-3.9.1"
                experimental = "pymosaics/assets/runtimes/darwin-arm64/mosaics-experimental-2026-07-21"
                three_point = "pymosaics/assets/forcefields/kb_3pt/par_3pt_prot_na.prm"
                self.assertIn("pymosaics/__init__.py", names)
                self.assertIn(stable, names)
                self.assertIn(experimental, names)
                self.assertIn(three_point, names)
                self.assertEqual((bundle.getinfo(stable).external_attr >> 16) & 0o777, 0o755)
                self.assertEqual((bundle.getinfo(experimental).external_attr >> 16) & 0o777, 0o755)
                self.assertFalse(any(name.endswith((".cpp", ".h", ".o")) for name in names))


if __name__ == "__main__":
    unittest.main()
