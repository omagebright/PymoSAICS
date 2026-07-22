import math
import shutil
import subprocess
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

from pymosaics.core.builder import default_settings, generate_mcmc_input, stage_force_field
from pymosaics.core.catalog import analysis_preset, force_field_profile, runtime_profile
from pymosaics.core.landscape import iter_pdb_frames
from pymosaics.core.trajectory import glycosidic_chi


FIXTURE = Path(__file__).parent / "fixtures" / "single_gud_cblc.pdb"


def _frame_atoms(frame):
    atoms = {}
    for raw in frame:
        line = raw.ljust(80)
        try:
            atoms[line[12:16].strip()] = (
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            )
        except ValueError:
            continue
    return atoms


def _unwrapped_span(values):
    unwrapped = [values[0]]
    for value in values[1:]:
        delta = (value - unwrapped[-1] + 180.0) % 360.0 - 180.0
        unwrapped.append(unwrapped[-1] + delta)
    return max(unwrapped) - min(unwrapped)


class NucleicRuntimeTests(unittest.TestCase):
    def test_experimental_runtime_loads_every_nucleic_profile(self):
        runtime = runtime_profile("mosaics-experimental-2026-07-21")
        executable = runtime.executable()
        if not runtime.available() or executable is None:
            self.skipTest("bundled experimental runtime is not executable on this host")

        preset = analysis_preset("single-point")
        profile_ids = (
            "bs0-standard",
            "bsc1-ol3-standard",
            "ol15-ol3-standard",
            "ol15-ol3-terminal",
            "ol21-ol3-standard",
            "ol21-ol3-terminal",
            "ol24-ol3-standard",
            "ol24-ol3-terminal",
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for identifier in profile_ids:
                with self.subTest(profile=identifier):
                    project = root / identifier
                    project.mkdir()
                    shutil.copy2(FIXTURE, project / "structure.pdb")
                    profile = force_field_profile(identifier)
                    stage_force_field(project, profile)
                    input_path = project / "mcmc.input"
                    input_path.write_text(
                        generate_mcmc_input(
                            preset,
                            profile,
                            "structure.pdb",
                            "simulation",
                            default_settings(preset),
                        ),
                        encoding="utf-8",
                    )
                    result = subprocess.run(
                        (str(executable), str(input_path)),
                        cwd=str(project),
                        capture_output=True,
                        text=True,
                        timeout=15,
                    )
                    self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_stable_runtime_exposes_and_moves_glycosidic_chi(self):
        runtime = runtime_profile("mosaics-3.9.1")
        executable = runtime.executable()
        if not runtime.available() or executable is None:
            self.skipTest("bundled stable runtime is not executable on this host")

        preset = analysis_preset("cblc-regression")
        profile = force_field_profile("bs0-standard")
        settings = replace(
            default_settings(preset),
            replica_number=19,
            energy_gap=math.pow(1000.0 / 300.0, 1.0 / 19.0),
        )
        with tempfile.TemporaryDirectory() as temporary:
            project = Path(temporary)
            shutil.copy2(FIXTURE, project / "structure.pdb")
            stage_force_field(project, profile)
            input_path = project / "mcmc.input"
            input_path.write_text(
                generate_mcmc_input(preset, profile, "structure.pdb", "simulation", settings),
                encoding="utf-8",
            )
            result = subprocess.run(
                (str(executable), str(input_path)),
                cwd=str(project),
                capture_output=True,
                text=True,
                timeout=30,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

            torsion_rows = [
                line.split()
                for line in (project / "simulation.torsions.dat").read_text().splitlines()
                if line.startswith("MODEL ")
            ]
            self.assertTrue(torsion_rows)
            self.assertEqual(len(torsion_rows[0]) - 2, 9)

            with (project / "simulation.trajectory.pdb").open() as handle:
                frames = tuple(iter_pdb_frames(handle))
            self.assertGreater(len(frames), 2)
            chi_values = [glycosidic_chi(_frame_atoms(frame)) for frame in frames]
            self.assertGreater(_unwrapped_span(chi_values), 5.0)


if __name__ == "__main__":
    unittest.main()
