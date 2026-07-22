import hashlib
import tempfile
import unittest
from unittest import mock
from pathlib import Path

from pymosaics.core.builder import default_settings, generate_mcmc_input, stage_force_field
from pymosaics.core.catalog import (
    ANALYSIS_PRESETS,
    FORCE_FIELD_PROFILES,
    RUNTIME_PROFILES,
    force_field_profile,
    host_machine,
    runtime_profile,
    runtime_supports_force_field,
)


class CatalogAndBuilderTests(unittest.TestCase):
    def test_catalog_identifiers_are_unique(self):
        for entries in (RUNTIME_PROFILES, FORCE_FIELD_PROFILES, ANALYSIS_PRESETS):
            identifiers = [item.identifier for item in entries]
            self.assertEqual(len(identifiers), len(set(identifiers)))

    def test_every_force_field_profile_is_complete(self):
        for profile in FORCE_FIELD_PROFILES:
            with self.subTest(profile=profile.identifier):
                self.assertTrue(all(path.is_file() for path in profile.all_paths()))

    def test_catalog_exposes_legacy_modern_terminal_and_protein_profiles(self):
        identifiers = {profile.identifier for profile in FORCE_FIELD_PROFILES}
        self.assertTrue(
            {
                "bs0-standard",
                "bsc1-ol3-standard",
                "ol15-ol3-standard",
                "ol15-ol3-terminal",
                "ol21-ol3-standard",
                "ol21-ol3-terminal",
                "ol24-ol3-standard",
                "ol24-ol3-terminal",
                "ff14sb-protein",
                "kb-3pt-protein",
            }.issubset(identifiers)
        )
        self.assertIn("AMBER99", force_field_profile("bs0-standard").label)
        self.assertIn("protein", force_field_profile("ff14sb-protein").label)

    def test_bundled_runtime_hashes_are_fixed(self):
        expected = {
            "mosaics-3.9.1": "a65d34474ba51c479352566928423c9560cb0776c1c37f17cae3bab59e9ab5ad",
            "mosaics-experimental-2026-07-21": "8a23058f859ff27409353b47c263280816f84c51c7d605e94897572ac1a43a15",
        }
        for identifier, digest in expected.items():
            path = runtime_profile(identifier).executable()
            self.assertTrue(path.is_file())
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), digest)

    @mock.patch("pymosaics.core.catalog.platform.machine", return_value="x86_64")
    @mock.patch("pymosaics.core.catalog.subprocess.run")
    def test_rosetta_python_detects_apple_silicon_host(self, run, machine):
        run.return_value = mock.Mock(returncode=0, stdout="1\n")
        self.assertEqual(host_machine("darwin"), "arm64")
        self.assertTrue(runtime_profile("mosaics-3.9.1").available("darwin"))

    def test_all_compatible_preset_profile_combinations_generate_explicit_input(self):
        for preset in ANALYSIS_PRESETS:
            for profile in FORCE_FIELD_PROFILES:
                if preset.chemistry not in ("any", profile.chemistry):
                    continue
                with self.subTest(preset=preset.identifier, profile=profile.identifier):
                    text = generate_mcmc_input(
                        preset,
                        profile,
                        "structure.pdb",
                        "simulation",
                        default_settings(preset),
                    )
                    self.assertIn("\\pos_init_file{structure.pdb}", text)
                    expected_files = {
                        "mol_parm_file": profile.rtf,
                        "bond_database_file": profile.bond,
                        "bend_database_file": profile.bend,
                        "tors_database_file": profile.torsion,
                        "onfo_database_file": profile.one_four,
                        "inter_database_file": profile.nonbonded,
                    }
                    for directive, path in expected_files.items():
                        self.assertIn(
                            "\\{}{{forcefield/{}}}".format(directive, Path(path).name),
                            text,
                        )
                    self.assertIn("\\simulation_typ{" + preset.simulation_type + "}", text)

    def test_landscape_ladder_reaches_documented_top_temperature(self):
        preset = next(item for item in ANALYSIS_PRESETS if item.identifier == "landscape-pilot")
        settings = default_settings(preset)
        top = settings.temperature * settings.energy_gap ** settings.replica_number
        self.assertAlmostEqual(top, preset.top_temperature, places=8)

    def test_visible_proposal_widths_are_written_verbatim(self):
        preset = ANALYSIS_PRESETS[0]
        profile = FORCE_FIELD_PROFILES[0]
        settings = default_settings(preset)
        settings = settings.__class__(
            settings.temperature,
            settings.total_steps,
            settings.statistics_frequency,
            settings.random_seed,
            settings.closure_sigma,
            settings.replica_number,
            settings.energy_gap,
            "2e-5",
            "3e-5",
            "4e-6",
        )
        text = generate_mcmc_input(preset, profile, "structure.pdb", "simulation", settings)
        self.assertIn("\\prop_tors_sig{2e-5}", text)
        self.assertIn("\\prop_trans_sig{3e-5}", text)
        self.assertIn("\\prop_rot_sig{4e-6}", text)

    def test_region_input_writes_explicit_single_region_propagation(self):
        preset = ANALYSIS_PRESETS[0]
        profile = FORCE_FIELD_PROFILES[0]
        text = generate_mcmc_input(
            preset,
            profile,
            "structure.pdb",
            "simulation",
            default_settings(preset),
            "region/region.data",
        )
        self.assertIn("\\prop_regions_type{superimpose}", text)
        self.assertIn("\\region_database_file{region/region.data}", text)

    def test_three_point_profile_writes_coarse_grained_model_and_energy_terms(self):
        preset = next(item for item in ANALYSIS_PRESETS if item.identifier == "three-point-natural-moves")
        profile = force_field_profile("kb-3pt-protein")
        text = generate_mcmc_input(
            preset,
            profile,
            "structure.pdb",
            "simulation",
            default_settings(preset),
            "region/region.data",
        )
        self.assertIn("\\cgres_model{KB_3pt}", text)
        for term in ("bond", "bend", "tors", "onfo", "inter"):
            self.assertIn("\\energy_term{" + term + "}", text)
        self.assertIn("\\prop_regions_type{onebyone}", text)

    def test_measured_runtime_force_field_compatibility(self):
        self.assertTrue(runtime_supports_force_field("mosaics-3.9.1", "bsc1-ol3-standard"))
        self.assertTrue(runtime_supports_force_field("mosaics-3.9.1", "ol15-ol3-standard"))
        self.assertFalse(runtime_supports_force_field("mosaics-3.9.1", "ol24-ol3-standard"))
        self.assertFalse(runtime_supports_force_field("mosaics-3.9.1", "ff14sb-protein"))
        self.assertTrue(runtime_supports_force_field("mosaics-experimental-2026-07-21", "ff14sb-protein"))

    def test_selected_force_field_is_staged_with_manifest(self):
        profile = force_field_profile("ol24-ol3-standard")
        with tempfile.TemporaryDirectory() as temporary:
            destination = stage_force_field(Path(temporary), profile)
            self.assertEqual(
                {path.name for path in destination.iterdir()},
                {path.name for path in profile.all_paths()} | {"PROFILE.txt"},
            )
            self.assertIn("identifier: ol24-ol3-standard", (destination / "PROFILE.txt").read_text())

    def test_staging_force_field_preserves_unrelated_user_files(self):
        profile = force_field_profile("ol24-ol3-standard")
        with tempfile.TemporaryDirectory() as temporary:
            forcefield = Path(temporary) / "forcefield"
            forcefield.mkdir()
            custom = forcefield / "my-custom-parameter.dat"
            custom.write_text("user data\n", encoding="utf-8")
            stage_force_field(Path(temporary), profile)
            self.assertEqual(custom.read_text(encoding="utf-8"), "user data\n")

    def test_input_rejects_invalid_numeric_or_embedded_path_values(self):
        preset = ANALYSIS_PRESETS[0]
        profile = FORCE_FIELD_PROFILES[0]
        settings = default_settings(preset)
        with self.assertRaisesRegex(ValueError, "filename, not a path"):
            generate_mcmc_input(preset, profile, "../structure.pdb", "simulation", settings)
        invalid = settings.__class__(
            settings.temperature,
            settings.total_steps,
            settings.statistics_frequency,
            settings.random_seed,
            "not-a-number",
            settings.replica_number,
            settings.energy_gap,
        )
        with self.assertRaisesRegex(ValueError, "closure sigma"):
            generate_mcmc_input(preset, profile, "structure.pdb", "simulation", invalid)


if __name__ == "__main__":
    unittest.main()
