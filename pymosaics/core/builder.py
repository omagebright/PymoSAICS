"""Generate explicit, reviewable MOSAICS input files from named presets."""

from dataclasses import dataclass
import hashlib
import math
import os
from pathlib import Path
import shutil
import tempfile

from .catalog import AnalysisPreset, FORCEFIELD_ROOT, ForceFieldProfile


@dataclass(frozen=True)
class InputSettings:
    temperature: float
    total_steps: int
    statistics_frequency: int
    random_seed: int
    closure_sigma: str
    replica_number: int = 0
    energy_gap: float = 1.1


def default_settings(preset: AnalysisPreset) -> InputSettings:
    return InputSettings(
        temperature=preset.temperature,
        total_steps=preset.total_steps,
        statistics_frequency=preset.statistics_frequency,
        random_seed=-7143580450,
        closure_sigma=preset.closure_sigma,
        replica_number=max(0, preset.replica_count - 1),
        energy_gap=(
            (preset.top_temperature / preset.temperature) ** (1.0 / (preset.replica_count - 1))
            if preset.replica_count > 1 and preset.temperature > 0
            else 1.1
        ),
    )


def _asset_reference(relative: str) -> str:
    return "forcefield/" + Path(relative).name


def validate_profile_files(profile: ForceFieldProfile) -> None:
    missing = [path for path in profile.all_paths() if not path.is_file()]
    if missing:
        raise ValueError("force-field profile is incomplete: {}".format(", ".join(map(str, missing))))


def stage_force_field(project: Path, profile: ForceFieldProfile) -> Path:
    """Copy one selected six-file profile into a transparent project directory."""

    validate_profile_files(profile)
    destination = project.expanduser().resolve() / "forcefield"
    destination.mkdir(parents=True, exist_ok=True)
    manifest = [
        "PymoSAICS force-field profile",
        "identifier: {}".format(profile.identifier),
        "label: {}".format(profile.label),
        "validation: {}".format(profile.validation),
        "",
        "files:",
    ]
    for source in profile.all_paths():
        target = destination / source.name
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(
                dir=str(destination), prefix=".pymosaics-", suffix=".tmp", delete=False
            ) as handle:
                temporary = Path(handle.name)
            shutil.copy2(str(source), str(temporary))
            os.replace(str(temporary), str(target))
        finally:
            if temporary is not None and temporary.exists():
                temporary.unlink()
        digest = hashlib.sha256(target.read_bytes()).hexdigest()
        manifest.append("{}  {}".format(digest, target.name))
    manifest_target = destination / "PROFILE.txt"
    temporary_manifest = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=str(destination),
            prefix=".pymosaics-profile-",
            suffix=".tmp",
            delete=False,
        ) as handle:
            handle.write("\n".join(manifest) + "\n")
            temporary_manifest = Path(handle.name)
        os.replace(str(temporary_manifest), str(manifest_target))
    finally:
        if temporary_manifest is not None and temporary_manifest.exists():
            temporary_manifest.unlink()
    return destination


def generate_mcmc_input(
    preset: AnalysisPreset,
    force_field: ForceFieldProfile,
    structure_filename: str,
    output_stem: str,
    settings: InputSettings,
    region_filename: str = "",
) -> str:
    """Return a complete MOSAICS input whose scientific choices remain visible."""

    if preset.chemistry not in ("any", force_field.chemistry):
        raise ValueError(
            "preset {} is not compatible with {}".format(preset.label, force_field.label)
        )
    if settings.total_steps < 0 or settings.statistics_frequency < 1:
        raise ValueError("steps must be non-negative and statistics frequency must be positive")
    if settings.temperature <= 0 or not math.isfinite(settings.temperature):
        raise ValueError("temperature must be a positive finite number")
    if settings.replica_number < 0 or settings.energy_gap < 1 or not math.isfinite(settings.energy_gap):
        raise ValueError("replica count and temperature-ladder energy gap are invalid")
    try:
        closure_sigma = float(settings.closure_sigma)
    except (TypeError, ValueError) as exc:
        raise ValueError("closure sigma must be numeric") from exc
    if closure_sigma < 0 or not math.isfinite(closure_sigma):
        raise ValueError("closure sigma must be a non-negative finite number")

    def validate_entry_value(label: str, value: str, allow_directory: bool = False) -> None:
        if not value.strip():
            raise ValueError("a {} is required".format(label))
        if any(character in value for character in "{}\r\n"):
            raise ValueError("{} contains a character that is not valid in MOSAICS input".format(label))
        if not allow_directory and ("/" in value or "\\" in value):
            raise ValueError("{} must be a filename, not a path".format(label))

    validate_entry_value("structure filename", structure_filename)
    validate_entry_value("output stem", output_stem)
    if region_filename:
        validate_entry_value("region filename", region_filename, allow_directory=True)

    ddd_values = {
        "off": ("1.0", "1.0", "0.4", "0.5", "6.0"),
        "DD0S": ("80.0", "4.0", "0.4", "0.5", "6.0"),
    }
    ddd_d, ddd_d0, ddd_s, ddd_c, ddd_e = ddd_values[preset.dielectric]

    def entry(name: str, value) -> str:
        return "   \\{}{{{}}}".format(name, value)

    general = [
        "~sim_gen_def[",
        entry("simulation_typ", preset.simulation_type),
        entry("minimize_tol", preset.minimize_tolerance),
        entry("minimize_type", preset.minimize_type),
        entry("minimize_report", preset.minimize_report),
        entry("energy_report", 2),
        entry("prop_type", preset.proposal_type),
        entry("prop_tors_sig", "1.e-5"),
        entry("prop_trans_sig", 0),
        entry("prop_rot_sig", 0),
        entry("prop_tors_type", preset.torsion_type),
        entry("prop_clos_sig", settings.closure_sigma),
        entry("replica_number", settings.replica_number),
        entry("prob_eemc_jump", 0.15),
        entry("eemc_disk_size", 10),
        entry("energy_gap", "{:.8g}".format(settings.energy_gap)),
        entry("total_step_mc", settings.total_steps),
        entry("local_step_md", 1),
        entry("time_step_md", 0.4),
        entry("statistics_freq", settings.statistics_frequency),
        entry("burn_in_B", 10),
        entry("burn_in_N", 10),
        entry("write_energy_unit", "kcal"),
        entry("temperature", "{:g}".format(settings.temperature)),
        entry("inter_list", "none"),
        entry("rinter_switch_length", 0),
        entry("rinter_exclude_length", 1000),
        entry("random_seed", settings.random_seed),
        entry("EEMC_Emin", -1.0),
        entry("EEMC_Emax", 0.0),
        "]",
    ]
    if preset.minimize_type == "stsamc":
        general[6:6] = [
            entry("stsamc_type", "trigonom"),
            entry("stsamc_period", preset.stsamc_period),
            entry("stsamc_ampl", "{:g}".format(preset.stsamc_amplitude)),
            entry("stsamc_shift", 0),
        ]
    molecule = [
        "~sim_mol_def[",
        entry("system_def", "residue"),
        entry("implicit_solvent", "off"),
        entry("ddd", preset.dielectric),
        entry("ddd_D", ddd_d),
        entry("ddd_D0", ddd_d0),
        entry("ddd_S", ddd_s),
        entry("ddd_c", ddd_c),
        entry("ddd_e", ddd_e),
        entry("neutralize", preset.neutralize),
        entry("mol_parm_file", _asset_reference(force_field.rtf)),
        entry("bond_database_file", _asset_reference(force_field.bond)),
        entry("bend_database_file", _asset_reference(force_field.bend)),
        entry("tors_database_file", _asset_reference(force_field.torsion)),
        entry("onfo_database_file", _asset_reference(force_field.one_four)),
        entry("inter_database_file", _asset_reference(force_field.nonbonded)),
    ]
    if region_filename:
        molecule.append(entry("region_database_file", region_filename))
    molecule.extend(
        [
            entry("pos_init_file", structure_filename),
            entry("pos_out_file", output_stem + ".pos_out.pdb"),
            entry("atom_pos_file", output_stem + ".trajectory.pdb"),
            entry("tors_pos_file", output_stem + ".torsions.dat"),
            entry("epot_file", output_stem + ".potential_energy.dat"),
            entry("einter_file", output_stem + ".interaction_energy.dat"),
            entry("hessian_file", output_stem + ".hessian.dat"),
            entry("eighess_file", output_stem + ".eighess.dat"),
            "]",
        ]
    )
    return "\n".join(general) + "\n\n" + "\n".join(molecule) + "\n"


def write_mcmc_input(path: Path, content: str) -> Path:
    path = path.expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return path
