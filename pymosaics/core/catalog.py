"""Bundled MOSAICS runtimes, force fields, and transparent analysis presets."""

import os
import platform
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
ASSET_ROOT = PACKAGE_ROOT / "assets"
FORCEFIELD_ROOT = ASSET_ROOT / "forcefields"


@dataclass(frozen=True)
class RuntimeProfile:
    identifier: str
    label: str
    description: str
    version: str
    source_commit: str
    relative_executable: Optional[str]
    supported_platform: Optional[str]
    supported_machine: Optional[str]

    def executable(self) -> Optional[Path]:
        if self.relative_executable is None:
            return None
        return ASSET_ROOT / self.relative_executable

    def available(self, platform_name: Optional[str] = None, machine: Optional[str] = None) -> bool:
        if self.relative_executable is None:
            return True
        platform_name = platform_name or sys.platform
        machine = (machine or host_machine(platform_name)).lower()
        machine = "arm64" if machine in ("arm64", "aarch64") else machine
        return (
            platform_name == self.supported_platform
            and machine == self.supported_machine
            and bool(self.executable() and self.executable().is_file())
        )


def host_machine(platform_name: Optional[str] = None) -> str:
    """Return host architecture, including Apple Silicon behind Rosetta PyMOL."""

    platform_name = platform_name or sys.platform
    machine = platform.machine().lower()
    if platform_name == "darwin" and machine == "x86_64":
        try:
            result = subprocess.run(
                ("/usr/sbin/sysctl", "-n", "hw.optional.arm64"),
                check=False,
                capture_output=True,
                text=True,
                timeout=2,
            )
        except (OSError, subprocess.SubprocessError):
            pass
        else:
            if result.returncode == 0 and result.stdout.strip() == "1":
                return "arm64"
    return machine


@dataclass(frozen=True)
class ForceFieldProfile:
    identifier: str
    label: str
    description: str
    chemistry: str
    topology_profile: str
    rtf: str
    bond: str
    bend: str
    torsion: str
    one_four: str
    nonbonded: str
    validation: str

    def path(self, relative: str) -> Path:
        return FORCEFIELD_ROOT / relative

    def all_paths(self) -> Tuple[Path, ...]:
        return tuple(
            self.path(value)
            for value in (self.rtf, self.bond, self.bend, self.torsion, self.one_four, self.nonbonded)
        )


@dataclass(frozen=True)
class AnalysisPreset:
    identifier: str
    label: str
    description: str
    chemistry: str
    proposal_type: str
    torsion_type: str
    total_steps: int
    statistics_frequency: int
    temperature: float
    closure_sigma: str
    dielectric: str
    neutralize: str
    recommended_header: str
    simulation_type: str = "PT"
    minimize_type: str = "bfgs"
    minimize_tolerance: str = "1.e-3"
    minimize_report: int = 0
    replica_count: int = 1
    top_temperature: float = 300.0
    stsamc_period: int = 0
    stsamc_amplitude: float = 0.0


RUNTIME_PROFILES = (
    RuntimeProfile(
        identifier="mosaics-3.9.1",
        label="MOSAICS 3.9.1 (stable)",
        description="Reference MOSAICS 3.9.1 executable used for backward-comparison runs.",
        version="3.9.1",
        source_commit="legacy-version.3.9.1_bgq",
        relative_executable="runtimes/darwin-arm64/mosaics-3.9.1",
        supported_platform="darwin",
        supported_machine="arm64",
    ),
    RuntimeProfile(
        identifier="mosaics-experimental-2026-07-21",
        label="MOSAICS experimental (validated stack, 2026-07-21)",
        description=(
            "Validated Apple-Silicon build combining general ff14SB/disulfides, CMAP, "
            "and standards-compliant PDB output."
        ),
        version="experimental-2026-07-21",
        source_commit="94120f4+6341e5e+370a8b5",
        relative_executable="runtimes/darwin-arm64/mosaics-experimental-2026-07-21",
        supported_platform="darwin",
        supported_machine="arm64",
    ),
    RuntimeProfile(
        identifier="custom",
        label="Custom executable",
        description="Select a MOSAICS executable compiled for this computer.",
        version="custom",
        source_commit="user-supplied",
        relative_executable=None,
        supported_platform=None,
        supported_machine=None,
    ),
)


FORCE_FIELD_PROFILES = (
    ForceFieldProfile(
        "bsc1-ol3-standard",
        "AMBER99 + parmbsc1 / OL3 — DNA/RNA",
        "Validated baseline for canonical internalized DNA, RNA, and DNA/RNA hybrids.",
        "nucleic_acid",
        "standard",
        "panel/top_all99-bsc1_prot_hybrid_chidef.rtf",
        "panel/db_bsc1/mosaics_amber99-bsc1.bond",
        "panel/db_bsc1/mosaics_amber99-bsc1.bend",
        "panel/db_bsc1/mosaics_amber99-bsc1.tors_and_impr",
        "panel/db_bsc1/mosaics_amber99-bsc1.onfo",
        "panel/db_bsc1/mosaics_amber99-bsc1.vdw",
        "Validated baseline one-step regression.",
    ),
    ForceFieldProfile(
        "bs0-standard",
        "AMBER99 + parmbsc0 — legacy DNA/RNA",
        "Validated historical comparator for canonical internalized nucleic acids.",
        "nucleic_acid",
        "standard",
        "panel/top_all99-bs0_prot_hybrid_chidef.rtf",
        "panel/db_bs0/mosaics_amber99-bs0.bond",
        "panel/db_bs0/mosaics_amber99-bs0.bend",
        "panel/db_bs0/mosaics_amber99-bs0.tors_and_impr",
        "panel/db_bs0/mosaics_amber99-bs0.onfo",
        "panel/db_bs0/mosaics_amber99-bs0.vdw",
        "Validated baseline one-step regression.",
    ),
    ForceFieldProfile(
        "ol15-ol3-standard",
        "AMBER OL15 / OL3 — DNA/RNA",
        "Validated OL15 DNA plus OL3 RNA profile for internalized systems.",
        "nucleic_acid",
        "standard",
        "panel/top_openmm-ol15_ol3_hybrid_chidef.rtf",
        "panel/db_ol15_ol3/mosaics_openmm-ol15_ol3.bond",
        "panel/db_ol15_ol3/mosaics_openmm-ol15_ol3.bend",
        "panel/db_ol15_ol3/mosaics_openmm-ol15_ol3.tors_and_impr",
        "panel/db_ol15_ol3/mosaics_openmm-ol15_ol3.onfo",
        "panel/db_ol15_ol3/mosaics_openmm-ol15_ol3.vdw",
        "1EFS, reduced-system, and terminal validation passed.",
    ),
    ForceFieldProfile(
        "ol15-ol3-terminal",
        "AMBER OL15 / OL3 — DNA/RNA true termini",
        "Validated for canonical chains with explicit 5-prime and 3-prime terminal chemistry.",
        "nucleic_acid",
        "terminal",
        "panel/top_openmm-ol15_ol3_terminal_hybrid_chidef.rtf",
        "panel/db_ol15_ol3_terminal/mosaics_openmm-ol15_ol3.bond",
        "panel/db_ol15_ol3_terminal/mosaics_openmm-ol15_ol3.bend",
        "panel/db_ol15_ol3_terminal/mosaics_openmm-ol15_ol3.tors_and_impr",
        "panel/db_ol15_ol3_terminal/mosaics_openmm-ol15_ol3.onfo",
        "panel/db_ol15_ol3_terminal/mosaics_openmm-ol15_ol3.vdw",
        "Pure-DNA and pure-RNA terminal controls passed.",
    ),
    ForceFieldProfile(
        "ol21-ol3-standard",
        "AMBER OL21 / OL3 — DNA/RNA",
        "Validated current AMBER DNA/RNA combination for internalized heteroduplexes.",
        "nucleic_acid",
        "standard",
        "panel/top_openmm-ol21_ol3_hybrid_chidef.rtf",
        "panel/db_ol21_ol3/mosaics_openmm-ol21_ol3.bond",
        "panel/db_ol21_ol3/mosaics_openmm-ol21_ol3.bend",
        "panel/db_ol21_ol3/mosaics_openmm-ol21_ol3.tors_and_impr",
        "panel/db_ol21_ol3/mosaics_openmm-ol21_ol3.onfo",
        "panel/db_ol21_ol3/mosaics_openmm-ol21_ol3.vdw",
        "1EFS, reduced-system, and terminal validation passed.",
    ),
    ForceFieldProfile(
        "ol21-ol3-terminal",
        "AMBER OL21 / OL3 — DNA/RNA true termini",
        "Validated for canonical chains with explicit 5-prime and 3-prime terminal chemistry.",
        "nucleic_acid",
        "terminal",
        "panel/top_openmm-ol21_ol3_terminal_hybrid_chidef.rtf",
        "panel/db_ol21_ol3_terminal/mosaics_openmm-ol21_ol3.bond",
        "panel/db_ol21_ol3_terminal/mosaics_openmm-ol21_ol3.bend",
        "panel/db_ol21_ol3_terminal/mosaics_openmm-ol21_ol3.tors_and_impr",
        "panel/db_ol21_ol3_terminal/mosaics_openmm-ol21_ol3.onfo",
        "panel/db_ol21_ol3_terminal/mosaics_openmm-ol21_ol3.vdw",
        "Pure-DNA and pure-RNA terminal controls passed.",
    ),
    ForceFieldProfile(
        "ol24-ol3-standard",
        "AMBER OL24 / OL3 — DNA/RNA",
        "Validated OL24 DNA plus OL3 RNA profile for internalized heteroduplexes.",
        "nucleic_acid",
        "standard",
        "panel/top_amber-ol24_ol3_hybrid_chidef.rtf",
        "panel/db_ol24_ol3/mosaics_amber-ol24_ol3.bond",
        "panel/db_ol24_ol3/mosaics_amber-ol24_ol3.bend",
        "panel/db_ol24_ol3/mosaics_amber-ol24_ol3.tors_and_impr",
        "panel/db_ol24_ol3/mosaics_amber-ol24_ol3.onfo",
        "panel/db_ol24_ol3/mosaics_amber-ol24_ol3.vdw",
        "1EFS, reduced-system, and terminal validation passed.",
    ),
    ForceFieldProfile(
        "ol24-ol3-terminal",
        "AMBER OL24 / OL3 — DNA/RNA true termini",
        "Validated for canonical chains with explicit 5-prime and 3-prime terminal chemistry.",
        "nucleic_acid",
        "terminal",
        "panel/top_amber-ol24_ol3_terminal_hybrid_chidef.rtf",
        "panel/db_ol24_ol3_terminal/mosaics_amber-ol24_ol3.bond",
        "panel/db_ol24_ol3_terminal/mosaics_amber-ol24_ol3.bend",
        "panel/db_ol24_ol3_terminal/mosaics_amber-ol24_ol3.tors_and_impr",
        "panel/db_ol24_ol3_terminal/mosaics_amber-ol24_ol3.onfo",
        "panel/db_ol24_ol3_terminal/mosaics_amber-ol24_ol3.vdw",
        "Pure-DNA and pure-RNA terminal controls passed.",
    ),
    ForceFieldProfile(
        "ff14sb-protein",
        "AMBER ff14SB — protein",
        "Standalone canonical-protein profile; CMAP is not required for ff14SB.",
        "protein",
        "protein",
        "ff14sb/top_openmm-ff14sb_protein.rtf",
        "ff14sb/mosaics_openmm-ff14sb.bond",
        "ff14sb/mosaics_openmm-ff14sb.bend",
        "ff14sb/mosaics_openmm-ff14sb.tors_and_impr",
        "ff14sb/mosaics_openmm-ff14sb.onfo",
        "ff14sb/mosaics_openmm-ff14sb.vdw",
        "Validated with the 1HHK protein control and side-chain natural moves.",
    ),
)


STABLE_391_FORCE_FIELDS = {
    "bsc1-ol3-standard",
    "bs0-standard",
    "ol15-ol3-standard",
}


ANALYSIS_PRESETS = (
    AnalysisPreset(
        "single-point",
        "Single-point energy check",
        "Reads the structure and force field, reports the initial energy, and performs no Monte Carlo steps.",
        "any",
        "cart",
        "full",
        0,
        1,
        300.0,
        "0.001",
        "off",
        "off",
        "none",
        simulation_type="MIN",
        minimize_type="samc",
        minimize_report=1,
    ),
    AnalysisPreset(
        "local-minimum",
        "Local energy minimization (BFGS)",
        "Relaxes the supplied conformation to a nearby local minimum. This does not claim to find the global minimum.",
        "any",
        "cart",
        "full",
        1,
        1,
        300.0,
        "0.001",
        "off",
        "off",
        "none",
        simulation_type="MIN",
        minimize_type="bfgs",
        minimize_tolerance="1.e-7",
        minimize_report=2,
    ),
    AnalysisPreset(
        "stsamc-minimum-search",
        "Broad minimum search (simulated tempering)",
        "A 200,000-step torsional simulated-tempering/annealing search for multiple low-energy basins; no finite search guarantees the global minimum.",
        "nucleic_acid",
        "tors",
        "full",
        200000,
        2000,
        300.0,
        "0.001",
        "DD0S",
        "nucl",
        "regular",
        simulation_type="MIN",
        minimize_type="stsamc",
        minimize_tolerance="1.e-10",
        minimize_report=2,
        stsamc_period=50000,
        stsamc_amplitude=800.0,
    ),
    AnalysisPreset(
        "cblc-regression",
        "Regular CBLC regression",
        "10,000 torsional proposals for a short reproducibility and acceptance-rate comparison.",
        "nucleic_acid",
        "tors",
        "full",
        10000,
        200,
        300.0,
        "0.001",
        "DD0S",
        "nucl",
        "regular",
    ),
    AnalysisPreset(
        "scblc-regression",
        "Successive CBLC regression",
        "10,000 successive-CBLC proposals; uses the full-system closure width validated on 1EFS.",
        "nucleic_acid",
        "tors",
        "full",
        10000,
        200,
        300.0,
        "0.0001",
        "DD0S",
        "nucl",
        "successive",
    ),
    AnalysisPreset(
        "landscape-pilot",
        "Nucleic-acid landscape (adaptive PT)",
        "Parallel tempering with a visible geometric temperature ladder. The initial ladder is sized from the loaded structure and must be checked against move and exchange acceptance before conclusions are drawn.",
        "nucleic_acid",
        "tors",
        "full",
        1000000,
        10000,
        300.0,
        "0.0001",
        "DD0S",
        "nucl",
        "successive",
        replica_count=6,
        top_temperature=500.0,
    ),
    AnalysisPreset(
        "protein-landscape-pilot",
        "Protein landscape (adaptive PT)",
        "Parallel tempering over ff14SB side-chain natural moves with a visible geometric ladder. Validate move and exchange acceptance before interpreting basin populations.",
        "protein",
        "tors",
        "side_chain",
        1000000,
        10000,
        300.0,
        "0.001",
        "off",
        "off",
        "none",
        replica_count=6,
        top_temperature=529.0,
    ),
    AnalysisPreset(
        "protein-side-chain",
        "Protein side-chain natural moves",
        "300-step ff14SB side-chain smoke test matching the validated MHC workflow.",
        "protein",
        "tors",
        "side_chain",
        300,
        3,
        300.0,
        "0.001",
        "off",
        "off",
        "none",
    ),
)


def runtime_profile(identifier: str) -> RuntimeProfile:
    return next(item for item in RUNTIME_PROFILES if item.identifier == identifier)


def force_field_profile(identifier: str) -> ForceFieldProfile:
    return next(item for item in FORCE_FIELD_PROFILES if item.identifier == identifier)


def analysis_preset(identifier: str) -> AnalysisPreset:
    return next(item for item in ANALYSIS_PRESETS if item.identifier == identifier)


def runtime_supports_force_field(runtime_id: str, force_field_id: str) -> bool:
    if runtime_id == "mosaics-3.9.1":
        return force_field_id in STABLE_391_FORCE_FIELDS
    return runtime_id in ("mosaics-experimental-2026-07-21", "custom")


def make_bundled_executable_runnable(profile: RuntimeProfile) -> None:
    path = profile.executable()
    if path is not None and path.is_file() and os.name != "nt":
        path.chmod(path.stat().st_mode | 0o111)
