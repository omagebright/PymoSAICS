# Release validation

Validation was performed on 2026-07-21 and repeated for the 0.3.0 interface on
2026-07-22 on Apple Silicon macOS. These checks
establish software/runtime compatibility; they do not establish scientific
convergence for a production study.

## Automated checks

- Python compilation completed without errors.
- The complete core test suite passed on the release source and extracted ZIP.
- The Plugin Manager ZIP contained all assets, no C/C++ source or object files,
  and executable mode metadata for both runtimes.
- The Qt plugin loaded under PyMOL's bundled Python 3.10 / Qt 5.15 runtime.
- All tabs rendered without native light-background leakage under forced light
  and dark host palettes; the minimum-width Build form retained every control.
- Build and Setup exposed the same ten profiles. Selecting each profile
  populated all six corresponding force-field directives in `mcmc.input`.
- Every combo popup used an explicit dark Qt view with readable active,
  selected, and disabled text under the PyMOL runtime.
- The region workbench generated a valid center by default, rejected empty or
  overlapping-pair definitions, and wrote explicit single-region propagation.
- The Run tab retained non-overlapping validation, execution-plan, action,
  automatic-loading, and live-output regions at both 1092×874 and 900×710.
- Project-directory changes discovered all `.input` and `.inp` files, preferred
  an existing `mcmc.input`, and supplied a planned `mcmc.input` when none existed.
- New logs were written visibly under `logs/`; the complete persisted bytes were
  reproduced in the Run tab, and legacy hidden logs remained discoverable.
- All three Analysis pages rendered without clipping at 900×710: energy and
  acceptance panels, structural-map controls/plot/representatives, and files/logs.

## Runtime and force-field checks

Each check used a prepared all-atom structure, the exact six selected files in
`forcefield/`, and a generated initialization input.

| Runtime | Profile | Result |
|---|---|---|
| MOSAICS 3.9.1 | bsc1/OL3 standard | Completed |
| MOSAICS 3.9.1 | bsc0 standard | Completed |
| MOSAICS 3.9.1 | OL15/OL3 standard | Completed |
| MOSAICS 3.9.1 | KB_3pt | Completed |
| Experimental validated stack | All ten selectable profiles | Completed |

The measured 3.9.1 incompatibilities with terminal, OL21/OL24, and ff14SB
decks are encoded in the selector and validator.

Fixtures covered a 1EFS DNA/RNA heteroduplex, true-terminal DNA ACGT, and the
1HHK MHC/peptide ff14SB system with three disulfides.

The previously failing RCSB 1A6Z protein path was repeated through the GUI.
PDB2PQR/PROPKA at pH 7 prepared 6,021 atoms in 371 residues, preserved three
disulfides, and passed the ff14SB atom validator with zero mismatches. The GUI
retained the selected input, PQR, converted PDB, and preparation log.

## Preset checks

Every preset was parsed and executed with a shortened smoke setting. The full
BFGS local-minimum path also completed. One-step parallel-tempering tests only
verify initialization and proposal wiring; the displayed production lengths
must still be run and assessed for move/exchange acceptance and convergence.

The graphical `region.data` path completed a two-step CBLC run with residue
centers and paired residues.

## Portable and three-point workflow checks

Tom's original 7QPJ directory was copied to an isolated project and imported
without manual path editing. PymoSAICS made 14 compatibility/path rewrites,
preserved the original input, validated every reference, and completed a
one-step MOSAICS 3.9.1 run. The run produced the requested PDB, trajectory,
torsion, potential-energy, and interaction-energy files under `output/`, plus
the runtime-owned `sim_param.out` in the project directory.

Starting independently from the RCSB 7QPJ PDB, PymoSAICS generated 2,466 sites
for 822 residues across five relabeled chains, exactly matching the site and
residue counts of Tom's supplied three-point structure. All 822 CA and 822 O
coordinates matched Tom's file exactly. One RCSB residue lacking carbonyl O was
omitted and recorded in `structure.mapping.tsv`. A generated whole-chain region
deck then completed a one-step MOSAICS 3.9.1 KB_3pt run.

The generated CMA coordinates use a transparent geometric side-chain centroid
and are not asserted to reproduce Tom's separate side-chain preparation. Tom's
exact CMA coordinates, STRIDE boundaries, and 17 custom regions are preserved
when his original project is imported. Likewise, a fetched Mm-cpn PDB cannot
replace missing historical 2D image, orientation, or three-level region files.

## End-to-end PyMOL check

The public GUI completed this sequence:

1. load a structure and synchronize a PyMOL coordinate edit;
2. prepare and strictly validate the PDB/RTF pairing;
3. write a short-path force-field deck and visible `mcmc.input`;
4. launch MOSAICS through `QProcess` and capture its persistent log;
5. discover energy and output files;
6. load final and trajectory PDBs into PyMOL;
7. build an RMSD landscape with matched energy coloring; and
8. switch PyMOL to a selected representative frame.

An HTTPS fetch of RCSB entry 1BNA also completed and produced the expected two
chains. RCSB structures without all required hydrogens are intentionally
stopped by topology validation before execution.
