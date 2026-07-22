# Release validation

Validation was performed on 2026-07-21 on Apple Silicon macOS. These checks
establish software/runtime compatibility; they do not establish scientific
convergence for a production study.

## Automated checks

- Python compilation completed without errors.
- The complete core test suite passed on the release source and extracted ZIP.
- The Plugin Manager ZIP contained all assets, no C/C++ source or object files,
  and executable mode metadata for both runtimes.
- The Qt plugin loaded under PyMOL's bundled Python 3.10 / Qt 5.15 runtime.

## Runtime and force-field checks

Each check used a prepared all-atom structure, the exact six selected files in
`forcefield/`, and a generated initialization input.

| Runtime | Profile | Result |
|---|---|---|
| MOSAICS 3.9.1 | bsc1/OL3 standard | Completed |
| MOSAICS 3.9.1 | bsc0 standard | Completed |
| MOSAICS 3.9.1 | OL15/OL3 standard | Completed |
| Experimental validated stack | All nine selectable profiles | Completed |

The measured 3.9.1 incompatibilities with terminal, OL21/OL24, and ff14SB
decks are encoded in the selector and validator.

Fixtures covered a 1EFS DNA/RNA heteroduplex, true-terminal DNA ACGT, and the
1HHK MHC/peptide ff14SB system with three disulfides.

## Preset checks

Every preset was parsed and executed with a shortened smoke setting. The full
BFGS local-minimum path also completed. One-step parallel-tempering tests only
verify initialization and proposal wiring; the displayed production lengths
must still be run and assessed for move/exchange acceptance and convergence.

The graphical `region.data` path completed a two-step CBLC run with residue
centers and paired residues.

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
