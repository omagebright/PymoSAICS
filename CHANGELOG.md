# Changelog

## 0.2.0 — 2026-07-21

- Added a project builder for local, RCSB, and live PyMOL structures.
- Added live coordinate synchronization, model/chain selection, CBLC/SCBLC headers, and graphical region files.
- Bundled authorized Apple-Silicon builds of MOSAICS 3.9.1 and the validated experimental stack.
- Added validated bsc0, bsc1, OL15/OL3, OL21/OL3, OL24/OL3, terminal, and ff14SB profiles.
- Added runtime/profile compatibility enforcement and strict PDB-to-RTF atom validation.
- Added protein disulfide discovery, CYX processing, and PyMOL visualization.
- Added visible presets for initialization, minimization, CBLC/SCBLC checks, protein side chains, and parallel-tempering landscape pilots.
- Added persistent logs, editable pre-run files, project-file inspection, automatic PDB/trajectory loading, energy plots, acceptance tables, and clickable RMSD landscape representatives.
- Added a purpose-built Qt interface and comprehensive cross-platform regression tests.

## 0.1.1 — 2026-07-21

- Corrected Windows test fixtures to use a native executable name.
- Normalized Windows temporary paths before portability assertions.
- Updated GitHub Actions to their Node.js 24-compatible releases.

## 0.1.0 — 2026-07-21

- Replaced the legacy Tk/Pmw interface with PyMOL's Qt compatibility layer.
- Added the modern PyMOL plugin entry point.
- Removed runtime downloading and 32-bit platform assumptions.
- Added an external-runtime setup and validation workflow.
- Added shell-free asynchronous execution, live logging, and cancellation.
- Added portable project and force-field placeholders.
- Added project input-reference validation and PyMOL output loading.
- Added standard per-user configuration paths.
- Added Windows, macOS, and Linux automated core tests.
