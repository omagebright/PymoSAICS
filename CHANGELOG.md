# Changelog

## 0.3.0 — 2026-07-22

- Added lossless MOSAICS input parsing and surgical scalar-setting updates that preserve repeated, unknown, region, and cryo-EM directives.
- Added recursive project discovery and safe portable import of foreign absolute paths into a managed `mcmc.input`, with ambiguity detection and preservation of the source deck.
- Added the historical KB_3pt topology/parameter profile to the stable and experimental runtime matrices.
- Added automatic canonical-protein conversion to CA/O/CMA coordinates, deterministic chain relabeling and residue numbering, an auditable mapping table, and explicit handling of incomplete residues.
- Added conservative whole-chain segment-region generation and a temperature-modulated three-point natural-move preset.
- Added direct UI controls to load an existing input into Build, apply visible values without regeneration, save reviewed changes, and make a project portable on demand.
- Completed one-step MOSAICS 3.9.1 runs for both Tom's imported 7QPJ deck and a fresh 7QPJ project generated from the RCSB PDB.

## 0.2.4 — 2026-07-22

- Moved new MOSAICS run logs from hidden application state into the visible project `logs/` directory while retaining discovery of legacy logs.
- Reloaded the complete persisted log into the Run tab after execution and removed the text-viewer size limit for logs.
- Automatically discovered `.input` and `.inp` files whenever a project directory is selected, with an explicit `mcmc.input` fallback for empty projects.
- Rebuilt the Analysis pages into responsive energy/acceptance, structural-map, and files/logs workspaces.
- Corrected Qt mnemonic rendering, clipped landscape controls, truncated acceptance headers, and minimum-window geometry regressions.

## 0.2.3 — 2026-07-22

- Rebuilt the Run tab so validation and execution-plan diagnostics cannot overlap subsequent controls.
- Placed diagnostics side by side, grouped all run actions, and gave live MOSAICS output the expanding workspace.
- Added minimum-window geometry assertions for diagnostics, actions, automatic loading, and live output.

## 0.2.2 — 2026-07-22

- Replaced host-native combo popups with readable, platform-independent Qt popup views.
- Added optional PDB2PQR/PROPKA preparation for heavy-atom proteins, with AMBER naming, explicit pH, retained PQR/log provenance, and mandatory ff14SB revalidation.
- Replaced the region dialog with a scientific workbench for membership, required rotation centers, non-overlapping residue pairs, WP2 move-width presets, units, live validation, and PyMOL selection/visualization.
- Removed the unsupported dependent residue-region choice and made single-region `superimpose` propagation explicit in generated inputs.
- Reproduced the previous 1A6Z failure path successfully: 6,021 prepared atoms, 371 residues, three disulfides, and zero ff14SB topology mismatches.

## 0.2.1 — 2026-07-22

- Rebuilt the Qt visual system with deterministic dark surfaces, readable labels, consistent spacing, and native-theme-independent plot canvases.
- Fixed clipped Build controls at the minimum supported window size and made Setup vertically scrollable on smaller displays.
- Added synchronized force-field profile selectors to Build and Setup.
- Made all nine bundled modern, terminal, protein, and AMBER99-based legacy profiles visible as complete sets.
- Added an exact preview of the six force-field directives populated by the selected profile.
- Added GUI checks for profile synchronization, complete directive generation, viewport width, and host-palette isolation.

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
