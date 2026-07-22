# PymoSAICS

PymoSAICS is a transparent PyMOL workbench for preparing, running, and
inspecting [MOSAICS](https://www.cs.ox.ac.uk/mosaics/) simulations. Structures,
headers, region definitions, force-field files, `mcmc.input`, commands, logs,
and outputs remain visible to the user.

## Install in PyMOL

1. Open **Plugin → Plugin Manager → Install New Plugin**.
2. Paste this URL into **Install from PyMOLWiki or any URL**:

   ```text
   https://github.com/omagebright/PymoSAICS/releases/download/v0.2.4/PymoSAICS-0.2.4.zip
   ```

3. Choose **Install**, restart PyMOL, then open **Plugin → PymoSAICS**.

## Included runtimes

The package includes two authorized compiled executables for Apple-Silicon
macOS. It does not include MOSAICS source code.

| Runtime | Selectable force-field profiles |
|---|---|
| MOSAICS 3.9.1 | bsc1/OL3 standard, bsc0 standard, OL15/OL3 standard |
| Experimental validated stack | Every bundled profile, including terminal OL15/21/24 and ff14SB |
| Custom executable | User-selected Windows, macOS, or Linux build; compatibility must be validated |

Windows, Linux, and Intel-macOS users should obtain a compatible executable
from the [official MOSAICS Downloads page](https://www.cs.ox.ac.uk/mosaics/Downloads.php),
then select **Custom executable** in **Setup**.

## Force-field profiles

Choose one complete profile under **Build → Runtime, force field, and analysis**
or **Setup → Force-field profile**. Both selectors remain synchronized. One
selection writes the topology, bond, angle, torsion/improper, 1–4 (`onfo`),
and nonbonded database directives into the visible `mcmc.input`; individual
files do not need to be assembled manually.

The bundled profiles are:

- AMBER99 + parmbsc1/OL3 for DNA/RNA;
- AMBER99 + parmbsc0 as the legacy DNA/RNA comparator;
- AMBER OL15/OL3, OL21/OL3, and OL24/OL3 for DNA/RNA;
- true-terminal variants of OL15/OL3, OL21/OL3, and OL24/OL3; and
- AMBER ff14SB for proteins, including the validated disulfide workflow.

Profiles that have not passed compatibility testing with the selected runtime
are disabled rather than silently mixed with an incompatible executable.

## Start a project

1. Choose a project directory.
2. Load a local PDB, fetch a four-character RCSB identifier, or use a live
   PyMOL object.
3. Select the model, chains, CBLC/SCBLC header, runtime, force field, and
   analysis preset.
4. If required, define a region graphically or confirm detected disulfides.
5. Generate and read the visible `mcmc.input`.
6. Choose **Prepare project**, inspect validation, then **Run MOSAICS**.
7. Follow the live log and inspect results under **Analysis**.

PyMOL coordinate edits are synchronized before preparation. The prepared PDB
is checked residue-by-residue against the selected all-atom RTF; missing
hydrogens or incorrect names stop the run with explicit diagnostics.

For an RCSB or other heavy-atom protein, choose **Automatic AMBER preparation**.
PymoSAICS uses PDB2PQR with AMBER naming and PROPKA at the selected pH, keeps
the PQR and log, and then requires an exact ff14SB match. If PDB2PQR is not
found automatically, [install it using its official instructions](https://pdb2pqr.readthedocs.io/en/stable/getting.html)
and select its executable under **Setup → Protein preparation**.

The region workbench defines one residue-level natural-move region. It exposes
members, required rotation centers, non-overlapping residue pairs, documented
WP2 proposal-width presets, Å/radian units, a live `region.data` preview, and
PyMOL selection/visualization. Invalid regions cannot be accepted.

## Presets

The interface provides reviewable starting points for:

- initial energy and force-field checks;
- BFGS local minimization;
- simulated-tempering minimum searches;
- regular and successive CBLC regressions;
- protein side-chain natural moves; and
- nucleic-acid or protein parallel-tempering landscape pilots.

Every parameter remains editable. A finite minimization or sampling run cannot
guarantee the global minimum. The clickable landscape is a two-dimensional map
of aligned pairwise RMSD—not, by itself, a converged free-energy surface.

## Files and analysis

A prepared project contains short, parser-safe relative paths:

```text
project/
├── structure.pdb
├── mcmc.input
├── forcefield/        # the six selected files plus checksums
├── region/            # only when enabled
├── analysis/          # exported structural-map coordinates
└── logs/              # visible, persistent, complete run logs
```

Opening or entering a project directory automatically discovers its `.input`
and `.inp` files in the Run selector. An existing `mcmc.input` is preferred;
when no input exists, PymoSAICS starts with the planned `mcmc.input` path and
the current visible defaults.

Before running, `mcmc.input` and `structure.pdb` can be opened and edited as
text inside PymoSAICS. After running, complete logs are visible under `logs/`
and are reloaded in full in the Run tab and text viewer; PDB files and
trajectories load directly into PyMOL.
Energy series, natural-move acceptance, and clickable representative
structures are available under **Analysis**.

## Credits and terms

MOSAICS was created by **Peter Minary**. Cite the relevant work listed on the
[official publications page](https://www.cs.ox.ac.uk/mosaics/Publications.php)
and follow the official licensing terms.

The original basic PyMOSAICS GUI is credited to **Konrad Krawczyk**. This
Python 3 / Qt redesign and current PymoSAICS release are by
**Folorunsho Bright Omage**.

PymoSAICS source code is MIT licensed. Bundled MOSAICS executables and
force-field data are separate works with their own provenance and terms; see
[`pymosaics/assets/ASSET_NOTICE.md`](pymosaics/assets/ASSET_NOTICE.md).

## Test and package

```bash
python -m compileall -q pymosaics
python -m unittest discover -v
python scripts/build_release.py
```

Report reproducible plugin problems through the
[issue tracker](https://github.com/omagebright/PymoSAICS/issues).

The release evidence and tested compatibility matrix are recorded in
[`docs/VALIDATION.md`](docs/VALIDATION.md).
