# PymoSAICS

PymoSAICS is a transparent PyMOL workbench for preparing, running, and
inspecting [MOSAICS](https://www.cs.ox.ac.uk/mosaics/) simulations. Structures,
headers, region definitions, force-field files, `mcmc.input`, commands, logs,
and outputs remain visible to the user.

## Install in PyMOL

1. Open **Plugin → Plugin Manager → Install New Plugin**.
2. Paste this URL into **Install from PyMOLWiki or any URL**:

   ```text
   https://github.com/omagebright/PymoSAICS/releases/download/v0.2.0/PymoSAICS-0.2.0.zip
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
└── .pymosaics/logs/   # persistent run logs
```

Before running, `mcmc.input` and `structure.pdb` can be opened and edited as
text inside PymoSAICS. After running, logs and text outputs can be opened from
the same interface; PDB files and trajectories load directly into PyMOL.
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
