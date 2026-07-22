# PymoSAICS

PymoSAICS is a transparent PyMOL workbench for preparing, running, and
inspecting [MOSAICS](https://www.cs.ox.ac.uk/mosaics/) simulations. Structures,
headers, region definitions, force-field files, `mcmc.input`, commands, logs,
and outputs remain visible to the user.

## Install in PyMOL

1. Open **Plugin → Plugin Manager → Install New Plugin**.
2. Paste this URL into **Install from PyMOLWiki or any URL**:

   ```text
   https://github.com/omagebright/PymoSAICS/releases/download/v0.4.1/PymoSAICS-0.4.1.zip
   ```

3. Choose **Install**, restart PyMOL, then open **Plugin → PymoSAICS**.

## Included runtimes

The package includes two authorized compiled executables for Apple-Silicon
macOS. It does not include MOSAICS source code.

| Runtime | Selectable force-field profiles |
|---|---|
| MOSAICS 3.9.1 | bsc1/OL3 standard, bsc0 standard, OL15/OL3 standard, KB_3pt |
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
- true-terminal variants of OL15/OL3, OL21/OL3, and OL24/OL3;
- AMBER ff14SB for proteins, including the validated disulfide workflow; and
- the historical MOSAICS KB_3pt protein/protein-complex topology and potential.

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

**Run MOSAICS** executes the selected input without silently rebuilding it.
During a run, the reported Monte Carlo step and `total_step_mc` drive a live
percentage, elapsed time, measured step rate, countdown, and estimated finish
time. The estimate begins after MOSAICS reports its first step and updates as
new reports arrive.

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

### Build DNA and RNA directly

Under **Build → Build DNA / RNA**, choose a single DNA or RNA nucleotide/strand,
a DNA or RNA duplex, or a DNA:RNA hybrid. Enter strand 1 in 5′→3′ order;
PymoSAICS validates the optional second strand or constructs its antiparallel
Watson–Crick reverse complement. A- or B-form geometry is explicit for each
strand. Mixed-form constructions are labeled as hypotheses that require
relaxation.

The builder uses PyMOL's nucleic-acid geometry, creates and correctly names all
hydrogens, applies the required single-chain or true-terminal chemistry, and
then performs an exact residue-by-residue check against the selected MOSAICS
RTF before reporting success. Measured sugar pseudorotation phases are shown
immediately and colored in PyMOL (cyan A-like/C3′-endo, orange
B-like/C2′-endo).

### Existing and historical MOSAICS projects

Point PymoSAICS at a project directory or browse to an existing `.input` or
`.inp` file. Input discovery is recursive. When a single deck contains foreign
absolute paths and every referenced basename has one unambiguous local match,
PymoSAICS automatically creates a managed `mcmc.input` using
`${PROJECT_DIR}`. The source deck is preserved. Output paths are placed under
`output/`; the unsupported historical `param_out_file` directive is removed
because MOSAICS writes `sim_param.out` in the working directory.

Selecting a Run input loads its actual temperature, step count, statistics
frequency, seed, closure width, replica count, energy gap, and output prefix
into Build. **Apply visible settings** performs a surgical update of only those
existing directives. Repeated `energy_term` entries, cryo-EM options, comments,
and every unknown scientific option remain unchanged. Review the complete text,
then choose **Save input**. **Make portable** performs the same path import on
demand and reports ambiguous or missing files instead of guessing.

### Build a KB_3pt starting project from a PDB

Select **MOSAICS KB_3pt** and **Three-point protein natural moves**, then load a
local PDB, fetch an RCSB identifier, or use a PyMOL object. Preparation writes:

- `structure.pdb` with CA, carbonyl O, and a geometric centroid of the
  side-chain heavy atoms for each canonical protein residue;
- `structure.mapping.tsv`, including chain relabeling, residue renumbering,
  centroid provenance, and any incomplete residue omitted because CA or O was
  absent;
- the authentic `top_3pt_prot_na.rtf` and `par_3pt_prot_na.prm`; and
- one explicit segment region per chain as a conservative, editable starting
  hierarchy.

Glycine CMA is placed 0.01 Å from CA to avoid coincident sites, following the
historical convention. A fetched PDB is enough to create and execute a valid
KB_3pt starter project, but it cannot reveal experiment-specific STRIDE
boundaries, domain hierarchy, cryo-EM images, orientations, or energy weights.
Import the original deck and associated files when exact reproduction is the
goal—for example Tom's 7QPJ regions or the three Mm-cpn refinement levels.

## Presets

The interface provides reviewable starting points for:

- initial energy and force-field checks;
- BFGS local minimization;
- simulated-tempering minimum searches;
- regular and successive CBLC regressions;
- protein side-chain natural moves;
- three-point protein/protein-complex natural moves; and
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
structures are available under **Analysis**. The trajectory workspace aligns
frames to frame 1, plots RMSD, ranks start-to-end residue changes, measures DNA
and RNA sugar-pucker phases, and flags invariant atoms at 5′/3′ (or protein
N/C) termini. It can play an aligned PyMOL movie, compare cyan start and
magenta end structures, and color final A/B-like sugar states. The input and
protocol page shows the full current deck and records comparable evidence for
single-replica, fluctuating-temperature, and parallel-tempering runs. An
acceptance ratio from 0.20 through 0.50 is labeled as the current target; this
label is a tuning aid, not evidence of convergence.

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
