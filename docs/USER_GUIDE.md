# User guide

## Build

Choose a project directory and one structure source:

- **Local PDB** opens an existing file.
- **RCSB PDB** downloads a four-character RCSB entry over HTTPS.
- **Live PyMOL object** reads the selected object's current coordinates.

When live synchronization is enabled, edits made in PyMOL are read again
before **Prepare project**. Select the model and chains explicitly. PymoSAICS
then writes the displayed regular CBLC (`CBLC >...`), successive CBLC
(`CBLC ~...`), or no-header choice.

Select a compatible runtime, force field, and preset. All generated parameters
are shown in the editable `mcmc.input` preview. A graphical region is optional.

For proteins, SG–SG pairs at or below 2.5 Å are listed as disulfide candidates.
Selected cysteines are written as `CYX`, thiol hydrogen is removed, and the
pair can be shown as sticks in PyMOL. One-to-one pairs are selected by default;
ambiguous sulfur neighborhoods are displayed for manual choice but left
unchecked to avoid assigning one cysteine to multiple bonds.

## Preparation checks

Preparation creates `structure.pdb`, `mcmc.input`, and a short-path
`forcefield/` directory containing the selected six-file profile plus
checksums. Unrelated files already in that directory are preserved. If enabled,
`region/region.data` is also written.

The prepared PDB is compared residue-by-residue with the selected RTF. Missing,
unexpected, or duplicate atoms stop preparation. This is especially important
for RCSB X-ray structures, which commonly lack the hydrogens required by the
bundled all-atom profiles.

## Run

The Run tab shows the exact executable, input argument, working directory, and
validation result. MOSAICS is started without a shell. Output appears live and
is retained under `.pymosaics/logs/`. **Stop** first requests normal termination
and then kills the process if it does not exit within three seconds.

Use **View / edit** for `mcmc.input`, **View / edit input PDB** for the prepared
structure, and **View latest log** for the newest run record.

## Analysis

- **Energy & acceptance** plots discovered energy files and reports MOSAICS
  natural-move acceptance counts.
- **Structural landscape** aligns trajectory frames, computes pairwise RMSD,
  projects the distances into two dimensions, and selects representative
  frames. Clicking a point shows that frame in PyMOL.
- **Files & logs** lists project outputs. Text files open in PymoSAICS; PDB
  files and trajectories load directly into PyMOL.

The structural map is descriptive. It is not automatically a converged
free-energy surface, and representative frequency is not interpretable as an
equilibrium population without an adequate sampling and reweighting protocol.
