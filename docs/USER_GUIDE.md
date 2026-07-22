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

For a heavy-atom protein, select **Automatic AMBER preparation** and set the
protonation pH. PymoSAICS runs PDB2PQR with AMBER output names and PROPKA,
retains the selected input, PQR, converted PDB, and log under
`.pymosaics/protein-preparation/`, then applies disulfides and revalidates every
residue against ff14SB. **Strict mode** is for an already prepared all-atom PDB.
PyMOL's generic hydrogen addition is not an ff14SB preparation method.

## Region workbench

The editor creates one independent residue-level region:

- **Move** includes a residue; **Center** selects a whole-region rotation pivot.
  MOSAICS requires at least one center.
- **Residue pairs** are coupled units such as nucleic-acid base pairs. Each
  residue can occur in only one pair.
- **Whole region**, **Free residues**, and **Residue pairs** have separate
  translation widths in Å and rotation widths in radians.
- **WP2 balanced pilot** enables documented free-residue and pair widths;
  **WP2 paired-residue motion** holds free residues at zero.

Use **Use current PyMOL selection** to populate membership. **Show region in
PyMOL** displays members in cyan, paired residues in magenta, and rotation
centers in yellow. The generated `region.data` remains visible and invalid
configurations cannot be accepted. Run a short pilot and tune proposal widths
from acceptance before production.

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
