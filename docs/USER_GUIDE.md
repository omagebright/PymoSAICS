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
validation result. Choosing or entering a project directory automatically
populates the input selector with every `.input` and `.inp` file. An existing
`mcmc.input` is preferred; if none exists, the planned path defaults to
`mcmc.input` for preparation from the visible Build settings.

MOSAICS is started without a shell. Output appears live and is retained under
the visible `logs/` directory. The complete persisted log is loaded back into
the Run tab after execution and can be opened without truncation in the text
viewer. Legacy `.pymosaics/logs/` records remain discoverable. **Stop** first
requests normal termination and then kills the process if it does not exit
within three seconds.

Use **View / edit** for `mcmc.input`, **View / edit input PDB** for the prepared
structure, and **View latest log** for the newest run record.

## Portable import and two-way input synchronization

Project scanning is recursive. If a directory contains one legacy input with
foreign absolute paths, PymoSAICS creates a portable `mcmc.input` when every
referenced filename has one unique match under that project. The original is
never overwritten. Use **Make portable** for an explicit import. Missing or
duplicate basename matches stop the import and are reported for manual choice.

Selecting an input loads its supported scalar values into Build without
regenerating the deck. **Apply visible settings** updates those existing
directives in the preview; **Save input** persists the reviewed text. Options
without a corresponding control—including repeated energy terms, segment
regions, and MOSAICS-EM settings—are retained exactly. Direct edits in the
text viewer are reloaded into Build when it closes.

## KB_3pt protein projects

The **MOSAICS KB_3pt** profile is compatible with the bundled stable 3.9.1 and
experimental runtimes. It stages the historical topology and parameter pair.
The **Three-point protein natural moves** preset generates a visible
temperature-modulated input with `cgres_model{KB_3pt}` and all five required
energy terms.

For an ordinary protein PDB, preparation selects CA and carbonyl O and computes
the geometric centroid of side-chain heavy atoms as CMA. Chains are relabeled
A, B, C… in selected PDB order and residues are renumbered from one; the exact
mapping and centroid method are written to `structure.mapping.tsv`. Glycine
uses CA plus a 0.01 Å x-offset. Non-polymer HETATM records are ignored.
Incomplete residues lacking CA or O are omitted and explicitly recorded rather
than silently repaired.

With regions enabled, PymoSAICS writes one whole-chain segment region per chain.
This is a runnable baseline, not an inferred biological hierarchy. Edit STRIDE
and `region.data` to encode domain boundaries, flexible closure loops, TCR CDRs,
Mm-cpn refinement levels, or other system-specific hypotheses. Exact historical
MOSAICS-EM reproduction additionally requires the original image,
`orientation.data`, region levels, and other experimental inputs.

## Analysis

- **Energy & acceptance** places the selected energy trace and complete-log
  acceptance counts in resizable, side-by-side scientific panels.
- **Structural landscape** aligns trajectory frames, computes pairwise RMSD,
  projects the distances into two dimensions, and selects representative
  frames. Clicking a point shows that frame in PyMOL.
- **Files & logs** lists project outputs. Text files open in PymoSAICS; PDB
  files and trajectories load directly into PyMOL.

The structural map is descriptive. It is not automatically a converged
free-energy surface, and representative frequency is not interpretable as an
equilibrium population without an adequate sampling and reweighting protocol.
