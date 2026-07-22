# Using the AMBER-family force-field decks in MOSAICS

This guide explains how to run a canonical DNA, RNA, or DNA/RNA hybrid system
in MOSAICS with the validated AMBER-family parameter decks: `parmbsc1`,
`parmbsc0`, `OL21/OL3`, `OL15/OL3`, and `OL24/OL3`. It is written for a user who has just
cloned the repository and has no local knowledge of how the decks were built.

For where the decks come from and their redistribution terms, see
[PROVENANCE.md](PROVENANCE.md). For the per-deck file layout and the list of
bug fixes baked into each deck, see [README.md](README.md).

## 1. Pick a force field

| Force field | DNA | RNA | Use it for |
|---|---|---|---|
| `parmbsc1` (`db_bsc1`) | parmbsc1 | parm99 + chi corr | general DNA/RNA hybrid baseline |
| `parmbsc0` (`db_bs0`) | parmbsc0 | parm99 | baseline comparator |
| `OL21/OL3` (`db_ol21_ol3`) | OL21 | OL3 | current AMBER DNA + RNA combination |
| `OL15/OL3` (`db_ol15_ol3`) | OL15 | OL3 | OL15 DNA + OL3 RNA combination |
| `OL24/OL3` (`db_ol24_ol3`) | OL24 | OL3 | OL24 DNA + OL3 RNA combination from official AMBER OL24 source |

All five are validated at the single-point energy level against AMBER `sander`
and OpenMM. The OL-family decks differ primarily in the DNA torsion set; the
RNA half (OL3) is shared.

## 2. Pick a profile: standard vs terminal

Each OL force field ships in two profiles because AMBER can resolve
different torsion terms for the same atom-type tuple
depending on whether a residue is internal or at a true chain terminus.

- **Standard profile** (`db_ol21_ol3`, `db_ol15_ol3`, `db_ol24_ol3`, and the
  matching hybrid RTF): use this for matched `1EFS`-style
  systems and any heteroduplex whose chain ends are treated as internalized
  (no explicit 5'-OH / 3'-OH terminal chemistry). This is the common case for
  the heteroduplex landscape work.

- **Terminal profile** (`db_ol21_ol3_terminal`, `db_ol15_ol3_terminal`,
  `db_ol24_ol3_terminal`, and the matching terminal RTF): use this for true
  pure DNA or pure RNA chains with explicit 5'-prime and 3'-prime terminal
  residues.

Do not mix: a true-terminal chain run under the standard profile (or vice
versa) will give the wrong torsion energy. If you are unsure, the terminal
profile is correct only when your input PDB actually contains terminal residue
chemistry (5'-OH / 3'-phosphate-free ends).

## 3. Supported and unsupported chemistry

Supported in this release:

- canonical DNA bases A, C, G, T;
- canonical RNA bases A, C, G, U;
- canonical DNA/RNA hybrids, pure DNA, pure RNA;
- bulges and mismatches, as long as every residue is a canonical base and the
  topology is valid.

Not supported (out of scope for these decks; do not expect correct energies):

- modified or noncanonical bases (5mC, 5hmC, phosphorothioate, etc.);
- ligands and protein-ligand complexes;
- proteins under the OL-family nucleic-acid decks;
- explicit solvent as mobile particles, and mobile ions;
- polarizable force fields;
- `CHARMM36+sHBfix` (the local `db_charmm36/` directory is excluded,
  not a release deck; see README.md).

## 4. Get your system into MOSAICS naming

MOSAICS reads a PDB whose residue and atom names match the hybrid RTF
(`top_openmm-ol*_hybrid_chidef.rtf`). The hybrid RTF uses the project's
DNA/RNA hybrid residue names (for example the hybrid DNA residues ADD/THD/CYD/GUD
and RNA ADE/URA/CYT/GUA) with the chi-definition atom naming.

The public one-command route is `tools/pdb_to_mosaics.py`. It is a pure Python
renamer for canonical DNA, RNA, and DNA/RNA hybrid PDB files. It does not call
`tleap`, ParmEd, OpenMM, GROMACS, AMBER `sander`, or any other MD engine.

For an internalized heteroduplex or a 1EFS-style system, use the standard
profile:

```bash
python tools/pdb_to_mosaics.py \
    --input my_hybrid_standard_names.pdb \
    --output my_hybrid_mosaics.pdb \
    --profile standard \
    --rtf params/top_openmm-ol21_ol3_hybrid_chidef.rtf \
    --report my_hybrid_mosaics_mapping.tsv
```

For a true-terminal pure DNA or pure RNA chain with explicit 5-prime and
3-prime end chemistry, use the terminal profile and the matching terminal RTF:

```bash
python tools/pdb_to_mosaics.py \
    --input my_terminal_dna.pdb \
    --output my_terminal_dna_mosaics.pdb \
    --profile terminal \
    --rtf params/top_openmm-ol21_ol3_terminal_hybrid_chidef.rtf \
    --report my_terminal_dna_mosaics_mapping.tsv
```

The renamer maps these canonical residue families:

| Input style | Standard output | Terminal output at 5-prime end | Terminal output at 3-prime end |
|---|---|---|---|
| DNA `DA/DC/DG/DT` and `DA5/DA3` variants | `ADD/CYD/GUD/THD` | `AD5/CD5/GD5/TD5` | `AD3/CD3/GD3/TD3` |
| RNA `A/C/G/U`, `RA/RC/RG/RU`, and terminal variants | `ADE/CYT/GUA/URA` | `AR5/CR5/GR5/UR5` | `AR3/CR3/GR3/UR3` |

It also applies the atom-name aliases used by the RTF: `OP1 -> O1P`,
`OP2 -> O2P`, RNA `H2' -> H2''`, RNA `HO2' -> H2'`, and thymine
`C7/H71/H72/H73 -> C5M/H51/H52/H53`. A computed `CBLC >...` chain line is
prepended by default. Use `--no-cblc` only if your downstream workflow has a
specific reason to omit it.

The optional `--rtf` check is recommended. By default it fails fast if the
renamed PDB contains atoms that are not defined by the selected MOSAICS residue
templates. It allows missing atoms because some validated inputs omit hydrogens
that are present in the RTF templates. Add `--strict-rtf` when you want missing
template atoms to fail as well. This is the main guard against accidentally
using the standard profile for a true-terminal chain, or the terminal profile
for an internalized heteroduplex.

The committed `examples/1efs_heteroduplex/1efs.pdb` is already in MOSAICS
naming and remains the reference template for hybrid residue/atom names.

## 5. Run a single-point energy and check it

Worked, runnable example using only committed files:

```bash
cd examples/1efs_heteroduplex
# replace with the path to your built MOSAICS serial binary
bash smoke_test_forcefield_panel.sh /path/to/mosaics_serial
```

This runs `OL21/OL3`, `OL15/OL3`, and `OL24/OL3` single points on the committed 1EFS
heteroduplex and compares the parsed MOSAICS initial energies against the
committed reference tables `expected_ol21_ol3.tsv`, `expected_ol15_ol3.tsv`,
and `expected_ol24_ol3.tsv`.
Exit code 0 means your build reproduces the validated energies.

Terminal-profile worked example using only committed files:

```bash
cd examples/pure_terminal_controls
# replace with the path to your built MOSAICS serial binary
bash smoke_test_terminal_profiles.sh /path/to/mosaics_serial
```

This runs `dna_acgt` and `rna_acgu` under `OL21/OL3`, `OL15/OL3`, and
`OL24/OL3` terminal profiles. Exit code 0 means the terminal profile can run
true 5-prime/3-prime canonical pure DNA and pure RNA controls and reproduce
the committed initial-energy references.

To run a single preset by hand:

```bash
cd examples/1efs_heteroduplex
/path/to/mosaics_serial mcmc_ol21_ol3_vac.input > my_run.log
python3 ../../tools/compare_mosaics_energy.py \
    --mosaics-log my_run.log \
    --reference-energy expected_ol21_ol3.tsv \
    --output-dir my_cmp --reference-label expected --fail-on-mismatch
```

The example input files run in vacuum with no distance-dependent dielectric and
no neutralization, matching the single-point validation conditions. For
production sampling you will change the simulation settings, but keep the
`\mol_parm_file` / `\*_database_file` block pointed at the deck and RTF for
your chosen force field and profile.
