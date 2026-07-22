# Bundled scientific assets

The PymoSAICS MIT licence covers the plugin code, not the bundled MOSAICS
executables or third-party force-field parameters.

The project maintainer confirmed authorization on 2026-07-21 to redistribute
the compiled executables and force-field/topology files in this package. No
MOSAICS source code is included.

## Executables

| Runtime | Platform | Provenance | SHA-256 |
|---|---|---|---|
| MOSAICS 3.9.1 | macOS arm64 | `version.3.9.1_bgq` | `a65d34474ba51c479352566928423c9560cb0776c1c37f17cae3bab59e9ab5ad` |
| MOSAICS experimental 2026-07-21 | macOS arm64 | validated stack: ff14SB/disulfides `94120f4`, CMAP `6341e5e`, PDB output `370a8b5` | `8a23058f859ff27409353b47c263280816f84c51c7d605e94897572ac1a43a15` |

Windows and Linux users must select a compatible custom executable until
native bundled builds are published.

## Force fields

The selectable panel contains parmbsc0, parmbsc1, OL15/OL3, OL21/OL3,
OL24/OL3, and ff14SB MOSAICS representations. Standard and true-terminal
nucleic-acid profiles are separate because their resolved torsion terms differ.
See `forcefields/panel/PROVENANCE.md` and `forcefields/panel/FORCEFIELD_USAGE.md`
for sources, validation boundaries, and citations.
