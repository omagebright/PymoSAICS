# Architecture

## Core and GUI

`pymosaics.core` contains catalog, input generation, project preparation,
structure/topology validation, output discovery, and analysis code. Most of the
core uses only the Python standard library; structural landscape projection
uses NumPy supplied by normal PyMOL distributions.

The interface imports Qt only through `pymol.Qt`, so PyMOL chooses its bundled
Qt binding. PyMOL coordinates are read through `cmd.get_pdbstr`; structures and
trajectory states are returned through `cmd.load` and `cmd.frame`.

## Process boundary

MOSAICS starts through `QProcess.setProgram` and `setArguments`. No shell string
is constructed. Standard output and standard error are merged, displayed live,
and written to a timestamped log.

## Scientific assets

The plugin carries authorized compiled Apple-Silicon executables and validated
force-field/topology profiles. It carries no MOSAICS source. Runtime hashes and
provenance are recorded in `pymosaics/assets/ASSET_NOTICE.md`.

Each generated input references only its selected six-file force-field profile.
Short relative paths avoid legacy MOSAICS parser limits and make the exact run
deck portable and inspectable. Restaging never deletes unrelated user files.

## Scientific boundary

Presets are explicit starting points, not automatic scientific conclusions.
PymoSAICS checks file compatibility and reproduces configured calculations; it
does not establish equilibration, convergence, force-field suitability, a
global energy minimum, or a thermodynamic free-energy landscape.
