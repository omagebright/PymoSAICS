# Architecture

## Dependency boundary

`pymosaics.core` uses only the Python standard library. It contains all file,
path, configuration, validation, and command-construction logic. Importing it
does not require PyMOL or Qt.

The GUI imports Qt only through `pymol.Qt`. This lets PyMOL select its bundled
Qt binding and avoids asking users to install PyQt or PySide separately.

## Process boundary

The GUI starts MOSAICS with `QProcess.setProgram` and `setArguments`. No shell,
redirection operator, Perl interpreter, `chmod`, `rm`, `zip`, or `unzip` command
is used. Standard output and standard error are merged, shown live, and written
to a per-run log.

## Filesystem boundary

The plugin writes only to:

- the operating system's per-user configuration directory; and
- `.pymosaics` inside the selected project directory.

It does not write into the PyMOL application, its plugin installation directory,
the MOSAICS installation, or the force-field directory.

## Scientific boundary

The legacy plugin embedded a dated AMBER 99-bs0 parameter generator and made
scientific choices while building input. This redesign excludes that generator
from the public core. Reintroducing a profile requires a versioned template,
provenance, regression fixtures, and comparison against an accepted manual run.
