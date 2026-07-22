# Public PymoSAICS redesign

## Goal

Ship one installable PyMOL plugin package for Windows, macOS, and Linux without
redistributing MOSAICS itself. The plugin must make execution transparent: the
user supplies the MOSAICS executable, force-field directory, project directory,
and parameter input file.

## Settled decisions

- Use PyMOL's supported Qt compatibility layer (`pymol.Qt`), not Tk/Pmw.
- Use the modern `__init_plugin__` entry point and `addmenuitemqt`.
- Do not bundle or download MOSAICS. Its upstream repository has no declared
  licence or platform releases as of 2026-07-21.
- Do not silently generate scientific force-field inputs. Version 0.1 runs a
  manually reviewable MOSAICS parameter file.
- Never invoke a shell. Pass an argument list directly to `QProcess` so spaces
  and shell metacharacters in user paths are safe.
- Store configuration in the operating system's user configuration directory,
  never in the PyMOL installation or plugin directory.
- Support `${PYMOSAICS_FORCEFIELD_DIR}` and `${PROJECT_DIR}` in input files.
  Resolve these into a generated copy and leave the source file unchanged.
- Keep all runtime, configuration, and validation logic independent of PyMOL
  and Qt so it can be tested on every platform.

## User workflow

1. Install `PymoSAICS-0.1.0.zip` through PyMOL's Plugin Manager.
2. Open **Plugin > PymoSAICS**.
3. Select the locally installed MOSAICS executable and force-field directory.
4. Select a MOSAICS parameter input file.
5. Review validation results and the exact command preview.
6. Run or stop MOSAICS while viewing the live log.
7. Load a generated PDB trajectory into PyMOL.

## Architecture

- `pymosaics/core/config.py`: portable, atomic JSON configuration.
- `pymosaics/core/runtime.py`: runtime validation and command construction.
- `pymosaics/core/project.py`: placeholder expansion, input reference checks,
  prepared-run paths, and output discovery.
- `pymosaics/gui.py`: Qt interface and asynchronous `QProcess` lifecycle.
- `pymosaics/plugin.py`: one-window lifecycle for PyMOL.
- `pymosaics/__init__.py`: PyMOL plugin registration.

## Public support boundary

The plugin UI is designed for current Qt-based PyMOL on Windows 10/11 x64,
macOS Intel/Apple Silicon, and Ubuntu 22.04/24.04 x64. Actual simulation support
also requires a compatible MOSAICS executable and force-field installation for
the user's platform. The plugin cannot make an incompatible MOSAICS binary
portable.

## Validation

- Unit tests for Windows, macOS, and Linux configuration paths.
- Unit tests for executable and force-field validation.
- Unit tests for safe argument construction and placeholder resolution.
- Unit tests proving source inputs are not modified.
- Unit tests for output discovery and plugin registration.
- GitHub Actions matrix across all three operating systems.
- Local smoke test with the installed PyMOL application.

## Deferred work

- Scientifically validated input-file builders for specific force-field profiles.
- Automatic runtime installation, pending a MOSAICS licence and official builds.
- Remote/HPC job submission.
- Project deletion or other destructive file management.
