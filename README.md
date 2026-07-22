# PymoSAICS

PymoSAICS is a PyMOL plugin for validating, running, stopping, and visualizing
simulations performed by a separately installed MOSAICS executable.

The public package is deliberately small and transparent:

- one Plugin Manager ZIP works on Windows, macOS, and Linux;
- the interface uses PyMOL's supported Qt layer, with no Tk, Pmw, or XQuartz;
- commands are executed without a shell;
- paths containing spaces are supported;
- parameter files remain visible and manually reviewable;
- the source parameter file is never modified;
- logs and resolved portable inputs are retained with the project.

## Important scope

PymoSAICS does **not** contain or download MOSAICS, force-field data, or
scientific input templates. Users need a MOSAICS executable compiled for their
operating system and a compatible force-field directory.

The [upstream MOSAICS repository](https://github.com/pminary/MOSAICS) did not
declare a licence or publish platform releases when this package was prepared.
Keeping the runtime external avoids making unsupported redistribution claims.

## Supported interface

| Platform | Public plugin target | External requirement |
|---|---|---|
| Windows | Windows 10/11 x64 | Compatible `mosaics.exe` |
| macOS | Intel and Apple Silicon | Compatible native or translated executable |
| Linux | Ubuntu 22.04/24.04 x64 | Compatible executable |

The interface targets Qt-based PyMOL 2.6 and later, including PyMOL 3.x. The
included automated core test matrix is configured for Python 3.9 through 3.12
on all three operating systems. Actual PyMOL distributions bundle their own
Python and Qt versions.

## Install

1. In PyMOL, open **Plugin > Plugin Manager > Install New Plugin**.
2. Paste this URL into **Install from PyMOLWiki or any URL**:

   ```text
   https://github.com/omagebright/PymoSAICS/releases/download/v0.1.1/PymoSAICS-0.1.1.zip
   ```

3. Select **Install**, restart PyMOL, then open **Plugin > PymoSAICS**.

If the old `PymoSAICS` plugin is installed, disable its **Load on startup** entry
and enable the new lowercase `pymosaics` entry. No extraction, administrator
privileges, XQuartz installation, or separate Python packages are required.

## First-time setup

In the **Setup** tab, select:

1. your external MOSAICS executable;
2. the root of your force-field data (normally containing `top_database` or
   `pot_database`);
3. an optional default project directory.

Choose **Validate and save**. The configuration is written to the normal user
configuration location:

- Windows: `%APPDATA%\PymoSAICS\config.json`
- macOS: `~/Library/Application Support/PymoSAICS/config.json`
- Linux: `${XDG_CONFIG_HOME:-~/.config}/pymosaics/config.json`

## Run a simulation

1. Open the **Run** tab.
2. Select a manually reviewed MOSAICS parameter input file.
3. Select **Validate** and read every diagnostic.
4. Confirm the displayed program, argument, and working directory.
5. Select **Run MOSAICS**.
6. Watch the live output or use **Stop** to terminate the process.
7. Load a detected output PDB into PyMOL.

PymoSAICS invokes the equivalent of:

```text
PROGRAM: /path/to/mosaics
ARGUMENT 1: /path/to/project/parameters.input
WORKING DIRECTORY: /path/to/project
```

These values are passed as separate process arguments. They are not assembled
into a shell command.

## Portable parameter files

Two placeholders can remove machine-specific absolute paths:

```text
\mol_parm_file{${PYMOSAICS_FORCEFIELD_DIR}/top_database/amber/profile.rtf}
\pos_init_file{${PROJECT_DIR}/start.pdb}
```

Immediately before execution, PymoSAICS creates a resolved copy under:

```text
PROJECT/.pymosaics/resolved/
```

The original parameter file is unchanged. Run logs are stored under:

```text
PROJECT/.pymosaics/logs/
```

## Build and test from source

No third-party development dependencies are needed for the core tests:

```bash
python -m unittest discover -v
python scripts/build_release.py
python -m zipfile -l dist/PymoSAICS-0.1.1.zip
```

The build command creates the Plugin Manager archive in `dist/`.

## Scientific responsibility

Passing PymoSAICS validation means paths and runtime configuration are
structurally usable. It does not establish that a force field, input structure,
simulation method, or result is scientifically valid. Review those choices with
the relevant MOSAICS documentation and project maintainers.

## Licence

PymoSAICS plugin code is released under the MIT Licence. MOSAICS and force-field
data are separate works and are not covered or distributed by this licence.

## Support

Report reproducible plugin problems through the
[PymoSAICS issue tracker](https://github.com/omagebright/PymoSAICS/issues).
Questions about MOSAICS itself, its scientific methods, or force-field data
belong with the corresponding MOSAICS maintainers.
