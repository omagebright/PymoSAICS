# User guide

## What PymoSAICS controls

PymoSAICS controls only four explicit values:

1. the external MOSAICS program;
2. the parameter input passed as its single argument;
3. the project directory used as the process working directory;
4. the output log captured from the process.

It does not alter the selected parameter input, select a force field, or infer a
scientific protocol.

## Runtime validation

The Setup tab checks that:

- the executable exists;
- it has an executable extension on Windows or executable permission on
  macOS/Linux;
- the force-field directory exists;
- the usual `top_database` or `pot_database` layout is present, when applicable.

A nonstandard force-field layout produces a warning rather than an error.

## Project validation

The Run tab reads the parameter input as UTF-8 and checks recognized input
references such as `mol_parm_file`, `bond_database_file`,
`region_database_file`, and `pos_init_file`. Missing referenced inputs block a
run. Output paths are not required to exist before a run.

Unrecognized MOSAICS options remain untouched. This is intentional: PymoSAICS
does not reinterpret the full MOSAICS input language.

## Portable paths

Use `${PYMOSAICS_FORCEFIELD_DIR}` for the force-field root configured in Setup,
and `${PROJECT_DIR}` for the directory containing the selected parameter input.
PymoSAICS replaces both with absolute paths using forward slashes and writes a
resolved copy. Any other unresolved `${NAME}` placeholder is a validation error.

## Logs and reproducibility

Every run has a UTC-stamped log under `.pymosaics/logs`. If placeholders were
used, the exact resolved input is retained under `.pymosaics/resolved`. Together
with the source input, these files show precisely what was executed.

## Stopping a run

Stop first requests normal termination. If the program does not exit within
three seconds, the plugin kills it. MOSAICS may leave partial output files after
either action; PymoSAICS never reports a stopped process as successful.

## Loading output

You may explicitly select an output PDB. Otherwise, PymoSAICS searches the
project directory for `simulation.pdb`, `simulation_result.pdb`, and
`*.pos_out.pdb`, then loads the newest match. It does not rename or rewrite the
file.
