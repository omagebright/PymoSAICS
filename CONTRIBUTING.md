# Contributing

Keep the public plugin transparent and platform-independent.

Before submitting a change:

```bash
python -m compileall -q pymosaics
python -m unittest discover -v
python scripts/build_release.py
```

Do not add:

- a shell command assembled from user input;
- a platform-specific dependency outside PyMOL;
- a MOSAICS binary or force-field data without explicit redistribution rights;
- a scientific input template without provenance and regression evidence;
- a write to the PyMOL or MOSAICS installation directory.

Changes to parameter parsing should include fixtures for paths containing spaces
and for all supported operating-system path conventions.
