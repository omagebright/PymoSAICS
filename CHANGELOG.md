# Changelog

## 0.1.1 — 2026-07-21

- Corrected Windows test fixtures to use a native executable name.
- Normalized Windows temporary paths before portability assertions.
- Updated GitHub Actions to their Node.js 24-compatible releases.

## 0.1.0 — 2026-07-21

- Replaced the legacy Tk/Pmw interface with PyMOL's Qt compatibility layer.
- Added the modern PyMOL plugin entry point.
- Removed runtime downloading and 32-bit platform assumptions.
- Added an external-runtime setup and validation workflow.
- Added shell-free asynchronous execution, live logging, and cancellation.
- Added portable project and force-field placeholders.
- Added project input-reference validation and PyMOL output loading.
- Added standard per-user configuration paths.
- Added Windows, macOS, and Linux automated core tests.
