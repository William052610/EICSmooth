# Changelog

## 2026-08-10 — Installation and dependency checks

- Unified runtime dependency declarations in `pyproject.toml` and made
  `requirements.txt` install the local EICSmooth project as well as its
  dependencies.
- Pinned the tested runtime set: NumPy 1.26.4, Matplotlib 3.9.3, and
  Pyteomics 5.0.1 with XML support.
- Made a project-local virtual environment the default Windows installation
  path so packages already present in Anaconda `base` cannot cause ABI or
  unrelated dependency conflicts.
- Added an explicit dependency consistency check to the documented
  verification workflow and CI.

## 2026-08-09 — GUI and environment updates

- Simplified EICSmooth to a single GUI program.
- Removed the optional batch runner, batch manifest, and related command-line interface.
- Restored Python 3.10-3.12 support and pinned NumPy 1.26.4 to prevent binary
  incompatibility with scientific packages compiled against NumPy 1.x.
- Converted remaining GUI, validation, status, and error messages to English.
- Improved dependency-import errors so the original cause and active Python
  interpreter are shown instead of reporting every failure as a missing package.
- Removed local virtual-environment and IDE-generated files from the release archive.

## 2026-08-09 — Windows installation updates

- Raised the NumPy lower bound to 2.1 for Python 3.13 wheel compatibility.
- Updated the Matplotlib and Pyteomics lower bounds to tested modern releases.
- Added a Windows one-click installer that uses an isolated `.venv` and binary wheels.
- Added a Python 3.10-3.13 compatibility check and clearer installation guidance.

## 2026-08-09 — Initial release

- Released the EICSmooth GUI and optional batch workflow.
- Reproduced five MS-DIAL smoothing methods and two MZmine smoothing methods.
- Added fixed raw-EIC-derived peak boundaries.
- Added baseline-corrected height, signed area, FWHM, apex RT, and ideal-slope metrics.
- Added raw-relative log2 retention/ratio values, apex shift, and ideal-slope change.
- Added synthetic example data, automated tests, citation metadata, and CI configuration.
