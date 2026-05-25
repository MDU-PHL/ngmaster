# Changelog

All notable changes to this project will be documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.1] - 2026-05-25

### Fixed

- Removed duplicate fixture definitions in `tests/conftest.py` that caused pytest to
  fail at collection with `ValueError: duplicate fixture` before any test ran (#55).
- Pinned `python=${{ matrix.python-version }}` in the `mamba create` CI step so the
  conda environment uses the correct Python version from the matrix, preventing `pip install`
  from silently refusing to install due to a `python_requires` mismatch, the root cause of
  `ModuleNotFoundError: requests` in CI (#55).
- Added `-c conda-forge` channel to `mamba create` to resolve transitive dependencies
  (e.g. `rauth`) that bioconda alone cannot satisfy.
- Standardised all CI steps to `shell: bash -el {0}` so `conda activate` failures are
  caught immediately rather than silently continuing in the base environment.

### Changed

- CI Python matrix expanded to `3.10`, `3.11`, `3.12`, `3.13`, `3.14` with `fail-fast: false`.
- `python_requires` lowered to `>=3.10` in `setup.py` to align with the CI matrix.
- Removed legacy `test.yml` workflow (tested EOL Python 3.7–3.10, never ran `pytest`).
- Publish workflow (`publish.yml`) now triggers only on version tags (`v*.*.*`), not on
  every push to `master`.
- Renamed `publsh.yml` → `publish.yml` (fixed filename typo).

## [2.0.0] - 2026-03-24

See the [v2.0.0 release notes](https://github.com/MDU-PHL/ngmaster/releases/tag/v2.0.0).
