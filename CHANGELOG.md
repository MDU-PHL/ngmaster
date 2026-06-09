# Changelog

All notable changes to this project will be documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.7] - 2026-06-09

### Fixed

- Added preflight FASTA validation so empty or no-sequence inputs exit with a clear error before calling `mlst`. This removes noisy `any2fasta`/`CalledProcessError` output for empty FASTA files ([#48]).

## [2.0.6] - 2026-06-09

### Fixed
- Updated `match_porb` to load NG-STAR porB sequences into memory once and break after first match, improving performance by eliminating redundant file I/O and unnecessary iterations. Pre-compiled regex pattern for efficiency. Always parallelise scheme execution regardless of `--threads` setting to leverage concurrent subprocess calls.
- Parallelised scheme execution by default (`max_workers=2`) to ensure NG-MAST and NG-STAR run concurrently, improving overall runtime for typical use cases.
- Addressed slow runtime reported in [#50].

## [2.0.5] - 2026-05-26

### Added

- `$NGMASTER_DB` environment variable support: if set, `ngmaster` uses it as the database path without requiring `--db` on every call. `--db` on the command line takes precedence; if neither is set, the bundled database is used. Closes [#56].
- `--mincov` and `--minid` now validate that the supplied value is within 0–100 and exit with a clear error message if not.
- `--printseq` now prints an informational message to stderr explaining that only novel (`~n`) allele sequences are saved, and that no file is created when all alleles are exact matches. Closes [#42].

### Changed

- Default `--mincov` changed from 10 to 50, matching the `mlst` default.
- `--mincov` and `--minid` help text now states the default value and valid range (0–100).
- `--printseq` help text updated to clarify novel/partial-only behaviour.
- `--db` help text updated to document `$NGMASTER_DB` precedence rules.
- `--version` output format changed to the standard `ngmaster <version>` on the first line (was multi-line `ngmaster version:\n<version>`), making it machine-parseable. Closes [#49].

### Fixed

- Removed duplicate "Updating the allele databases" section from README (one used `ngmaster.py`, the other `ngmaster`); consolidated into a single correct section using `ngmaster`. Closes [#54].
- All remaining `ngmaster.py` references in source code and docs replaced with `ngmaster`.
- Fixed typo "unsucessful" → "unsuccessful" in test failure message.

## [2.0.4] - 2026-05-26

### Fixed

- Corrected DB update failure (`argument of type 'NoneType' is not iterable`) caused by a rauth 0.7.3 bug: `OAuth1Session.get()` was called without `params={}`, which triggered a `TypeError` in rauth's `_parse_optional_params`. All OAuth requests now pass `params={}` explicitly ([#60]).
- Replaced the deferred "test auth only once" pattern in `PubMLSTAuth` with upfront auth-mode selection: personal API key (stored via `mlstdb connect --db pubmlst --api-key`) is tried first, OAuth session tokens second, unauthenticated last.
- Fallback messaging now directs users to `mlstdb connect --db pubmlst --api-key` (recommended) or `mlstdb connect --db pubmlst` (OAuth) when no credentials are found.
- Added `User-Agent: ngmaster/<version>` header to all BIGSdb requests (API key and OAuth modes).
- Bumped `mlstdb` dependency from `==1.0.0` to `>=1.2.0` to pick up the upstream rauth fix and personal API key support.

## [2.0.3] - 2026-05-26

### Fixed

- Corrected `--minid` and `--mincov` argparse definitions: removed `nargs=1` and added `type=int` so both flags always produce a scalar integer. Previously, `nargs=1` caused argparse to return a list when the user supplied a value (e.g. `--mincov 20` produced `[20]`), which `str()` serialised as `"['20']"`, an argument that mlst rejects with a non-zero exit code ([#39]).
- Improved error reporting when an mlst subprocess fails: mlst's own stderr output is now relayed to the user via `msg()` before ngmaster exits, so the root cause is visible rather than just the command-line dump and exit code ([#57]).
- Added `TestMinidMincov` integration tests covering scalar value passing, combined flag usage, invalid type rejection, and mlst stderr surfacing.

## [2.0.2] - 2026-05-25

### Security

- Bumped `requests` minimum to `>=2.33.0` to address CVEs covering proxy header leakage, `verify=False` bypass propagation, `.netrc` credential exposure, and insecure temp file reuse ([#53]).
- Added explicit `urllib3>=2.7.0` floor to address eight CVEs covering redirect handling, decompression bombs, and cookie/header stripping ([#53]).
- Added explicit `certifi>=2023.7.22` floor to remove compromised root certificates from the trust store ([#53]).
- Added explicit `idna>=3.15` floor to address denial-of-service via specially crafted inputs to `idna.encode()` ([#53]).
- Bumped `numpy` minimum to `>=1.26.5` to address buffer overflow and null pointer CVEs ([#53]).
- Note: the `biopython` XXE advisory (via `Bio.Entrez`) has no upstream fix released yet; it does not affect ngmaster as `Bio.Entrez` is not used anywhere in the codebase ([#53]).

## [2.0.1] - 2026-05-25

### Fixed

- Removed duplicate fixture definitions in `tests/conftest.py` that caused pytest to fail at collection with `ValueError: duplicate fixture` before any test ran ([#55]).
- Pinned `python=${{ matrix.python-version }}` in the `mamba create` CI step so the conda environment uses the correct Python version from the matrix, preventing `pip install` from silently refusing to install due to a `python_requires` mismatch, the root cause of `ModuleNotFoundError: requests` in CI ([#55]).
- Added `-c conda-forge` channel to `mamba create` to resolve transitive dependencies (e.g. `rauth`) that bioconda alone cannot satisfy.
- Standardised all CI steps to `shell: bash -el {0}` so `conda activate` failures are caught immediately rather than silently continuing in the base environment.

### Changed

- CI Python matrix expanded to `3.10`, `3.11`, `3.12`, `3.13`, `3.14` with `fail-fast: false`.
- `python_requires` lowered to `>=3.10` in `setup.py` to align with the CI matrix.
- Removed legacy `test.yml` workflow (tested EOL Python 3.7–3.10, never ran `pytest`).
- Publish workflow (`publish.yml`) now triggers only on version tags (`v*.*.*`), not on every push to `master`.
- Renamed `publsh.yml` → `publish.yml` (fixed filename typo).

## [2.0.0] - 2026-03-24

See the [v2.0.0 release notes](https://github.com/MDU-PHL/ngmaster/releases/tag/v2.0.0).

[#39]: https://github.com/MDU-PHL/ngmaster/issues/39
[#42]: https://github.com/MDU-PHL/ngmaster/issues/42
[#48]: https://github.com/MDU-PHL/ngmaster/issues/48
[#49]: https://github.com/MDU-PHL/ngmaster/issues/49
[#50]: https://github.com/MDU-PHL/ngmaster/issues/50
[#53]: https://github.com/MDU-PHL/ngmaster/issues/53
[#54]: https://github.com/MDU-PHL/ngmaster/issues/54
[#55]: https://github.com/MDU-PHL/ngmaster/issues/55
[#56]: https://github.com/MDU-PHL/ngmaster/issues/56
[#57]: https://github.com/MDU-PHL/ngmaster/issues/57
[#60]: https://github.com/MDU-PHL/ngmaster/issues/60

[2.0.6]: https://github.com/MDU-PHL/ngmaster/releases/tag/v2.0.6
[2.0.7]: https://github.com/MDU-PHL/ngmaster/releases/tag/v2.0.7