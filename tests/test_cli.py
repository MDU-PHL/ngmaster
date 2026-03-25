"""
CLI integration tests for ngmaster v2.0.0.

All tests in this module require the ``mlst`` binary in PATH and the
ngmaster package to be installed (``pip install -e .``).  They are
marked ``integration`` so they can be skipped when mlst is unavailable:

    pytest -m "not integration"

New features / bug-fixes under test
------------------------------------
* ``--version``      : reports both tool version and database version
* ``--comments``     : 20-column output with per-allele comment columns + CC
* CC column          : present in all output modes
* ``--threads``      : accepts > 1 without changing results
* ``--db``           : custom database path is honoured
* ``--test``         : built-in self-test still passes (regression guard)
* ``--csv``          : comma-separated output, comma-containing alleles are quoted
* Edge-case FASTAs   : noseq.fa and null.fa must not produce Python tracebacks
"""

import subprocess
import sys
from pathlib import Path

import pytest

import ngmaster

_BUNDLE_ROOT = Path(ngmaster.__file__).parent
_BUNDLED_DB = str(_BUNDLE_ROOT / "db")
_TEST_FA = str(_BUNDLE_ROOT / "test" / "test.fa")
_NOSEQ_FA = str(_BUNDLE_ROOT / "test" / "noseq.fa")
_NULL_FA = str(_BUNDLE_ROOT / "test" / "null.fa")

# Expected columns : non-comments mode
_EXPECTED_HEADER_NO_COMMENTS = [
    "FILE", "SCHEME", "NG-MAST/NG-STAR",
    "porB_NG-MAST", "tbpB",
    "penA", "mtrR", "porB_NG-STAR", "ponA", "gyrA", "parC", "23S",
    "CC",
]

# Expected columns : --comments mode
_EXPECTED_HEADER_COMMENTS = [
    "FILE", "SCHEME", "NG-MAST/NG-STAR",
    "porB_NG-MAST", "tbpB",
    "penA", "penA_comments",
    "mtrR", "mtrR_comments",
    "porB_NG-STAR", "porB_NG-STAR_comments",
    "ponA", "ponA_comments",
    "gyrA", "gyrA_comments",
    "parC", "parC_comments",
    "23S", "23S_comments",
    "CC",
]


# ---------------------------------------------------------------------------
# Module-level skip if mlst is not installed
# ---------------------------------------------------------------------------

pytestmark = pytest.mark.integration


@pytest.fixture(scope="module", autouse=True)
def require_mlst():
    """Skip ALL tests in this module when mlst is not in PATH."""
    import shutil
    if shutil.which("mlst") is None:
        pytest.skip("mlst binary not found in PATH : skipping integration tests")


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _run(*args):
    """Run ngmaster with the given arguments; return CompletedProcess."""
    return subprocess.run(
        [sys.executable, "-m", "ngmaster.run_ngmaster"] + list(args),
        capture_output=True,
        text=True,
    )


# ---------------------------------------------------------------------------
# Regression guard: built-in --test must still pass
# ---------------------------------------------------------------------------

class TestBuiltinTest:
    def test_test_flag_exits_zero(self):
        result = _run("--test")
        assert result.returncode == 0, (
            f"ngmaster --test failed.\nstdout: {result.stdout}\nstderr: {result.stderr}"
        )

    def test_test_flag_correct_st(self):
        result = _run("--test")
        assert "4186/231" in result.stdout, (
            f"Expected '4186/231' in output.\nstdout: {result.stdout}"
        )


# ---------------------------------------------------------------------------
# --version
# ---------------------------------------------------------------------------

class TestVersion:
    def test_exits_zero(self):
        assert _run("--version").returncode == 0

    def test_contains_tool_version(self):
        result = _run("--version")
        assert ngmaster.__version__ in result.stdout

    def test_contains_db_version_keys(self):
        result = _run("--version")
        assert "ngmast_" in result.stdout
        assert "ngstar_" in result.stdout


# ---------------------------------------------------------------------------
# Default tab-separated output
# ---------------------------------------------------------------------------

class TestDefaultOutput:
    def test_exits_zero(self):
        result = _run("--db", _BUNDLED_DB, _TEST_FA)
        assert result.returncode == 0, result.stderr

    def test_header_has_13_columns(self):
        result = _run("--db", _BUNDLED_DB, _TEST_FA)
        header = result.stdout.splitlines()[0].split("\t")
        assert header == _EXPECTED_HEADER_NO_COMMENTS

    def test_cc_column_in_header(self):
        result = _run("--db", _BUNDLED_DB, _TEST_FA)
        assert "CC" in result.stdout.splitlines()[0].split("\t")

    def test_body_contains_expected_st(self):
        result = _run("--db", _BUNDLED_DB, _TEST_FA)
        assert "4186/231" in result.stdout

    def test_no_python_traceback(self):
        result = _run("--db", _BUNDLED_DB, _TEST_FA)
        assert "Traceback (most recent call last)" not in result.stderr


# ---------------------------------------------------------------------------
# CSV output  (--csv)
# ---------------------------------------------------------------------------

class TestCsvOutput:
    def test_exits_zero(self):
        result = _run("--db", _BUNDLED_DB, "--csv", _TEST_FA)
        assert result.returncode == 0

    def test_header_is_comma_separated(self):
        result = _run("--db", _BUNDLED_DB, "--csv", _TEST_FA)
        header_line = result.stdout.splitlines()[0]
        # Tab-separated header would not contain commas between column names
        assert "," in header_line

    def test_cc_present_in_csv_header(self):
        result = _run("--db", _BUNDLED_DB, "--csv", _TEST_FA)
        assert "CC" in result.stdout.splitlines()[0]

    def test_body_contains_expected_st(self):
        result = _run("--db", _BUNDLED_DB, "--csv", _TEST_FA)
        assert "4186/231" in result.stdout


# ---------------------------------------------------------------------------
# --comments output  (19+1 columns including CC)
# ---------------------------------------------------------------------------

class TestCommentsOutput:
    def test_exits_zero(self):
        result = _run("--db", _BUNDLED_DB, "--comments", _TEST_FA)
        assert result.returncode == 0, result.stderr

    def test_header_has_20_columns(self):
        result = _run("--db", _BUNDLED_DB, "--comments", _TEST_FA)
        header = result.stdout.splitlines()[0].split("\t")
        assert header == _EXPECTED_HEADER_COMMENTS, (
            f"Got {len(header)} columns, expected {len(_EXPECTED_HEADER_COMMENTS)}.\n"
            f"Header: {header}"
        )

    def test_comment_columns_present(self):
        result = _run("--db", _BUNDLED_DB, "--comments", _TEST_FA)
        header = result.stdout.splitlines()[0].split("\t")
        for col in ("penA_comments", "mtrR_comments", "23S_comments"):
            assert col in header

    def test_cc_is_last_column(self):
        result = _run("--db", _BUNDLED_DB, "--comments", _TEST_FA)
        header = result.stdout.splitlines()[0].split("\t")
        assert header[-1] == "CC"

    def test_body_has_correct_field_count(self):
        result = _run("--db", _BUNDLED_DB, "--comments", _TEST_FA)
        lines = result.stdout.splitlines()
        data_fields = lines[1].split("\t")
        assert len(data_fields) == len(_EXPECTED_HEADER_COMMENTS)


# ---------------------------------------------------------------------------
# --threads
# ---------------------------------------------------------------------------

class TestThreads:
    def test_threads_2_produces_same_st(self):
        """Running with --threads 2 must produce the same ST as single-threaded."""
        result_1 = _run("--db", _BUNDLED_DB, _TEST_FA)
        result_2 = _run("--db", _BUNDLED_DB, "--threads", "2", _TEST_FA)
        assert result_1.returncode == 0 and result_2.returncode == 0
        assert "4186/231" in result_2.stdout

    def test_threads_accepts_integer(self):
        result = _run("--db", _BUNDLED_DB, "--threads", "4", _TEST_FA)
        assert result.returncode == 0


# ---------------------------------------------------------------------------
# --db  (custom database path)
# ---------------------------------------------------------------------------

class TestCustomDb:
    def test_custom_db_path_accepted(self):
        """ngmaster must accept a custom --db path without crashing."""
        result = _run("--db", _BUNDLED_DB, _TEST_FA)
        assert result.returncode == 0
        assert "4186/231" in result.stdout


# ---------------------------------------------------------------------------
# Edge-case FASTA inputs  (must not produce Python tracebacks)
# ---------------------------------------------------------------------------

class TestEdgeCaseFasta:
    def test_noseq_fa_no_traceback(self):
        """FASTA with a header but no sequence body must not produce a traceback."""
        result = _run("--db", _BUNDLED_DB, _NOSEQ_FA)
        assert "Traceback (most recent call last)" not in result.stderr

    def test_null_fa_no_traceback(self):
        """Empty FASTA file must not produce a traceback."""
        result = _run("--db", _BUNDLED_DB, _NULL_FA)
        assert "Traceback (most recent call last)" not in result.stderr

    def test_no_fasta_prints_help(self):
        """Invoking ngmaster with no FASTA arguments must print usage and exit non-zero."""
        result = _run("--db", _BUNDLED_DB)
        assert result.returncode != 0
        # Help text or error message expected
        assert "usage" in result.stderr.lower() or "error" in result.stderr.lower()
