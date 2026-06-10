"""
CLI integration tests for ngmaster v2.0.0.

Modelled after the mlst test/test.sh bats suite:
  – exit code is checked FIRST in every test
  – bad-option / /dev/null / no-args paths are covered before the happy path
  – no Python tracebacks must escape to the user under any input

Test categories
---------------
1.  Syntax         : package importability (equiv. of ``perl -c``)
2.  Help / version : always exit 0; output content checks
3.  Error inputs   : always exit non-zero; no Python traceback
4.  Happy path     : exit 0, correct ST (4186/231), header + data rows
5.  Output formats : tab (13 cols), CSV (comma-sep), --comments (20 cols)
6.  v2 features    : --threads, --db custom path, --printseq
7.  Regression     : --test built-in flag (issue #37/#38)

All tests require ``mlst`` in PATH (provided by the pixi env) and the
ngmaster package installed (``pip install -e .``).  They are marked
``integration`` so they can be skipped without mlst:

    pytest -m "not integration"
    pytest tests/test_cli.py -v          # full CLI suite
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
# Helpers
# ---------------------------------------------------------------------------

def _run(*args):
    """Run ngmaster with the given arguments; return CompletedProcess."""
    return subprocess.run(
        [sys.executable, "-m", "ngmaster.run_ngmaster"] + list(args),
        capture_output=True,
        text=True,
    )


def _no_traceback(result):
    """Assert no Python traceback appeared in stderr."""
    assert "Traceback (most recent call last)" not in result.stderr, (
        f"Python traceback in stderr:\n{result.stderr}"
    )


# ---------------------------------------------------------------------------
# 1. Syntax / importability  (equiv. of mlst's `perl -c`)
# ---------------------------------------------------------------------------

class TestSyntax:
    def test_package_imports_cleanly(self):
        """All three ngmaster modules must import without errors."""
        result = subprocess.run(
            [sys.executable, "-c",
             "import ngmaster; import ngmaster.utils; import ngmaster.run_ngmaster"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, f"Import error:\n{result.stderr}"


# ---------------------------------------------------------------------------
# 2. --help and --version  (must always exit 0)
# ---------------------------------------------------------------------------

class TestHelpAndVersion:
    def test_help_exits_zero(self):
        assert _run("--help").returncode == 0

    def test_help_mentions_key_options(self):
        out = _run("--help").stdout
        for flag in ("--db", "--csv", "--comments", "--threads", "--updatedb"):
            assert flag in out, f"'{flag}' missing from --help output"

    def test_help_mentions_threads(self):
        """--threads is a v2.0.0 addition; must appear in --help (mirrors mlst check)."""
        assert "--threads" in _run("--help").stdout

    def test_version_exits_zero(self):
        assert _run("--version").returncode == 0

    def test_version_contains_tool_version(self):
        assert ngmaster.__version__ in _run("--version").stdout

    def test_version_contains_db_version_keys(self):
        out = _run("--version").stdout
        assert "ngmast_" in out
        assert "ngstar_" in out


# ---------------------------------------------------------------------------
# 3. Error inputs  (all must exit non-zero, no Python traceback)
#    modelled after mlst's bad-option, null-input and /dev/null tests
# ---------------------------------------------------------------------------

class TestErrorInputs:
    def test_no_arguments_exits_nonzero(self):
        """Bare `ngmaster` with no arguments must exit non-zero."""
        result = _run()
        assert result.returncode != 0
        _no_traceback(result)

    def test_bad_option_exits_nonzero(self):
        """Unknown option must exit non-zero (mlst: `run ! $EXE --doesnotexist`)."""
        result = _run("--doesnotexist")
        assert result.returncode != 0
        _no_traceback(result)

    def test_dev_null_exits_nonzero(self):
        """Empty /dev/null input must exit non-zero (mlst: `run ! $EXE /dev/null`)."""
        result = _run("--db", _BUNDLED_DB, "/dev/null")
        assert result.returncode != 0
        _no_traceback(result)

    def test_nonexistent_fasta_exits_nonzero(self):
        result = _run("--db", _BUNDLED_DB, "/tmp/does_not_exist_xyz.fa")
        assert result.returncode != 0
        _no_traceback(result)

    def test_nonexistent_db_exits_nonzero(self):
        result = _run("--db", "/nonexistent/path", _TEST_FA)
        assert result.returncode != 0
        _no_traceback(result)

    def test_null_fa_no_traceback(self):
        """Completely empty FASTA file must not produce a Python traceback."""
        _no_traceback(_run("--db", _BUNDLED_DB, _NULL_FA))

    def test_noseq_fa_no_traceback(self):
        """FASTA with a header but no sequence must not produce a Python traceback."""
        _no_traceback(_run("--db", _BUNDLED_DB, _NOSEQ_FA))

    def test_empty_fasta_reports_clean_error(self, tmp_path):
        fasta = tmp_path / "empty.fasta"
        fasta.touch()

        result = subprocess.run(
            [sys.executable, "-m", "ngmaster.run_ngmaster", str(fasta)],
            capture_output=True,
            text=True,
        )

        assert result.returncode != 0
        assert "ERROR: Input FASTA is empty or contains no readable sequence" in result.stderr
        assert str(fasta) in result.stderr
        assert "CalledProcessError" not in result.stderr
        assert "Command '[" not in result.stderr
        assert "COMMAND FAILED" not in result.stderr
        assert "Traceback" not in result.stderr


# ---------------------------------------------------------------------------
# 4. Happy path  (exit 0, correct ST, correct structure)
# ---------------------------------------------------------------------------

class TestHappyPath:
    def test_exits_zero(self):
        result = _run("--db", _BUNDLED_DB, _TEST_FA)
        assert result.returncode == 0, f"stderr: {result.stderr}"

    def test_correct_st(self):
        """test.fa must type as NG-MAST 4186 / NG-STAR 231."""
        assert "4186/231" in _run("--db", _BUNDLED_DB, _TEST_FA).stdout

    def test_output_has_header_and_data_row(self):
        lines = _run("--db", _BUNDLED_DB, _TEST_FA).stdout.strip().splitlines()
        assert len(lines) >= 2

    def test_no_traceback_on_good_input(self):
        _no_traceback(_run("--db", _BUNDLED_DB, _TEST_FA))

    def test_multiple_fasta_files_produce_two_rows(self, tmp_path):
        """Two distinct FASTA paths must produce two data rows (mlst multi-file test)."""
        import shutil
        fasta2 = str(tmp_path / "test_copy.fa")
        shutil.copy(_TEST_FA, fasta2)
        result = _run("--db", _BUNDLED_DB, _TEST_FA, fasta2)
        assert result.returncode == 0, f"stderr: {result.stderr}"
        data_lines = result.stdout.strip().splitlines()[1:]  # skip header
        assert len(data_lines) == 2


# ---------------------------------------------------------------------------
# 5. Output formats
# ---------------------------------------------------------------------------

class TestTabOutput:
    def test_header_columns_exact(self):
        header = _run("--db", _BUNDLED_DB, _TEST_FA).stdout.splitlines()[0].split("\t")
        assert header == _EXPECTED_HEADER_NO_COMMENTS

    def test_cc_is_last_column(self):
        header = _run("--db", _BUNDLED_DB, _TEST_FA).stdout.splitlines()[0].split("\t")
        assert header[-1] == "CC"

    def test_data_row_field_count_matches_header(self):
        lines = _run("--db", _BUNDLED_DB, _TEST_FA).stdout.strip().splitlines()
        assert len(lines[1].split("\t")) == len(lines[0].split("\t"))


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
# 6. New v2.0.0 features
# ---------------------------------------------------------------------------

class TestThreads:
    def test_threads_2_exits_zero(self):
        assert _run("--db", _BUNDLED_DB, "--threads", "2", _TEST_FA).returncode == 0

    def test_threads_2_correct_st(self):
        """--threads 2 must produce the same ST as single-threaded."""
        assert "4186/231" in _run("--db", _BUNDLED_DB, "--threads", "2", _TEST_FA).stdout

    def test_threads_invalid_value_exits_nonzero(self):
        """Non-integer --threads value must exit non-zero gracefully."""
        result = _run("--db", _BUNDLED_DB, "--threads", "notanumber", _TEST_FA)
        assert result.returncode != 0
        _no_traceback(result)


class TestCustomDb:
    def test_explicit_bundled_db_path(self):
        """--db pointing to the bundled DB explicitly must produce the correct ST."""
        result = _run("--db", _BUNDLED_DB, _TEST_FA)
        assert result.returncode == 0
        assert "4186/231" in result.stdout

    def test_versioned_test_db_snapshot(self):
        """Versioned DB snapshot (test/test_db/db_v20250416) must work if present."""
        test_db = str(
            Path(__file__).parent.parent / "test" / "test_db" / "db_v20250416"
        )
        if not Path(test_db).exists():
            pytest.skip("test/test_db/db_v20250416 not present in working tree")
        result = _run("--db", test_db, _TEST_FA)
        assert result.returncode == 0


class TestPrintseq:
    def test_printseq_exits_zero(self, tmp_path):
        """--printseq must exit 0 even when all alleles are exact matches.

        mlst's --novel flag (used internally) only writes sequences for novel
        (~n) alleles. When all alleles are exact matches, no output file is
        created -- this is correct behaviour, not a bug (issue #42).
        """
        result = subprocess.run(
            [sys.executable, "-m", "ngmaster.run_ngmaster",
             "--db", _BUNDLED_DB, "--printseq", "alleles.fa", _TEST_FA],
            capture_output=True, text=True, cwd=str(tmp_path),
        )
        assert result.returncode == 0, f"stderr: {result.stderr}"
        _no_traceback(result)

    def test_printseq_note_in_stderr(self, tmp_path):
        """--printseq must emit an explanatory NOTE to stderr (issue #42)."""
        result = subprocess.run(
            [sys.executable, "-m", "ngmaster.run_ngmaster",
             "--db", _BUNDLED_DB, "--printseq", "alleles.fa", _TEST_FA],
            capture_output=True, text=True, cwd=str(tmp_path),
        )
        assert "NOTE" in result.stderr, (
            f"Expected NOTE message in stderr.\nstderr: {result.stderr}"
        )
        assert "novel" in result.stderr.lower()

    def test_printseq_no_files_for_exact_alleles(self, tmp_path):
        """No output files created when all alleles are exact matches (issue #42)."""
        subprocess.run(
            [sys.executable, "-m", "ngmaster.run_ngmaster",
             "--db", _BUNDLED_DB, "--printseq", "alleles.fa", _TEST_FA],
            capture_output=True, text=True, cwd=str(tmp_path),
        )
        created = list(tmp_path.iterdir())
        assert len(created) == 0, (
            f"Expected no files for exact-match alleles, but found: {[f.name for f in created]}"
        )


# ---------------------------------------------------------------------------
# 7. Regression: built-in --test flag  (issue #37 / #38)
# ---------------------------------------------------------------------------

class TestRegression:
    def test_builtin_test_flag_exits_zero(self):
        result = _run("--test")
        assert result.returncode == 0, (
            f"ngmaster --test failed.\nstdout: {result.stdout}\nstderr: {result.stderr}"
        )

    def test_builtin_test_flag_correct_st(self):
        assert "4186/231" in _run("--test").stdout

    def test_builtin_test_no_traceback(self):
        _no_traceback(_run("--test"))


# ---------------------------------------------------------------------------
# 8. Regression: --minid / --mincov flag passing  (issue #39 / #57)
# ---------------------------------------------------------------------------

class TestMinidMincov:
    """Regression tests for issue #39 (nargs=1 caused list-as-string to be passed to mlst)
    and issue #57 (mlst stderr not surfaced on failure).
    """

    def test_minid_flag_exits_zero(self):
        """--minid must be accepted as a scalar value and passed correctly to mlst."""
        result = _run("--minid", "85", _TEST_FA)
        assert result.returncode == 0, (
            f"--minid 85 failed.\nstdout: {result.stdout}\nstderr: {result.stderr}"
        )
        _no_traceback(result)

    def test_mincov_flag_exits_zero(self):
        """--mincov must be accepted as a scalar value and passed correctly to mlst."""
        result = _run("--mincov", "15", _TEST_FA)
        assert result.returncode == 0, (
            f"--mincov 15 failed.\nstdout: {result.stdout}\nstderr: {result.stderr}"
        )
        _no_traceback(result)

    def test_minid_mincov_combined_exits_zero(self):
        """Combining --minid and --mincov must produce the correct ST."""
        result = _run("--minid", "85", "--mincov", "15", _TEST_FA)
        assert result.returncode == 0, (
            f"--minid 85 --mincov 15 failed.\nstdout: {result.stdout}\nstderr: {result.stderr}"
        )
        _no_traceback(result)
        assert "4186/231" in result.stdout

    def test_minid_invalid_type_exits_nonzero(self):
        """Non-integer --minid must exit non-zero with no Python traceback (issue #57)."""
        result = _run("--minid", "notanumber", _TEST_FA)
        assert result.returncode != 0
        _no_traceback(result)

    def test_mincov_invalid_type_exits_nonzero(self):
        """Non-integer --mincov must exit non-zero with no Python traceback (issue #57)."""
        result = _run("--mincov", "notanumber", _TEST_FA)
        assert result.returncode != 0
        _no_traceback(result)

    def test_mlst_stderr_surfaced_on_failure(self):
        """When mlst fails, its stderr must appear in ngmaster stderr (issue #57).

        null.fa is an empty file that causes mlst to exit non-zero with a
        descriptive error message. ngmaster must relay that message to stderr
        rather than only showing the CalledProcessError exception summary.
        """
        result = _run(_NULL_FA)
        assert result.returncode != 0
        _no_traceback(result)
        # mlst writes its own error to stderr; ngmaster must relay it.
        assert result.stderr.strip() != "", (
            "Expected mlst error message in stderr, but stderr was empty."
        )

# ---------------------------------------------------------------------------
# 9. JSON Output Format
# ---------------------------------------------------------------------------

class TestJsonOutput:
    def test_exits_zero(self):
        result = _run("--db", _BUNDLED_DB, "--json", _TEST_FA)
        assert result.returncode == 0

    def test_is_valid_json(self):
        import json
        result = _run("--db", _BUNDLED_DB, "--json", _TEST_FA)
        data = json.loads(result.stdout)
        assert isinstance(data, list)
        assert len(data) == 1
        assert "CC" in data[0]
        assert data[0]["NG-MAST/NG-STAR"] == "4186/231"
        assert data[0]["FILE"].endswith("test.fa")

    def test_json_with_comments(self):
        import json
        result = _run("--db", _BUNDLED_DB, "--json", "--comments", _TEST_FA)
        data = json.loads(result.stdout)
        assert len(data) == 1
        assert "penA_comments" in data[0]
        assert "23S_comments" in data[0]
