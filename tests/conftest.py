"""
Shared pytest fixtures for the ngmaster test suite.
"""

import os
import shutil
import sys
from pathlib import Path

import pytest

import ngmaster

_BUNDLE_ROOT = Path(ngmaster.__file__).parent
_TEST_DATA = _BUNDLE_ROOT / "test"

# ---------------------------------------------------------------------------
# Ensure mlst is discoverable by subprocess calls.
# When pytest is launched via `pixi run pytest` the pixi env bin/ is already
# on PATH.  When launched via a raw python interpreter (e.g. inside VS Code)
# the pixi bin/ may be absent.  We add it here so that subprocess.run(['mlst'])
# works in both scenarios.
# ---------------------------------------------------------------------------

def _extend_path_with_pixi():
    pixi_bin = Path(sys.executable).parent
    current = os.environ.get("PATH", "")
    if str(pixi_bin) not in current.split(os.pathsep):
        os.environ["PATH"] = str(pixi_bin) + os.pathsep + current


_extend_path_with_pixi()


# ---------------------------------------------------------------------------
# Path fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def bundled_db_path():
    """Absolute path to the bundled ngmaster database directory."""
    return str(_BUNDLE_ROOT / "db")


@pytest.fixture(scope="session")
def test_fa_path():
    """Path to the bundled test.fa (expected result: NG-MAST 4186 / NG-STAR 231)."""
    return str(_TEST_DATA / "test.fa")


@pytest.fixture(scope="session")
def noseq_fa_path():
    """Path to noseq.fa – FASTA with a header but no sequence data."""
    return str(_TEST_DATA / "noseq.fa")


@pytest.fixture(scope="session")
def null_fa_path():
    """Path to null.fa – completely empty file."""
    return str(_TEST_DATA / "null.fa")


@pytest.fixture(scope="session")
def ngstar_profile_path(bundled_db_path):
    """Absolute path to the bundled ngstar.txt profile table."""
    return str(Path(bundled_db_path) / "pubmlst" / "ngstar" / "ngstar.txt")


@pytest.fixture(scope="session")
def ngmast_porb_path(bundled_db_path):
    """Absolute path to the bundled NG-MAST porB.tfa allele file."""
    return str(Path(bundled_db_path) / "pubmlst" / "ngmast" / "porB.tfa")


@pytest.fixture(scope="session")
def ngstar_porb_path(bundled_db_path):
    """Absolute path to the bundled NG-STAR porB.tfa allele file."""
    return str(Path(bundled_db_path) / "pubmlst" / "ngstar" / "porB.tfa")


@pytest.fixture(scope="session")
def comments_tsv_path(bundled_db_path):
    """Absolute path to the bundled allele_comments.tsv file."""
    return str(Path(bundled_db_path) / "pubmlst" / "ngstar" / "allele_comments.tsv")


# ---------------------------------------------------------------------------
# CLI fixture
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def ngmaster_cmd():
    """
    Path to the ngmaster CLI entry point in the active Python environment.
    Resolves via shutil.which first, then falls back to the bin/ sibling of the
    current Python interpreter (works in venvs, conda envs, and pixi envs).
    """
    cmd = shutil.which("ngmaster")
    if cmd is None:
        cmd = str(Path(sys.executable).parent / "ngmaster")
    return cmd
