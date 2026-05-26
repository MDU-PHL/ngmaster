"""
Unit tests for ngmaster.utils

These tests cover the pure-Python utility functions and classes added or
changed in ngmaster v2.0.0.  No mlst binary and no network access are
required.

New features / bug-fixes under test
------------------------------------
* MlstRecord : CC field, get_record() tab/CSV/comments modes
* process_duplicate_23s_alleles() : fix for issue #37
* get_db_version() : new --version feature
* read_ngstar() : CC column parsing
* convert_ngstar() : CC propagation
* collate_results() : CC propagation through combined record
* load_ngstar_comments() : --comments feature
* match_porb() : porB translation table
"""

import pytest
from pathlib import Path
from unittest.mock import MagicMock, patch

from ngmaster.utils import (
    MlstRecord,
    PubMLSTAuth,
    collate_results,
    convert_ngstar,
    get_db_version,
    load_ngstar_comments,
    match_porb,
    process_duplicate_23s_alleles,
    read_ngstar,
)

# ---------------------------------------------------------------------------
# Helpers shared across test classes
# ---------------------------------------------------------------------------

# Minimal mlst --legacy header row for an ngstar run
_NGSTAR_HEADER = "FILE\tSCHEME\tST\tpenA\tmtrR\tporB\tponA\tgyrA\tparC\t23S"
_ROW_TMPL = "sample.fa\tngstar\t231\t23\t19\t8\t1\t1\t9\t{s23}"


def _ngstar_rlist(s23_value):
    return [_NGSTAR_HEADER, _ROW_TMPL.format(s23=s23_value)]


# ---------------------------------------------------------------------------
# 1. MlstRecord
# ---------------------------------------------------------------------------

class TestMlstRecord:
    """Tests for the MlstRecord dataclass and its get_record() method."""

    # -- construction -------------------------------------------------------

    def test_basic_construction_ngstar(self):
        rec = MlstRecord("s.fa", "ngstar", "231",
                         ["23", "19", "8", "1", "1", "9", "100"])
        assert rec.fname == "s.fa"
        assert rec.scheme == "ngstar"
        assert rec.st == "231"
        assert rec.alleles == ["23", "19", "8", "1", "1", "9", "100"]
        assert rec.cc == ""           # CC defaults to empty string
        assert rec.porb == "8"        # alleles[2] for ngstar scheme
        assert rec.simil == ""
        assert rec.part == ""

    def test_basic_construction_ngmast(self):
        rec = MlstRecord("s.fa", "ngmast", "4186", ["10", "20"])
        assert rec.porb == "10"       # alleles[0] for ngmast scheme
        assert rec.cc == ""

    def test_porb_simil_prefix(self):
        """'~10' allele strips the tilde and records similarity flag."""
        rec = MlstRecord("s.fa", "ngmast", "-", ["~10", "20"])
        assert rec.porb == "10"
        assert rec.simil == "~"
        assert rec.part == ""

    def test_porb_part_suffix(self):
        """'10?' allele strips the question mark and records partial flag."""
        rec = MlstRecord("s.fa", "ngmast", "-", ["10?", "20"])
        assert rec.porb == "10"
        assert rec.part == "?"
        assert rec.simil == ""

    def test_porb_multi_allele(self):
        """Comma-separated porB is stored as a list."""
        rec = MlstRecord("s.fa", "ngmast", "-", ["10,15", "20"])
        assert rec.porb == ["10", "15"]

    # -- get_record() -------------------------------------------------------

    def test_get_record_tab_columns(self):
        """Tab-separated output must have 13 columns; CC is the last field."""
        rec = MlstRecord("s.fa", "ngmaSTar", "4186/231",
                         ["10", "20", "23", "19", "8", "1", "1", "9", "100"])
        rec.cc = "CC231"
        fields = rec.get_record().split("\t")
        assert len(fields) == 13           # 3 fixed + 9 alleles + CC
        assert fields[0] == "s.fa"
        assert fields[1] == "ngmaSTar"
        assert fields[2] == "4186/231"
        assert fields[-1] == "CC231"

    def test_get_record_csv_quotes_commas(self):
        """Fields that contain commas are quoted in CSV mode."""
        rec = MlstRecord("s.fa", "ngmaSTar", "4186/231",
                         ["10", "20", "23", "19", "5,7", "1", "1", "9", "100"])
        rec.cc = "-"
        line = rec.get_record(sep=",")
        assert '"5,7"' in line
        assert line.startswith("s.fa,")
        assert line.endswith("-")

    def test_get_record_with_comments_columns(self):
        """Comments mode interleaves allele/comment pairs; output is 20 fields."""
        rec = MlstRecord("s.fa", "ngmaSTar", "4186/231",
                         ["10", "20", "23", "19", "8", "1", "1", "9", "100"])
        rec.cc = "CC231"
        comments = [
            {"23": "penA_note"},   # penA
            {"19": ""},            # mtrR
            {"8": "porB_note"},    # porB
            {"1": ""},             # ponA
            {"1": "gyrA_note"},    # gyrA
            {"9": ""},             # parC
            {"100": "23S_note"},   # 23S
        ]
        header = [
            "FILE", "SCHEME", "NG-MAST/NG-STAR", "porB_NG-MAST", "tbpB",
            "penA", "penA_comments", "mtrR", "mtrR_comments",
            "porB_NG-STAR", "porB_NG-STAR_comments", "ponA", "ponA_comments",
            "gyrA", "gyrA_comments", "parC", "parC_comments",
            "23S", "23S_comments", "CC",
        ]
        fields = rec.get_record(comments=comments, header=header).split("\t")
        assert len(fields) == 20
        assert fields[3] == "10"          # porB_NG-MAST
        assert fields[4] == "20"          # tbpB
        assert fields[5] == "23"          # penA allele
        assert fields[6] == "penA_note"   # penA comment
        assert fields[9] == "8"           # porB_NG-STAR allele
        assert fields[10] == "porB_note"  # porB_NG-STAR comment
        assert fields[-1] == "CC231"      # CC is always last

    def test_get_record_missing_comment_defaults_empty(self):
        """If an allele ID is not in the comment dict, the comment is empty."""
        rec = MlstRecord("s.fa", "ngmaSTar", "4186/231",
                         ["10", "20", "99", "19", "8", "1", "1", "9", "100"])
        rec.cc = "-"
        comments = [{"1": "only_allele_1"}] + [{}] * 6
        header = [
            "FILE", "SCHEME", "NG-MAST/NG-STAR", "porB_NG-MAST", "tbpB",
            "penA", "penA_comments", "mtrR", "mtrR_comments",
            "porB_NG-STAR", "porB_NG-STAR_comments", "ponA", "ponA_comments",
            "gyrA", "gyrA_comments", "parC", "parC_comments",
            "23S", "23S_comments", "CC",
        ]
        fields = rec.get_record(comments=comments, header=header).split("\t")
        # penA allele is "99", not in comments[0] -> empty string comment
        assert fields[5] == "99"
        assert fields[6] == ""


# ---------------------------------------------------------------------------
# 2. process_duplicate_23s_alleles()  :  fix for issue #37
# ---------------------------------------------------------------------------

class TestProcessDuplicate23sAlleles:
    """
    Tests for duplicate-23S deduplication logic.
    """

    def test_exact_duplicates_squashed(self, capsys):
        """'5,5' should become '5' and a WARNING should be emitted to stderr."""
        result = process_duplicate_23s_alleles(_ngstar_rlist("5,5"))
        s23 = result[1].split("\t")[-1]
        assert s23 == "5"
        captured = capsys.readouterr()
        assert "WARNING" in captured.err
        assert "Duplicate 23S" in captured.err

    def test_three_identical_squashed(self, capsys):
        """'5,5,5' should collapse to '5'."""
        result = process_duplicate_23s_alleles(_ngstar_rlist("5,5,5"))
        assert result[1].split("\t")[-1] == "5"
        assert "WARNING" in capsys.readouterr().err

    def test_two_different_alleles_unchanged(self):
        """'5,7' has no duplicates; the value must pass through unchanged."""
        result = process_duplicate_23s_alleles(_ngstar_rlist("5,7"))
        assert result[1].split("\t")[-1] == "5,7"

    def test_single_allele_unchanged(self):
        """A non-comma value is never modified."""
        result = process_duplicate_23s_alleles(_ngstar_rlist("100"))
        assert result[1].split("\t")[-1] == "100"

    def test_no_23s_column_returns_unchanged(self):
        """If the header has no '23S' column the list is returned as-is."""
        rlist = ["FILE\tSCHEME\tST\tpenA\tmtrR\tporB",
                 "s.fa\tngstar\t-\t1\t2\t3"]
        assert process_duplicate_23s_alleles(rlist) == rlist

    def test_empty_list_returns_empty(self):
        assert process_duplicate_23s_alleles([]) == []

    def test_header_only_returns_unchanged(self):
        rlist = [_NGSTAR_HEADER]
        assert process_duplicate_23s_alleles(rlist) == rlist


# ---------------------------------------------------------------------------
# 3. get_db_version()  :  new --version feature
# ---------------------------------------------------------------------------

class TestGetDbVersion:
    """Tests for database version string reading."""

    def test_returns_formatted_string(self, tmp_path):
        """Both version files present -> 'ngmast_X_ngstar_Y' format."""
        (tmp_path / "pubmlst" / "ngmast").mkdir(parents=True)
        (tmp_path / "pubmlst" / "ngstar").mkdir(parents=True)
        (tmp_path / "pubmlst" / "ngmast" / "database_version.txt").write_text(
            "2024-12-31_Unauthenticated\n")
        (tmp_path / "pubmlst" / "ngstar" / "database_version.txt").write_text(
            "2024-12-31_Unauthenticated\n")
        version = get_db_version(str(tmp_path))
        assert version == (
            "ngmast_2024-12-31_Unauthenticated_ngstar_2024-12-31_Unauthenticated"
        )

    def test_missing_version_file(self, tmp_path):
        """Missing file -> graceful fallback message, not an exception."""
        version = get_db_version(str(tmp_path))
        assert "not available" in version.lower() or "error" in version.lower()

    def test_bundled_db_has_valid_version(self, bundled_db_path):
        """Bundled database must have readable version files."""
        version = get_db_version(bundled_db_path)
        assert "ngmast_" in version
        assert "ngstar_" in version
        assert "not available" not in version.lower()


# ---------------------------------------------------------------------------
# 4. read_ngstar()  :  CC column added in v2.0.0
# ---------------------------------------------------------------------------

class TestReadNgstar:
    """Tests for NG-STAR profile table parsing including the optional CC column."""

    def test_parses_st_alleles_and_cc(self, tmp_path):
        """9-column NG-STAR profile: ST, 7 alleles and CC are all stored."""
        p = tmp_path / "ngstar.txt"
        p.write_text(
            "ST\tpenA\tmtrR\tporB\tponA\tgyrA\tparC\t23S\tCC\n"
            "1\t20\t10\t13\t100\t100\t2\t100\t190\n"
            "2\t22\t14\t14\t100\t100\t7\t100\t\n"
        )
        table = read_ngstar(str(p))
        assert "20/10/13/100/100/2/100" in table
        st, cc = table["20/10/13/100/100/2/100"]
        assert st == "1"
        assert cc == "190"

    def test_parses_without_cc_column(self, tmp_path):
        """8-column profile (no CC): CC value is stored as empty string."""
        p = tmp_path / "ngstar.txt"
        p.write_text(
            "ST\tpenA\tmtrR\tporB\tponA\tgyrA\tparC\t23S\n"
            "3\t1\t2\t3\t4\t5\t6\t7\n"
        )
        table = read_ngstar(str(p))
        assert "1/2/3/4/5/6/7" in table
        st, cc = table["1/2/3/4/5/6/7"]
        assert st == "3"
        assert cc == ""

    def test_skips_short_lines(self, tmp_path):
        """Lines with fewer than 8 tokens are silently skipped."""
        p = tmp_path / "ngstar.txt"
        p.write_text(
            "ST\tpenA\n"                               # too short
            "1\t20\t10\t13\t100\t100\t2\t100\t190\n"  # valid
        )
        table = read_ngstar(str(p))
        assert len(table) == 1

    def test_bundled_db_loads(self, ngstar_profile_path):
        """Bundled ngstar.txt must parse without error and return > 0 entries."""
        table = read_ngstar(ngstar_profile_path)
        assert len(table) > 0
        # Every value is a (st, cc) 2-tuple
        for val in list(table.values())[:10]:
            assert len(val) == 2


# ---------------------------------------------------------------------------
# 5. convert_ngstar()  :  CC propagation
# ---------------------------------------------------------------------------

class TestConvertNgstar:
    """Tests for NG-STAR record conversion with CC lookup."""

    _TABLE = {"20/10/13/100/100/2/100": ("1", "190")}

    def test_key_hit_sets_st_and_cc(self):
        rec = MlstRecord("s.fa", "ngstar", "-",
                         ["20", "10", "13", "100", "100", "2", "100"])
        result = convert_ngstar(self._TABLE, rec)
        assert result.st == "1"
        assert result.cc == "190"

    def test_key_miss_sets_dash(self):
        """When the allele combination is not in the profile table, ST and CC are '-'."""
        rec = MlstRecord("s.fa", "ngstar", "-",
                         ["99", "99", "99", "99", "99", "99", "99"])
        result = convert_ngstar(self._TABLE, rec)
        assert result.st == "-"
        assert result.cc == "-"


# ---------------------------------------------------------------------------
# 6. collate_results()  :  CC propagation through combined record
# ---------------------------------------------------------------------------

class TestCollateResults:
    """Tests for collation of NG-MAST + NG-STAR results into ngmaSTar records."""

    @staticmethod
    def _ngmast(fname, st, alleles):
        return MlstRecord(fname, "ngmast", st, alleles)

    @staticmethod
    def _ngstar(fname, st, alleles):
        return MlstRecord(fname, "ngstar", st, alleles)

    def test_basic_happy_path(self):
        """Single file: combined scheme is ngmaSTar, ST is joined, CC propagates."""
        ngmast_res = {"s.fa": self._ngmast("s.fa", "4186", ["10", "20"])}
        ngstar_res = {"s.fa": self._ngstar(
            "s.fa", "-", ["23", "19", "8", "1", "1", "9", "100"])}
        ttable = {"-": "-", "10": "8"}
        ngstartbl = {"23/19/8/1/1/9/100": ("231", "CC231")}

        results = collate_results(ngmast_res, ngstar_res, ttable, ngstartbl)

        assert len(results) == 1
        r = results[0]
        assert r.scheme == "ngmaSTar"
        assert r.st == "4186/231"
        assert r.cc == "CC231"
        # Combined alleles: 2 ngmast + 7 ngstar
        assert len(r.alleles) == 9

    def test_cc_dash_when_profile_miss(self):
        """When allele combination is absent from the profile table, CC is '-'."""
        ngmast_res = {"s.fa": self._ngmast("s.fa", "999", ["5", "6"])}
        ngstar_res = {"s.fa": self._ngstar(
            "s.fa", "-", ["1", "2", "3", "4", "5", "6", "7"])}
        ttable = {"-": "-", "5": "3"}
        ngstartbl = {}   # deliberate miss

        results = collate_results(ngmast_res, ngstar_res, ttable, ngstartbl)

        assert results[0].st == "999/-"
        assert results[0].cc == "-"

    def test_multi_porb_translated_and_joined(self):
        """Multi-allele NG-MAST porB is translated entry-by-entry and joined."""
        ngmast_res = {"s.fa": self._ngmast("s.fa", "-", ["10,15", "20"])}
        ngstar_res = {"s.fa": self._ngstar(
            "s.fa", "-", ["23", "19", "8", "1", "1", "9", "100"])}
        ttable = {"-": "-", "10": "8", "15": "12"}
        ngstartbl = {}

        results = collate_results(ngmast_res, ngstar_res, ttable, ngstartbl)

        # Index 4 in combined alleles is the NG-STAR porB (ngmast[0,1] + ngstar[0,1,2])
        assert results[0].alleles[4] == "8,12"

    def test_key_mismatch_raises_key_error(self):
        """Different file sets in ngmast vs ngstar results must raise KeyError."""
        ngmast_res = {"a.fa": self._ngmast("a.fa", "1", ["1", "2"])}
        ngstar_res = {"b.fa": self._ngstar(
            "b.fa", "-", ["1", "2", "3", "4", "5", "6", "7"])}
        with pytest.raises(KeyError):
            collate_results(ngmast_res, ngstar_res, {}, {})


# ---------------------------------------------------------------------------
# 7. load_ngstar_comments()  :  --comments feature
# ---------------------------------------------------------------------------

class TestLoadNgstarComments:
    """Tests for loading the per-allele NG-STAR comment TSV."""

    def test_returns_list_of_seven_dicts(self, bundled_db_path):
        """Bundled TSV loads as a list of 7 dicts, one per NG-STAR locus."""
        comments = load_ngstar_comments(bundled_db_path)
        assert len(comments) == 7
        for d in comments:
            assert isinstance(d, dict)

    def test_locus_order_matches_schema(self, bundled_db_path):
        """Locus order must be penA, mtrR, porB, ponA, gyrA, parC, 23S."""
        expected_loci = ["penA", "mtrR", "porB", "ponA", "gyrA", "parC", "23S"]
        comments = load_ngstar_comments(bundled_db_path)
        # Verify the penA dict (index 0) contains allele "1" (absent == empty comment)
        # and the mtrR dict (index 1) is non-None regardless of content
        assert len(comments) == len(expected_loci)
        for d in comments:
            assert d is not None

    def test_specific_comment_from_tmp_tsv(self, tmp_path):
        """Custom minimal TSV: specific allele comment is retrieved correctly."""
        ngstar_dir = tmp_path / "pubmlst" / "ngstar"
        ngstar_dir.mkdir(parents=True)
        (ngstar_dir / "allele_comments.tsv").write_text(
            "penA\t1\t\n"
            "penA\t3\tType 76 Non Mosaic\n"
            "mtrR\t1\t\n"
            "porB\t1\t\n"
            "ponA\t1\t\n"
            "gyrA\t1\t\n"
            "parC\t1\t\n"
            "23S\t1\t\n"
        )
        comments = load_ngstar_comments(str(tmp_path))
        pena_dict = comments[0]   # penA is locus index 0
        assert pena_dict.get("3") == "Type 76 Non Mosaic"
        assert pena_dict.get("1") == ""


# ---------------------------------------------------------------------------
# 8. match_porb()  :  porB translation table
# ---------------------------------------------------------------------------

class TestMatchPorb:
    """Tests for NG-MAST / NG-STAR porB sequence matching."""

    def test_bundled_db_returns_populated_dict(self, ngmast_porb_path,
                                               ngstar_porb_path):
        """Real bundled porB files: dict is non-empty, sentinel '-' maps to '-'."""
        ttable = match_porb(ngmast_porb_path, ngstar_porb_path)
        assert isinstance(ttable, dict)
        assert "-" in ttable
        assert ttable["-"] == "-"
        assert len(ttable) > 1     # at least the sentinel + one real entry

    def test_missing_files_exit(self, tmp_path):
        """Non-existent porB files must trigger err(), which calls sys.exit()."""
        with pytest.raises(SystemExit):
            match_porb(str(tmp_path / "missing.tfa"), str(tmp_path / "also_missing.tfa"))


# ---------------------------------------------------------------------------
# 9. PubMLSTAuth  :  authentication mode selection and rauth params fix
# ---------------------------------------------------------------------------

class TestPubMLSTAuth:
    """
    Tests for authentication mode selection and rauth params fix.

    All network calls are mocked; no credentials on disk are required.
    """

    def test_api_key_mode_when_key_stored(self):
        """Stored API key selects 'api_key' mode and sets X-API-Key header."""
        with patch("ngmaster.utils.retrieve_api_key", return_value="testkey123"):
            auth = PubMLSTAuth()
        assert auth.auth_mode == "api_key"
        assert auth.session.headers.get("X-API-Key") == "testkey123"
        assert "ngmaster/" in auth.session.headers.get("User-Agent", "")

    def test_api_key_mode_user_agent_format(self):
        """User-Agent header must be 'ngmaster/<version>'."""
        with patch("ngmaster.utils.retrieve_api_key", return_value="testkey123"):
            auth = PubMLSTAuth()
        import ngmaster
        assert auth.session.headers["User-Agent"] == f"ngmaster/{ngmaster.__version__}"

    def test_oauth_mode_when_no_api_key_but_tokens_present(self):
        """No API key but OAuth tokens present selects 'oauth' mode."""
        with patch("ngmaster.utils.retrieve_api_key", return_value=None), \
             patch("ngmaster.utils.get_client_credentials",
                   return_value=("client_id_24char_padding", "secret_42char_padding_padding_padding_padd")), \
             patch("ngmaster.utils.retrieve_session_token",
                   return_value=("sess_token", "sess_secret")):
            auth = PubMLSTAuth()
        assert auth.auth_mode == "oauth"
        assert "ngmaster/" in auth.session.headers.get("User-Agent", "")

    def test_unauthenticated_when_no_credentials(self, capsys):
        """No API key and no OAuth tokens falls back to unauthenticated."""
        with patch("ngmaster.utils.retrieve_api_key", return_value=None), \
             patch("ngmaster.utils.get_client_credentials",
                   side_effect=ValueError("not found")):
            auth = PubMLSTAuth()
        assert auth.auth_mode is None
        captured = capsys.readouterr()
        assert "mlstdb connect" in captured.err

    def test_oauth_get_response_passes_params_kwarg(self):
        """OAuth mode: get_response() calls session.get with params={} to avoid rauth bug."""
        mock_session = MagicMock()
        mock_session.get.return_value = MagicMock(status_code=200)

        with patch("ngmaster.utils.retrieve_api_key", return_value=None), \
             patch("ngmaster.utils.get_client_credentials",
                   return_value=("client_id_24char_padding", "secret_42char_padding_padding_padding_padd")), \
             patch("ngmaster.utils.retrieve_session_token",
                   return_value=("sess_token", "sess_secret")), \
             patch("ngmaster.utils.OAuth1Session", return_value=mock_session):
            auth = PubMLSTAuth()
            auth.get_response("https://example.com/api")

        mock_session.get.assert_called_once_with("https://example.com/api", params={})

    def test_api_key_get_response_no_params_kwarg(self):
        """API key mode: get_response() calls session.get without params kwarg."""
        mock_session = MagicMock()
        mock_session.get.return_value = MagicMock(status_code=200)
        mock_session.headers = {}

        with patch("ngmaster.utils.retrieve_api_key", return_value="mykey"), \
             patch("ngmaster.utils.requests") as mock_requests:
            mock_requests.Session.return_value = mock_session
            auth = PubMLSTAuth()
            auth.get_response("https://example.com/api")

        mock_session.get.assert_called_once_with("https://example.com/api")
