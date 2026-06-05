"""Tests for denovowest.annotations.annotate_custom."""

import gzip
from pathlib import Path
from unittest.mock import MagicMock

import pandas as pd
import pytest
from click.testing import CliRunner

from denovowest.annotations.annotate_custom import (
    _parse,
    check_columns,
    extract_columns_indices,
)
from denovowest.cli import main

#############################
# Paths to fixtures         #
#############################

_ROOT = Path(__file__).parent.parent.parent
_INPUT = _ROOT / "examples" / "annotation" / "custom" / "input"
_INPUT_CADD = _ROOT / "examples" / "annotation" / "CADD" / "input"
_OUTPUT = _ROOT / "examples" / "annotation" / "custom" / "output"

VARIANTS = _INPUT / "OR4F5_dnms.tsv"
ANNOTATION = _INPUT_CADD / "OR4F5.cadd.tsv.gz"
CUSTOM_EXPECTED = _OUTPUT / "OR4F5_dnms_custom.tsv"


#############################
# Smoke tests               #
#############################


@pytest.mark.smoke
def test_annotate_custom_smoke(tmp_path):
    """CLI exits cleanly and appends the requested annotation column."""
    output = tmp_path / "out.tsv"
    runner = CliRunner()

    result = runner.invoke(
        main, ["annotate-custom", str(VARIANTS), str(ANNOTATION), str(output), "-c", "score"]
    )

    assert result.exit_code == 0, result.output
    df = pd.read_csv(output, sep="\t")
    assert "score" in df.columns
    assert len(df) == 11


#############################
# Regression tests          #
#############################


@pytest.mark.regression
def test_annotate_custom_regression(tmp_path):
    """Output matches the known-good reference."""
    output = tmp_path / "out.tsv"
    runner = CliRunner()
    runner.invoke(
        main, ["annotate-custom", str(VARIANTS), str(ANNOTATION), str(output), "-c", "score"]
    )

    actual = pd.read_csv(output, sep="\t", dtype=str)
    expected = pd.read_csv(CUSTOM_EXPECTED, sep="\t", dtype=str)
    pd.testing.assert_frame_equal(actual, expected)


#############################
# Unit – _parse             #
#############################


@pytest.mark.unit
def test_parse_extracts_chrom_pos_ref_alt():
    """_parse extracts coordinates and no extra columns when columns list is empty."""
    line = "chr1\t100\tA\tT\t0.5\t17.82"
    result = _parse(line, [])

    assert result == {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T"}


@pytest.mark.unit
def test_parse_extracts_requested_columns():
    """_parse adds the requested column indices to the returned dict."""
    line = "chr1\t100\tA\tT\t0.5\t17.82"
    result = _parse(line, [4, 5])

    assert result["4"] == "0.5"
    assert result["5"] == "17.82"


@pytest.mark.unit
def test_parse_casts_pos_to_int():
    """_parse stores pos as int, not string."""
    line = "1\t65518\tA\tT\t2.23"
    result = _parse(line, [])

    assert isinstance(result["pos"], int)


#######################################
# Unit – extract_columns_indices      #
#######################################


@pytest.mark.unit
def test_extract_columns_indices_resolves_names_to_positions():
    """Column names are correctly mapped to their tab-delimited indices."""
    mock_tabix = MagicMock()
    mock_tabix.header = ["#chrom\tpos\tref\talt\traw\tscore"]

    indices, names = extract_columns_indices(mock_tabix, "score,raw", "")

    assert names == ["raw", "score"]
    assert indices == [4, 5]


@pytest.mark.unit
def test_extract_columns_indices_gzip_fallback(tmp_path):
    """Falls back to reading the gzip header when the tabix header list is empty."""
    # Annotation file whose first line is NOT a # comment → tabix header is empty
    gz_file = tmp_path / "no_hash_header.tsv.gz"
    with gzip.open(gz_file, "wt") as f:
        f.write("chrom\tpos\tref\talt\tVEST4_score\n")
        f.write("1\t100\tA\tT\t0.9\n")

    mock_tabix = MagicMock()
    mock_tabix.header = []  # IndexError on header[0] triggers the fallback
    mock_tabix.filename = bytes(str(gz_file), "utf-8")

    indices, names = extract_columns_indices(mock_tabix, "VEST4_score", "")

    assert names == ["VEST4_score"]
    assert indices == [4]


@pytest.mark.unit
def test_extract_columns_indices_defaults_to_all_beyond_coords():
    """With no column filter, every column beyond the four coordinate columns is returned."""
    mock_tabix = MagicMock()
    mock_tabix.header = ["#chrom\tpos\tref\talt\traw\tscore\textra"]

    indices, names = extract_columns_indices(mock_tabix, "", "")

    assert names == ["raw", "score", "extra"]
    assert indices == [4, 5, 6]


#############################
# Unit – check_columns      #
#############################


@pytest.mark.unit
def test_check_columns_detects_name_clash():
    """Returns the set of column names present in both the variant file and annotation."""
    shared = check_columns(["gene_id", "chrom", "score"], ["score", "VEST4"])

    assert shared == {"score"}


@pytest.mark.unit
def test_check_columns_returns_empty_set_when_no_clash():
    """Returns an empty set when variant file and annotation have no columns in common."""
    shared = check_columns(["gene_id", "chrom", "pos"], ["score", "VEST4"])

    assert shared == set()
