"""Tests for denovowest.annotations.annotate_cadd."""

from pathlib import Path
from unittest.mock import MagicMock

import pandas as pd
import pytest
from click.testing import CliRunner

from denovowest.annotations.annotate_cadd import load_cadd
from denovowest.cli import main

#############################
# Paths to fixtures         #
#############################

_ROOT = Path(__file__).parent.parent.parent
_INPUT = _ROOT / "examples" / "annotation" / "CADD" / "input"
_OUTPUT = _ROOT / "examples" / "annotation" / "CADD" / "output"

VARIANTS = _INPUT / "OR4F5_dnms.tsv"
CADD = _INPUT / "OR4F5.cadd.tsv.gz"
CADD_EXPECTED = _OUTPUT / "OR4F5_dnms_cadd.tsv"

_EXPECTED_COLUMNS = ["gene_id", "chrom", "pos", "ref", "alt", "region_type", "raw", "score"]


#############################
# Smoke tests               #
#############################


@pytest.mark.smoke
def test_annotate_cadd_smoke(tmp_path):
    """CLI exits cleanly and produces the expected columns and row count."""
    output = tmp_path / "out.tsv"
    runner = CliRunner()

    result = runner.invoke(main, ["annotate-cadd", str(VARIANTS), str(CADD), str(output)])

    assert result.exit_code == 0, result.output
    df = pd.read_csv(output, sep="\t")
    assert list(df.columns) == _EXPECTED_COLUMNS
    assert len(df) == 11


#############################
# Regression tests          #
#############################


@pytest.mark.regression
def test_annotate_cadd_regression(tmp_path):
    """Output matches the known-good reference."""
    output = tmp_path / "out.tsv"
    runner = CliRunner()
    runner.invoke(main, ["annotate-cadd", str(VARIANTS), str(CADD), str(output)])

    actual = pd.read_csv(output, sep="\t")
    expected = pd.read_csv(CADD_EXPECTED, sep="\t")
    pd.testing.assert_frame_equal(actual, expected)


#############################
# Unit tests                #
#############################


@pytest.mark.unit
def test_load_cadd_parses_records_correctly():
    """load_cadd converts tab-delimited CADD lines into a correctly typed DataFrame."""
    mock_tabix = MagicMock()
    mock_tabix.fetch.return_value = [
        "1\t65518\tA\tT\t2.232995\t17.82",
        "1\t65549\tG\tT\t0.362594\t3.97",
    ]

    df = load_cadd(mock_tabix, "1", 65517, 65549)

    assert list(df.columns) == ["chrom", "pos", "ref", "alt", "raw", "score"]
    assert len(df) == 2
    assert df.iloc[0]["pos"] == 65518
    assert df.iloc[0]["raw"] == pytest.approx(2.232995)
    assert df.iloc[0]["score"] == pytest.approx(17.82)
    assert df.iloc[1]["ref"] == "G"


@pytest.mark.unit
def test_annotate_cadd_empty_input_exits_cleanly(tmp_path):
    """Empty rates file: exit 0, output has the expected columns but no data rows."""
    empty_input = tmp_path / "empty.tsv"
    empty_input.write_text("gene_id\tchrom\tpos\tref\talt\tregion_type\n")
    output = tmp_path / "out.tsv"
    runner = CliRunner()

    result = runner.invoke(main, ["annotate-cadd", str(empty_input), str(CADD), str(output)])

    assert result.exit_code == 0, result.output
    df = pd.read_csv(output, sep="\t")
    assert len(df) == 0
    assert "raw" in df.columns
    assert "score" in df.columns


@pytest.mark.unit
def test_annotate_cadd_no_chr_prefix(tmp_path):
    """Rates file without "chr" prefix: output chromosomes also have no prefix."""
    no_chr_input = tmp_path / "nochr.tsv"
    # Rewrite the example file stripping the chr prefix
    df = pd.read_csv(VARIANTS, sep="\t", dtype=str)
    df["chrom"] = df["chrom"].str.replace("chr", "", regex=False)
    df.to_csv(no_chr_input, sep="\t", index=False)
    output = tmp_path / "out.tsv"
    runner = CliRunner()

    result = runner.invoke(main, ["annotate-cadd", str(no_chr_input), str(CADD), str(output)])

    assert result.exit_code == 0, result.output
    out_df = pd.read_csv(output, sep="\t", dtype=str)
    assert not out_df["chrom"].str.startswith("chr").any()
