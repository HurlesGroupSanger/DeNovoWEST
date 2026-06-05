"""Tests for denovowest.annotations.annotate_dbnsfp."""

from pathlib import Path
from unittest.mock import MagicMock

import pandas as pd
import pytest
from click.testing import CliRunner

from denovowest.annotations.annotate_dbnsfp import (
    get_gene_based_score,
    get_transcript_based_score,
    get_transcripts,
    max_value,
    parse,
)
from denovowest.cli import main

#############################
# Paths to fixtures         #
#############################

_ROOT = Path(__file__).parent.parent.parent
_INPUT = _ROOT / "examples" / "annotation" / "dbNSFP" / "input"
_OUTPUT = _ROOT / "examples" / "annotation" / "dbNSFP" / "output"

VARIANTS = _INPUT / "OR4F5_dnms.tsv"
DBNSFP = _INPUT / "OR4F5_dbNSFP5.3.1a_grch38.gz"
GFF = _INPUT / "OR4F5.gff"
DBNSFP_EXPECTED = _OUTPUT / "OR4F5_dbnsfp.tsv"


#############################
# Smoke tests               #
#############################


@pytest.mark.smoke
def test_annotate_dbnsfp_smoke(tmp_path):
    """CLI exits cleanly and appends the AlphaMissense_score column."""
    output = tmp_path / "out.tsv"
    runner = CliRunner()

    result = runner.invoke(
        main,
        ["annotate-dbnsfp", str(VARIANTS), str(DBNSFP), str(output),
         "-c", "AlphaMissense_score", "--gff", str(GFF)],
    )

    assert result.exit_code == 0, result.output
    df = pd.read_csv(output, sep="\t")
    assert "AlphaMissense_score" in df.columns
    assert "transcript_in_dbnsfp" in df.columns
    assert len(df) == 11


@pytest.mark.smoke
def test_annotate_dbnsfp_smoke_no_gff(tmp_path):
    """CLI exits cleanly without a GFF (coordinate + gene match only, no transcript filtering)."""
    output = tmp_path / "out.tsv"
    runner = CliRunner()

    result = runner.invoke(
        main,
        ["annotate-dbnsfp", str(VARIANTS), str(DBNSFP), str(output),
         "-c", "AlphaMissense_score"],
    )

    assert result.exit_code == 0, result.output
    df = pd.read_csv(output, sep="\t")
    assert "AlphaMissense_score" in df.columns
    assert len(df) == 11


#############################
# Regression tests          #
#############################


@pytest.mark.regression
def test_annotate_dbnsfp_regression(tmp_path):
    """Output matches the known-good reference."""
    output = tmp_path / "out.tsv"
    runner = CliRunner()
    runner.invoke(
        main,
        ["annotate-dbnsfp", str(VARIANTS), str(DBNSFP), str(output),
         "-c", "AlphaMissense_score", "--gff", str(GFF)],
    )

    actual = pd.read_csv(output, sep="\t", dtype=str)
    expected = pd.read_csv(DBNSFP_EXPECTED, sep="\t", dtype=str)
    pd.testing.assert_frame_equal(actual, expected)


#########################################
# Helpers for unit tests                #
#########################################

# In test records, gene/transcript IDs are placed at these indices
_GENE_IDX = 13
_TRANSCRIPT_IDX = 14


def _make_record(
    chrom: str = "1",
    pos: str = "65518",
    ref: str = "A",
    alt: str = "T",
    gene_ids: str = "ENSG00000186092",
    transcript_ids: str = "ENST00000335137",
    extra: dict | None = None,
    n_fields: int = 20,
) -> str:
    """Build a minimal tab-separated dbNSFP record string."""
    fields = ["."] * n_fields
    fields[0] = chrom
    fields[1] = pos
    fields[2] = ref
    fields[3] = alt
    fields[_GENE_IDX] = gene_ids
    fields[_TRANSCRIPT_IDX] = transcript_ids
    if extra:
        for idx, val in extra.items():
            fields[idx] = val
    return "\t".join(fields)


def _make_mock_gff_db(transcript_ids: list[str]) -> MagicMock:
    """Return a minimal gffutils database mock with a single gene and its transcripts."""
    mock_transcript = MagicMock()
    mock_transcript.__getitem__ = MagicMock(return_value=transcript_ids)

    mock_db = MagicMock()
    mock_db.__getitem__ = MagicMock(return_value=MagicMock())  # gene feature
    mock_db.children = MagicMock(return_value=[mock_transcript])
    return mock_db


#########################################
# Unit – parse                          #
#########################################


@pytest.mark.unit
def test_parse_returns_record_for_matching_gene():
    """parse returns a populated dict when the gene ID matches."""
    record = _make_record()
    result = parse(record, [], "ENSG00000186092", ensembl_geneid_idx=_GENE_IDX, ensembl_transcriptid_idx=_TRANSCRIPT_IDX)

    assert result["chrom"] == "1"
    assert result["pos"] == 65518
    assert result["ref"] == "A"
    assert result["alt"] == "T"


@pytest.mark.unit
def test_parse_returns_empty_dict_for_mismatched_gene():
    """parse returns an empty dict when the record belongs to a different gene."""
    record = _make_record(gene_ids="ENSG00000999999")
    result = parse(record, [], "ENSG00000186092", ensembl_geneid_idx=_GENE_IDX, ensembl_transcriptid_idx=_TRANSCRIPT_IDX)

    assert result == {}


@pytest.mark.unit
def test_parse_strips_version_from_query_gene_id():
    """parse matches versioned query IDs (e.g. ENSG00000186092.7) against unversioned dbNSFP IDs."""
    record = _make_record(gene_ids="ENSG00000186092")
    result = parse(record, [], "ENSG00000186092.7", ensembl_geneid_idx=_GENE_IDX, ensembl_transcriptid_idx=_TRANSCRIPT_IDX)

    assert result != {}
    assert result["pos"] == 65518


@pytest.mark.unit
def test_parse_extracts_score_columns():
    """parse populates the dict with the requested column indices."""
    record = _make_record(extra={5: "0.95", 6: "0.80"})
    result = parse(record, [5, 6], "ENSG00000186092", ensembl_geneid_idx=_GENE_IDX, ensembl_transcriptid_idx=_TRANSCRIPT_IDX)

    assert result["5"] == "0.95"
    assert result["6"] == "0.80"


@pytest.mark.unit
def test_parse_uses_transcript_based_score_when_gff_provided():
    """With transcript filtering, parse picks the score for the matching transcript."""
    # Two transcripts; only the first is in the GFF
    record = _make_record(
        gene_ids="ENSG00000186092",
        transcript_ids="ENST00000335137;ENST00000999999",
        extra={5: "0.90;0.20"},  # scores align with transcript order
    )
    result = parse(record, [5], "ENSG00000186092", transcript_ids_gff=["ENST00000335137"], ensembl_geneid_idx=_GENE_IDX, ensembl_transcriptid_idx=_TRANSCRIPT_IDX)

    # max of the GFF-matched transcript score → 0.90
    assert result["5"] == pytest.approx(0.90)


@pytest.mark.unit
def test_parse_falls_back_to_gene_score_when_transcript_absent():
    """Falls back to a gene-level score when the GFF transcript is not in dbNSFP."""
    # dbNSFP stores one gene_id entry per transcript (separated by ";")
    record = _make_record(
        gene_ids="ENSG00000186092;ENSG00000186092",
        transcript_ids="ENST00000111111;ENST00000222222",  # neither in GFF
        extra={5: "0.70;0.80"},
    )
    result = parse(record, [5], "ENSG00000186092", transcript_ids_gff=["ENST00000335137"], ensembl_geneid_idx=_GENE_IDX, ensembl_transcriptid_idx=_TRANSCRIPT_IDX)

    # Falls back to gene-level → max across all positions for this gene = 0.80
    assert result["5"] == pytest.approx(0.80)


#########################################
# Unit – max_value                      #
#########################################


@pytest.mark.unit
def test_max_value_ignores_missing_sentinel():
    """max_value skips "." entries and returns the numeric maximum."""
    assert max_value(["1.5", ".", "2.3"]) == pytest.approx(2.3)


@pytest.mark.unit
def test_max_value_all_missing_returns_dot():
    """max_value returns "." when every entry is missing."""
    assert max_value([".", "."]) == "."


@pytest.mark.unit
def test_max_value_nonnumeric_joins_values():
    """max_value joins non-numeric entries with semicolons."""
    assert max_value(["pathogenic", "benign"]) == "pathogenic;benign"


#########################################
# Unit – get_transcript_based_score     #
#########################################


@pytest.mark.unit
def test_get_transcript_based_score_returns_max_for_matched_transcripts():
    scores = ["0.5", "0.9", "0.3"]
    gff_transcripts = ["ENST00000002"]
    dbnsfp_transcripts = ["ENST00000001", "ENST00000002", "ENST00000003"]

    result = get_transcript_based_score(scores, gff_transcripts, dbnsfp_transcripts)

    assert result == pytest.approx(0.9)


@pytest.mark.unit
def test_get_gene_based_score_returns_max_for_matched_gene():
    scores = ["0.5", "0.9", "0.3"]
    gene_ids_dbnsfp = ["ENSG00000001", "ENSG00000001", "ENSG00000002"]

    result = get_gene_based_score(scores, "ENSG00000001", gene_ids_dbnsfp)

    assert result == pytest.approx(0.9)


#########################################
# Unit – get_transcripts                #
#########################################


@pytest.mark.unit
def test_get_transcripts_versioned_gene_id():
    """Versioned gene IDs are looked up directly in the GFF database."""
    mock_db = _make_mock_gff_db(["ENST00000335137.5", "ENST00000414099.2"])

    result = get_transcripts("ENSG00000186092.7", mock_db, {})

    # Version suffix is stripped from transcript IDs
    assert result == ["ENST00000335137", "ENST00000414099"]
    mock_db.__getitem__.assert_called_once_with("ENSG00000186092.7")


@pytest.mark.unit
def test_get_transcripts_unversioned_gene_id_uses_mapping():
    """Unversioned gene IDs are resolved via the ensembl_gene_id_map_version dict."""
    mock_db = _make_mock_gff_db(["ENST00000335137.5"])
    version_map = {"ENSG00000186092": "ENSG00000186092.7"}

    result = get_transcripts("ENSG00000186092", mock_db, version_map)

    assert result == ["ENST00000335137"]
    mock_db.__getitem__.assert_called_once_with("ENSG00000186092.7")


@pytest.mark.unit
def test_get_transcripts_missing_gene_returns_empty_list():
    """Returns an empty list without raising when the gene is absent from the GFF."""
    mock_db = MagicMock()
    mock_db.__getitem__ = MagicMock(side_effect=KeyError("ENSG00000000000"))

    result = get_transcripts("ENSG00000000000", mock_db, {})

    assert result == []
