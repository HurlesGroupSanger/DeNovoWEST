"""Tests for denovowest.annotations.annotate_vcf."""

from pathlib import Path

import pandas as pd
import pytest
from click.testing import CliRunner

from denovowest.annotations.annotate_vcf import (
    _parse,
    build_gene_mapping_from_column,
    format_gene_annotation,
)
from denovowest.cli import main

#############################
# Paths to fixtures         #
#############################

_ROOT = Path(__file__).parent.parent.parent
_INPUT = _ROOT / "examples" / "annotation" / "vcf" / "input"
_OUTPUT = _ROOT / "examples" / "annotation" / "vcf" / "output"

VARIANTS = _INPUT / "OR4F5_dnms.tsv"
VCF = _INPUT / "OR4F5_popevemap.vcf.gz"
GFF = _INPUT / "OR4F5.gff"
VCF_EXPECTED = _OUTPUT / "OR4F5_vcf.tsv"

_COLUMNS = "-c popEVE,EVE,ESM1v".split()
_EXPECTED_COLUMNS = ["gene_id", "chrom", "pos", "ref", "alt", "region_type", "popEVE", "EVE", "ESM1v"]


#############################
# Smoke tests               #
#############################


@pytest.mark.smoke
def test_annotate_vcf_smoke(tmp_path):
    """CLI exits cleanly and appends the requested INFO columns."""
    output = tmp_path / "out.tsv"
    runner = CliRunner()

    result = runner.invoke(main, ["annotate-vcf", str(VARIANTS), str(VCF), str(output)] + _COLUMNS)

    assert result.exit_code == 0, result.output
    df = pd.read_csv(output, sep="\t")
    assert list(df.columns) == _EXPECTED_COLUMNS
    assert len(df) == 11


@pytest.mark.smoke
def test_annotate_vcf_match_gene_smoke(tmp_path):
    """--match-gene with a GFF exits cleanly and resolves versioned gene IDs to symbols."""
    output = tmp_path / "out.tsv"
    runner = CliRunner()

    result = runner.invoke(
        main,
        ["annotate-vcf", str(VARIANTS), str(VCF), str(output)]
        + _COLUMNS
        + ["--match-gene", "--gff", str(GFF)],
    )

    assert result.exit_code == 0, result.output
    df = pd.read_csv(output, sep="\t")
    assert list(df.columns) == _EXPECTED_COLUMNS
    assert len(df) == 11


#############################
# Regression tests          #
#############################


@pytest.mark.regression
def test_annotate_vcf_regression(tmp_path):
    """Coordinate-only output matches the known-good reference."""
    output = tmp_path / "out.tsv"
    runner = CliRunner()
    runner.invoke(main, ["annotate-vcf", str(VARIANTS), str(VCF), str(output)] + _COLUMNS)

    actual = pd.read_csv(output, sep="\t", dtype=str)
    expected = pd.read_csv(VCF_EXPECTED, sep="\t", dtype=str)
    pd.testing.assert_frame_equal(actual, expected)



#########################################
# Helpers                               #
#########################################


class _MockVCFRecord:
    """Minimal pysam VCF record stand-in."""

    def __init__(self, chrom, pos, ref, alt, info):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alts = (alt,)
        self.info = info


def _make_annotation_df(
    chroms=("1", "1"),
    positions=(100, 200),
    refs=("A", "G"),
    alts=("T", "C"),
    extra: dict | None = None,
) -> pd.DataFrame:
    """Build a minimal annotation block DataFrame."""
    df = pd.DataFrame(
        {"chrom": list(chroms), "pos": list(positions), "ref": list(refs), "alt": list(alts)}
    )
    if extra:
        for col, vals in extra.items():
            df[col] = vals
    return df


#########################################
# Unit – _parse                         #
#########################################


@pytest.mark.unit
def test_parse_extracts_coordinates_and_info_fields():
    """_parse turns a pysam record into a flat dict with coordinates + INFO fields."""
    record = _MockVCFRecord("1", 100, "A", "T", {"score": 1.5, "GENE": "BRCA1"})
    result = _parse(record)

    assert result["chrom"] == "1"
    assert result["pos"] == 100
    assert result["ref"] == "A"
    assert result["alt"] == "T"
    assert result["score"] == 1.5
    assert result["GENE"] == "BRCA1"


@pytest.mark.unit
def test_parse_uses_first_alt_allele():
    """_parse takes only the first ALT allele from a multi-allelic record."""
    record = _MockVCFRecord("1", 100, "A", "T", {})
    result = _parse(record)

    assert result["alt"] == "T"


#########################################
# Unit – build_gene_mapping_from_column #
#########################################


@pytest.mark.unit
def test_build_gene_mapping_from_column_basic():
    """Returns a dict mapping each unique gene_id to its symbol."""
    df = pd.DataFrame(
        {
            "gene_id": ["ENSG00000001", "ENSG00000001", "ENSG00000002"],
            "symbol": ["BRCA1", "BRCA1", "TP53"],
            "pos": [100, 200, 300],
        }
    )
    mapping = build_gene_mapping_from_column(df)

    assert mapping == {"ENSG00000001": "BRCA1", "ENSG00000002": "TP53"}


@pytest.mark.unit
def test_build_gene_mapping_from_column_deduplicates():
    """Each gene_id appears exactly once in the output mapping."""
    df = pd.DataFrame({"gene_id": ["ENSG1", "ENSG1"], "symbol": ["GeneA", "GeneA"]})
    mapping = build_gene_mapping_from_column(df)

    assert len(mapping) == 1


#########################################
# Unit – format_gene_annotation         #
#########################################


@pytest.mark.unit
def test_format_gene_annotation_adds_chr_prefix():
    """When chr_prefixed=True the chrom column gains the "chr" prefix."""
    blocks = [_make_annotation_df()]
    result = format_gene_annotation("ENSG00000001", blocks, False, {}, True, "", "")

    assert (result["chrom"] == "chr1").all()


@pytest.mark.unit
def test_format_gene_annotation_keeps_no_prefix():
    """When chr_prefixed=False the chrom column is left as-is."""
    blocks = [_make_annotation_df()]
    result = format_gene_annotation("ENSG00000001", blocks, False, {}, False, "", "")

    assert (result["chrom"] == "1").all()


@pytest.mark.unit
def test_format_gene_annotation_deduplicates_by_pos_ref_alt():
    """Duplicate (pos, ref, alt) rows are reduced to one entry."""
    dup_df = _make_annotation_df(
        chroms=("1", "1"),
        positions=(100, 100),
        refs=("A", "A"),
        alts=("T", "T"),
    )
    result = format_gene_annotation("ENSG00000001", [dup_df], False, {}, False, "", "")

    assert len(result) == 1


@pytest.mark.unit
def test_format_gene_annotation_filters_by_ensembl_id():
    """With match_gene=True and gene_content='ensembl_id', keeps only rows whose GENE field matches."""
    df = _make_annotation_df(extra={"GENE": ["ENSG00000001", "ENSG00000002"]})
    result = format_gene_annotation(
        "ENSG00000001", [df], True, {}, False, "GENE", "ensembl_id"
    )

    assert len(result) == 1
    assert result.iloc[0]["pos"] == 100


@pytest.mark.unit
def test_format_gene_annotation_filters_by_symbol():
    """With match_gene=True and gene_content='symbol', keeps rows whose GENE field matches the mapped symbol."""
    df = _make_annotation_df(extra={"GENE": ["BRCA1", "TP53"]})
    gene_mapping = {"ENSG00000001": "BRCA1"}
    result = format_gene_annotation(
        "ENSG00000001", [df], True, gene_mapping, False, "GENE", "symbol"
    )

    assert len(result) == 1
    assert result.iloc[0]["GENE"] == "BRCA1"


@pytest.mark.unit
def test_format_gene_annotation_returns_empty_df_when_no_blocks():
    """Returns an empty DataFrame when the block list is empty."""
    result = format_gene_annotation("ENSG00000001", [], False, {}, False, "", "")

    assert result.empty
