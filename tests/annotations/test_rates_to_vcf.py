"""Tests for denovowest.annotations.rates_to_vcf."""

from pathlib import Path

import gffutils
import pandas as pd
import pytest
from click.testing import CliRunner

from denovowest.annotations.rates_to_vcf import assign_indel_csq, cli, is_indel

#############################
# Paths to fixtures         #
#############################

_ROOT = Path(__file__).parent.parent.parent
_GFF = _ROOT / "examples" / "annotation" / "CADD" / "input" / "OR4F5.gff"


@pytest.fixture(scope="module")
def or4f5_db(tmp_path_factory):
    """Build a gffutils database for OR4F5 once per test session."""
    db_path = tmp_path_factory.mktemp("gff") / "OR4F5.db"
    gffutils.create_db(str(_GFF), str(db_path), merge_strategy="create_unique")
    return str(db_path)


#########################################
# Helpers                               #
#########################################


def _make_rates_file(tmp_path: Path, extra_columns: dict | None = None) -> Path:
    """Write a minimal rates TSV for one OR4F5 variant."""
    data = {
        "gene_id": ["ENSG00000186092.7"],
        "chrom": ["1"],
        "pos": [65569],
        "ref": ["A"],
        "alt": ["T"],
        "prob": [0.001],
    }
    if extra_columns:
        data.update(extra_columns)
    path = tmp_path / "rates.tsv"
    pd.DataFrame(data).to_csv(path, sep="\t", index=False)
    return path


def _make_bcftools_vcf(tmp_path: Path, consequence: str = "missense") -> Path:
    """Write a minimal bcftools csq-style VCF for one OR4F5 variant."""
    bcsq = f"{consequence}|OR4F5|ENST00000641515|protein_coding|+|K2M|65569A>T"
    content = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"1\t65569\t.\tA\tT\t.\t.\tGENE=ENSG00000186092.7;BCSQ={bcsq}\n"
    )
    path = tmp_path / "annotated.vcf"
    path.write_text(content)
    return path


#############################
# Unit – is_indel           #
#############################


@pytest.mark.unit
def test_is_indel_returns_false_for_snp():
    assert not is_indel(pd.Series({"REF": "A", "ALT": "T"}))


@pytest.mark.unit
def test_is_indel_returns_true_for_insertion():
    assert is_indel(pd.Series({"REF": "A", "ALT": "AT"}))


@pytest.mark.unit
def test_is_indel_returns_true_for_deletion():
    assert is_indel(pd.Series({"REF": "AT", "ALT": "A"}))


#############################
# Unit – assign_indel_csq   #
#############################


@pytest.mark.unit
def test_assign_indel_csq_keeps_frameshift_unchanged():
    """Consequences already classified as indel are returned as-is."""
    csq, csqs = assign_indel_csq("frameshift", "frameshift")
    assert csq == "frameshift"


@pytest.mark.unit
def test_assign_indel_csq_keeps_inframe_unchanged():
    csq, csqs = assign_indel_csq("inframe", "inframe")
    assert csq == "inframe"


@pytest.mark.unit
def test_assign_indel_csq_coerces_to_inframe_when_inframe_in_csqs():
    """When csqs contains an inframe term, consequence is coerced to it."""
    csq, csqs = assign_indel_csq("missense", "missense&inframe_insertion")
    assert csq == "inframe_insertion"


@pytest.mark.unit
def test_assign_indel_csq_coerces_to_frameshift_when_frameshift_in_csqs():
    csq, csqs = assign_indel_csq("missense", "missense&frameshift")
    assert csq == "frameshift"


#############################
# Smoke tests               #
#############################


@pytest.mark.smoke
def test_vcf_to_rates_smoke(tmp_path, or4f5_db):
    """CLI exits cleanly and adds a 'consequence' column to the rates file."""
    rates_file = _make_rates_file(tmp_path)
    vcf_file = _make_bcftools_vcf(tmp_path)
    output = tmp_path / "out.tsv"
    runner = CliRunner()

    result = runner.invoke(
        cli, ["vcf-to-rates", str(vcf_file), str(rates_file), or4f5_db, str(output)]
    )

    assert result.exit_code == 0, result.output
    out_df = pd.read_csv(output, sep="\t")
    assert "consequence" in out_df.columns
    assert out_df["consequence"].iloc[0] == "missense"


#############################
# Conflict handling         #
#############################


@pytest.mark.unit
def test_vcf_to_rates_renames_existing_consequence_column(tmp_path, or4f5_db):
    """When the rates file already has 'consequence', it is renamed to 'consequence_prev'
    and the new bcftools csq annotation takes the canonical 'consequence' name."""
    rates_file = _make_rates_file(tmp_path, extra_columns={"consequence": ["synonymous"]})
    vcf_file = _make_bcftools_vcf(tmp_path, consequence="missense")
    output = tmp_path / "out.tsv"
    runner = CliRunner()

    result = runner.invoke(
        cli, ["vcf-to-rates", str(vcf_file), str(rates_file), or4f5_db, str(output)]
    )

    assert result.exit_code == 0, result.output
    # Warning names the renamed column
    assert "consequence_prev" in result.output

    out_df = pd.read_csv(output, sep="\t")
    # New bcftools annotation occupies the canonical column name
    assert "consequence" in out_df.columns
    assert out_df["consequence"].iloc[0] == "missense"
    # Original column is preserved under its renamed form
    assert "consequence_prev" in out_df.columns
    assert out_df["consequence_prev"].iloc[0] == "synonymous"
