"""Tests for denovowest.dnm.filter_dnm."""

from pathlib import Path

import gffutils
import pandas as pd
import pytest
from click.testing import CliRunner

from denovowest.dnm.filter_dnm import (
    filter_dnm,
    filter_on_gene_list,
    filter_on_gff,
    format_gene_id,
    load_gene_list,
)
from denovowest.utils.params import CDS_OFFSET

#############################
# Paths to fixtures         #
#############################

_ROOT = Path(__file__).parent.parent.parent
_EXAMPLES = _ROOT / "examples" / "dnm" / "filter_dnm"

DNM_INPUT = _EXAMPLES / "input" / "DNM_b38_MANE.tsv"
GENE_LIST_INPUT = _EXAMPLES / "input" / "gene_list.tsv"
KEPT_REF = _EXAMPLES / "output" / "kept_dnm.tsv"
DISCARDED_REF = _EXAMPLES / "output" / "discarded_dnm.tsv"


#############################
# Unit – format_gene_id     #
#############################


@pytest.mark.unit
def test_format_gene_id_unversioned_dnm_matched_to_versioned_list():
    """Core use case: DNM has bare IDs, rates file has versioned ones."""
    result = format_gene_id(
        genes_in_dnm=["ENSG00000012048", "ENSG00000139618"],
        genes_list=["ENSG00000012048.19", "ENSG00000139618.7"],
    )
    assert result == ["ENSG00000012048.19", "ENSG00000139618.7"]


@pytest.mark.unit
def test_format_gene_id_both_versioned_no_change():
    result = format_gene_id(
        genes_in_dnm=["ENSG00000012048.19"],
        genes_list=["ENSG00000012048.19"],
    )
    assert result == ["ENSG00000012048.19"]


@pytest.mark.unit
def test_format_gene_id_both_unversioned_no_change():
    result = format_gene_id(
        genes_in_dnm=["ENSG00000012048"],
        genes_list=["ENSG00000012048"],
    )
    assert result == ["ENSG00000012048"]


@pytest.mark.unit
def test_format_gene_id_unversioned_not_in_list_kept_as_is():
    """Unversioned ID absent from the rates list stays unchanged."""
    result = format_gene_id(
        genes_in_dnm=["ENSG00000012048", "ENSG00000000001"],
        genes_list=["ENSG00000012048.19"],
    )
    assert result == ["ENSG00000012048.19", "ENSG00000000001"]


@pytest.mark.unit
def test_format_gene_id_strips_whitespace():
    result = format_gene_id(
        genes_in_dnm=["  ENSG00000012048  "],
        genes_list=["ENSG00000012048.19"],
    )
    assert result == ["ENSG00000012048.19"]


@pytest.mark.unit
def test_format_gene_id_empty_inputs():
    assert format_gene_id([], []) == []
    assert format_gene_id([], ["ENSG00000012048.19"]) == []


@pytest.mark.unit
def test_format_gene_id_non_ensembl_ids_unchanged():
    """Non-ENSEMBL gene identifiers (e.g. gene symbols) are passed through."""
    result = format_gene_id(
        genes_in_dnm=["BRCA1", "BRCA2"],
        genes_list=["BRCA1", "BRCA2"],
    )
    assert result == ["BRCA1", "BRCA2"]


#############################
# Unit – filter_on_gene_list
#############################


@pytest.mark.unit
def test_filter_on_gene_list_keeps_matching_genes():
    dnm = pd.DataFrame({
        "gene_id": ["G1", "G1", "G2", "G3"],
        "chrom": ["chr1"] * 4,
        "pos": [100, 200, 300, 400],
        "ref": ["A"] * 4,
        "alt": ["T"] * 4,
    })
    kept, discarded = filter_on_gene_list(dnm, ["G1", "G2"])

    assert set(kept["gene_id"]) == {"G1", "G2"}
    assert len(kept) == 3


@pytest.mark.unit
def test_filter_on_gene_list_discards_with_correct_reason():
    dnm = pd.DataFrame({
        "gene_id": ["G1", "G2"],
        "chrom": ["chr1", "chr1"],
        "pos": [100, 200],
        "ref": ["A", "A"],
        "alt": ["T", "T"],
    })
    _, discarded = filter_on_gene_list(dnm, ["G1"])

    assert list(discarded["gene_id"]) == ["G2"]
    assert list(discarded["reason"]) == ["not_in_gene_list"]


@pytest.mark.unit
def test_filter_on_gene_list_empty_gene_list_discards_all():
    dnm = pd.DataFrame({
        "gene_id": ["G1", "G2"],
        "chrom": ["chr1", "chr1"],
        "pos": [100, 200],
        "ref": ["A", "A"],
        "alt": ["T", "T"],
    })
    kept, discarded = filter_on_gene_list(dnm, [])

    assert len(kept) == 0
    assert len(discarded) == 2


@pytest.mark.unit
def test_filter_on_gene_list_returns_independent_copies():
    """Modifying the kept slice must not affect the discarded slice."""
    dnm = pd.DataFrame({
        "gene_id": ["G1", "G2"],
        "chrom": ["chr1", "chr1"],
        "pos": [100, 200],
        "ref": ["A", "A"],
        "alt": ["T", "T"],
    })
    kept, discarded = filter_on_gene_list(dnm, ["G1"])
    kept["new_col"] = "x"

    assert "new_col" not in discarded.columns


#############################
# Unit – filter_on_gff      #
#############################

# Gene structure used by the tests below:
#   gene   chr1:1000-2000
#   CDS    chr1:1200-1800   (effective range with CDS_OFFSET: 1150-1850)
#
# Variants placed at three positions:
#   pos 1100 → outside (< 1150)  → discarded
#   pos 1500 → inside             → kept
#   pos 1900 → outside (> 1850)  → discarded

_MINIMAL_GFF = """\
##gff-version 3
chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=ENSG00000000001
chr1\t.\ttranscript\t1000\t2000\t.\t+\t.\tID=ENST00000000001;Parent=ENSG00000000001
chr1\t.\tCDS\t1200\t1800\t.\t+\t0\tID=CDS0001;Parent=ENST00000000001
"""

_DNM_AROUND_CDS = pd.DataFrame({
    "gene_id": ["ENSG00000000001"] * 3,
    "chrom":   ["chr1"] * 3,
    "pos":     [1100, 1500, 1900],
    "ref":     ["A"] * 3,
    "alt":     ["T"] * 3,
})


@pytest.fixture(scope="module")
def minimal_gff_path(tmp_path_factory):
    """Write the minimal GFF to a temp file and return its path."""
    gff = tmp_path_factory.mktemp("gff") / "test.gff"
    gff.write_text(_MINIMAL_GFF)
    return gff


@pytest.mark.unit
def test_filter_on_gff_keeps_only_variant_inside_cds(minimal_gff_path):
    kept, discarded = filter_on_gff(_DNM_AROUND_CDS.copy(), minimal_gff_path)

    assert len(kept) == 1
    assert int(kept.iloc[0]["pos"]) == 1500


@pytest.mark.unit
def test_filter_on_gff_discards_variants_outside_cds(minimal_gff_path):
    kept, discarded = filter_on_gff(_DNM_AROUND_CDS.copy(), minimal_gff_path)

    assert len(discarded) == 2
    assert set(discarded["pos"]) == {1100, 1900}
    assert (discarded["reason"] == "not_in_cds").all()


@pytest.mark.unit
def test_filter_on_gff_cds_offset_is_respected(minimal_gff_path):
    """Variants just within the CDS_OFFSET extension must be retained."""
    # pos 1200 - CDS_OFFSET = 1200 - 50 = 1150  →  pos 1155 is inside the padded region
    just_inside = _DNM_AROUND_CDS.copy()
    just_inside = pd.DataFrame({
        "gene_id": ["ENSG00000000001"] * 2,
        "chrom":   ["chr1"] * 2,
        "pos":     [1200 - CDS_OFFSET, 1800 + CDS_OFFSET],  # exact boundary → kept
        "ref":     ["A"] * 2,
        "alt":     ["T"] * 2,
    })
    kept, discarded = filter_on_gff(just_inside, minimal_gff_path)

    assert len(kept) == 2
    assert len(discarded) == 0


@pytest.mark.unit
def test_filter_on_gff_gene_absent_from_gff_all_discarded(minimal_gff_path):
    dnm = pd.DataFrame({
        "gene_id": ["ENSG99999999999"],
        "chrom":   ["chr1"],
        "pos":     [1500],
        "ref":     ["A"],
        "alt":     ["T"],
    })
    kept, discarded = filter_on_gff(dnm, minimal_gff_path)

    assert len(kept) == 0
    assert list(discarded["reason"]) == ["gene_not_in_gff"]


#############################
# Unit – load_gene_list     #
#############################


@pytest.mark.unit
def test_load_gene_list_reads_identifiers(tmp_path):
    f = tmp_path / "genes.txt"
    f.write_text("ENSG00000012048.19\nENSG00000139618.7\n")

    result = load_gene_list(f)

    assert result == ["ENSG00000012048.19", "ENSG00000139618.7"]


@pytest.mark.unit
def test_load_gene_list_skips_blank_lines(tmp_path):
    f = tmp_path / "genes.txt"
    f.write_text("ENSG00000012048.19\n\nENSG00000139618.7\n\n")

    result = load_gene_list(f)

    assert result == ["ENSG00000012048.19", "ENSG00000139618.7"]


#############################
# CLI parameter validation  #
#############################


@pytest.fixture
def runner():
    return CliRunner()


@pytest.mark.unit
def test_missing_dnm_argument_exits_nonzero(runner):
    result = runner.invoke(filter_dnm, [])
    assert result.exit_code != 0


@pytest.mark.unit
def test_nonexistent_dnm_exits_nonzero(runner, tmp_path):
    result = runner.invoke(filter_dnm, [
        str(tmp_path / "missing.tsv"),
        str(GENE_LIST_INPUT),
        "--output_kept_dnm", str(tmp_path / "kept.tsv"),
        "--output_discarded_dnm", str(tmp_path / "disc.tsv"),
    ])
    assert result.exit_code != 0


#############################
# Functional / non-regression
#############################


@pytest.fixture(scope="module")
def filter_results(tmp_path_factory):
    """Run filter_dnm once on the example data and return (kept_df, discarded_df)."""
    out = tmp_path_factory.mktemp("filter_dnm")
    kept_path = out / "kept.tsv"
    disc_path = out / "discarded.tsv"

    result = CliRunner().invoke(filter_dnm, [
        str(DNM_INPUT),
        str(GENE_LIST_INPUT),
        "--output_kept_dnm", str(kept_path),
        "--output_discarded_dnm", str(disc_path),
    ])
    assert result.exit_code == 0, result.output

    return pd.read_csv(kept_path, sep="\t"), pd.read_csv(disc_path, sep="\t")


@pytest.mark.functional
def test_functional_kept_count(filter_results):
    kept, _ = filter_results
    reference = pd.read_csv(KEPT_REF, sep="\t")
    assert len(kept) == len(reference)


@pytest.mark.functional
def test_functional_discarded_count(filter_results):
    _, discarded = filter_results
    reference = pd.read_csv(DISCARDED_REF, sep="\t")
    assert len(discarded) == len(reference)


@pytest.mark.functional
def test_functional_discarded_reason_column(filter_results):
    _, discarded = filter_results
    assert "reason" in discarded.columns
    assert discarded["reason"].notna().all()


@pytest.mark.functional
def test_functional_no_overlap_between_kept_and_discarded(filter_results):
    kept, discarded = filter_results
    kept_ids = set(zip(kept.chrom, kept.pos, kept.ref, kept.alt))
    disc_ids = set(zip(discarded.chrom, discarded.pos, discarded.ref, discarded.alt))
    assert kept_ids.isdisjoint(disc_ids)


@pytest.mark.functional
def test_non_regression_kept(filter_results, tmp_path):
    kept, _ = filter_results
    out = tmp_path / "kept.tsv"
    kept.to_csv(out, sep="\t", index=False)
    assert out.read_bytes() == KEPT_REF.read_bytes()


@pytest.mark.functional
def test_non_regression_discarded(filter_results, tmp_path):
    _, discarded = filter_results
    out = tmp_path / "discarded.tsv"
    discarded.to_csv(out, sep="\t", index=False)
    assert out.read_bytes() == DISCARDED_REF.read_bytes()
