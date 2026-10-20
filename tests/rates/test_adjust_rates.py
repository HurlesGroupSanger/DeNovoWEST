"""Tests for denovowest.rates.adjust_rates."""

import json
from pathlib import Path

import pandas as pd
import pytest
from click import UsageError
from click.testing import CliRunner

from denovowest.rates.adjust_rates import (
    apply_sperm_selection_adjustment,
    _load_adjustment_json,
    main,
)
from denovowest.utils.helpers import ensembl_ids_map_version

#############################
# Paths to fixtures         #
#############################

_ROOT = Path(__file__).parent.parent.parent

INPUT = _ROOT / "examples" / "rates" / "adjust_rates" / "input"
OUTPUT = _ROOT / "examples" / "rates" / "adjust_rates" / "output"

RATES = INPUT / "test_rates.tsv"
JSON = INPUT / "sperm_adjustment.json"
EXPECTED = OUTPUT / "mutation_rates_adjusted.tsv"

# Gene present in the JSON (unversioned) and in the rates file (versioned)
GENE_IN_JSON = "ENSG00000117139"
GENE_IN_RATES = "ENSG00000117139.19"

# Gene in the rates file but absent from the JSON (control — should be untouched)
GENE_CONTROL = "ENSG00000179295.20"

# Site that is present in both the JSON and the rates file
SITE_ID = "chr12_112450389_A_G"
SITE_INFO = {"chrom": "chr12", "pos": 112450389, "ref": "A", "alt": "G", "multiplier": 486.82}

#############################
# Unit tests — JSON loading #
#############################


@pytest.mark.unit
def test_load_adjustment_json_returns_genes_and_sites():
    genes, sites = _load_adjustment_json(JSON)

    assert GENE_IN_JSON in genes
    assert SITE_ID in sites


@pytest.mark.unit
def test_load_adjustment_json_genes_only(tmp_path):
    path = tmp_path / "adj.json"
    path.write_text(json.dumps({"genes": {"ENSG00000117139": {"missense": 2.0}}}))

    genes, sites = _load_adjustment_json(path)

    assert "ENSG00000117139" in genes
    assert sites == {}


@pytest.mark.unit
def test_load_adjustment_json_sites_only(tmp_path):
    path = tmp_path / "adj.json"
    path.write_text(
        json.dumps(
            {
                "sites": {
                    "chr1_100_A_G": {
                        "chrom": "chr1",
                        "pos": 100,
                        "ref": "A",
                        "alt": "G",
                        "multiplier": 5.0,
                    }
                }
            }
        )
    )

    genes, sites = _load_adjustment_json(path)

    assert genes == {}
    assert "chr1_100_A_G" in sites


@pytest.mark.unit
def test_load_adjustment_json_missing_keys_raises(tmp_path):
    path = tmp_path / "adj.json"
    path.write_text(json.dumps({"other_key": {}}))

    with pytest.raises(UsageError, match="at least one of the keys"):
        _load_adjustment_json(path)


@pytest.mark.unit
def test_load_adjustment_json_site_missing_field_raises(tmp_path):
    path = tmp_path / "adj.json"
    path.write_text(
        json.dumps(
            {
                "sites": {
                    "bad_site": {
                        "chrom": "chr1",
                        "pos": 100,
                        # ref, alt and multiplier are missing
                    }
                }
            }
        )
    )

    with pytest.raises(UsageError, match="missing required fields"):
        _load_adjustment_json(path)


######################################
# Unit tests — version ID mapping    #
######################################


def _rates_gene_ids() -> list[str]:
    return pd.read_csv(RATES, sep="\t", usecols=["gene_id"]).gene_id.unique().tolist()


@pytest.mark.unit
def test_version_mapping_resolves_unversioned_to_versioned():
    """JSON has unversioned ENSG00000117139; rates file has ENSG00000117139.19."""
    genes, _ = _load_adjustment_json(JSON)
    mapping = ensembl_ids_map_version(list(genes), _rates_gene_ids())

    assert mapping[GENE_IN_JSON] == GENE_IN_RATES


@pytest.mark.unit
def test_version_mapping_absent_gene_not_in_result():
    """Genes that have no match in the rates file must be absent from the mapping."""
    genes, _ = _load_adjustment_json(JSON)
    mapping = ensembl_ids_map_version(list(genes), _rates_gene_ids())

    # All JSON genes except ENSG00000117139 are absent from the tiny rates fixture
    absent_genes = [g for g in genes if g != GENE_IN_JSON]
    for gene in absent_genes:
        assert gene not in mapping


@pytest.mark.unit
def test_version_mapping_handles_versioned_and_mismatched_versions():
    mapping = ensembl_ids_map_version(
        ["ENSG00000117139.19", "ENSG00000117139.5", "ENSG00000179295"],
        ["ENSG00000117139", "ENSG00000179295.20"],
    )

    assert mapping == {
        "ENSG00000117139.19": "ENSG00000117139",
        "ENSG00000117139.5": "ENSG00000117139",
        "ENSG00000179295": "ENSG00000179295.20",
    }


#############################################
# Unit tests — applying adjustments         #
#############################################


def _make_chunk(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(rows)


@pytest.mark.unit
def test_apply_gene_missense_multiplier():
    chunk = _make_chunk(
        [
            {
                "gene_id": "ENSG00000117139.19",
                "consequence": "missense",
                "prob": 1.0,
                "chrom": "chr1",
                "pos": 999,
                "ref": "A",
                "alt": "G",
            }
        ]
    )
    genes = {GENE_IN_JSON: {"missense": 4.88, "stop_gained": 27.92}}
    mapping = {GENE_IN_JSON: GENE_IN_RATES}

    result, _ = apply_sperm_selection_adjustment(chunk, genes, {}, mapping)

    assert pytest.approx(result.iloc[0]["prob"]) == 4.88


@pytest.mark.unit
def test_apply_gene_stop_gained_multiplier():
    chunk = _make_chunk(
        [
            {
                "gene_id": "ENSG00000117139.19",
                "consequence": "stop_gained",
                "prob": 1.0,
                "chrom": "chr1",
                "pos": 999,
                "ref": "A",
                "alt": "G",
            }
        ]
    )
    genes = {GENE_IN_JSON: {"missense": 4.88, "stop_gained": 27.92}}
    mapping = {GENE_IN_JSON: GENE_IN_RATES}

    result, _ = apply_sperm_selection_adjustment(chunk, genes, {}, mapping)

    assert pytest.approx(result.iloc[0]["prob"]) == 27.92


@pytest.mark.unit
def test_apply_leaves_non_target_consequence_unchanged():
    chunk = _make_chunk(
        [
            {
                "gene_id": "ENSG00000117139.19",
                "consequence": "3_prime_utr",
                "prob": 1.0,
                "chrom": "chr1",
                "pos": 999,
                "ref": "A",
                "alt": "G",
            }
        ]
    )
    genes = {GENE_IN_JSON: {"missense": 4.88}}
    mapping = {GENE_IN_JSON: GENE_IN_RATES}

    result, _ = apply_sperm_selection_adjustment(chunk, genes, {}, mapping)

    assert result.iloc[0]["prob"] == 1.0


@pytest.mark.unit
def test_apply_leaves_unmapped_gene_unchanged():
    """Gene present in the adjustment dict but absent from the mapping must be skipped."""
    chunk = _make_chunk(
        [
            {
                "gene_id": "ENSG00000117139.19",
                "consequence": "missense",
                "prob": 1.0,
                "chrom": "chr1",
                "pos": 999,
                "ref": "A",
                "alt": "G",
            }
        ]
    )
    genes = {GENE_IN_JSON: {"missense": 4.88}}
    # Empty mapping — gene_id is not resolved
    result, _ = apply_sperm_selection_adjustment(chunk, genes, {}, mapping_gene_ids={})

    assert result.iloc[0]["prob"] == 1.0


@pytest.mark.unit
def test_apply_site_multiplier():
    chunk = _make_chunk(
        [
            {
                "gene_id": GENE_CONTROL,
                "consequence": "missense",
                "prob": 1.0,
                "chrom": "chr12",
                "pos": 112450389,
                "ref": "A",
                "alt": "G",
            }
        ]
    )
    result, adjusted_genes = apply_sperm_selection_adjustment(chunk, {}, {SITE_ID: SITE_INFO}, {})

    assert pytest.approx(result.iloc[0]["prob"]) == 486.82
    assert adjusted_genes == [GENE_CONTROL]


@pytest.mark.unit
def test_apply_unmatched_site_adjusts_no_gene():
    chunk = _make_chunk(
        [
            {
                "gene_id": GENE_CONTROL,
                "consequence": "missense",
                "prob": 1.0,
                "chrom": "chr1",
                "pos": 1,
                "ref": "A",
                "alt": "G",
            }
        ]
    )
    result, adjusted_genes = apply_sperm_selection_adjustment(chunk, {}, {SITE_ID: SITE_INFO}, {})

    assert result.iloc[0]["prob"] == 1.0
    assert adjusted_genes == []


######################################
# Integration test — CLI round-trip  #
######################################


@pytest.mark.smoke
def test_cli_output_matches_expected(tmp_path):
    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            "--rates", str(RATES),
            "--sperm-selection-adjustment", str(JSON),
            "--outdir", str(tmp_path),
            "--output", "out.tsv",
        ],
    )

    assert result.exit_code == 0, result.output

    actual = pd.read_csv(tmp_path / "out.tsv", sep="\t")
    expected = pd.read_csv(EXPECTED, sep="\t")

    pd.testing.assert_frame_equal(actual.reset_index(drop=True), expected.reset_index(drop=True))
