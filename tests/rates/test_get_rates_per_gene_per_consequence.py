"""Tests for denovowest.rates.get_rates_per_gene_per_consequence."""

from pathlib import Path

import pandas as pd
import pytest
from click.testing import CliRunner

from denovowest.rates.get_rates_per_gene_per_consequence import (
    CONSEQUENCES,
    _iter_rates,
    compute_rates_per_gene_per_consequence,
    main,
)

#####################
# Regression data   #
#####################

# 10 synthetic genes: gene_id → size multiplier (relative to base prob sums)
_GENE_MULTIPLIERS: dict[str, float] = {
    "ENSG00000000001": 1.00,
    "ENSG00000000002": 1.20,
    "ENSG00000000003": 0.80,
    "ENSG00000000004": 1.50,
    "ENSG00000000005": 1.00,
    "ENSG00000000006": 1.10,
    "ENSG00000000007": 0.90,
    "ENSG00000000008": 1.30,
    "ENSG00000000009": 0.70,
    "ENSG00000000010": 1.00,
}

# Sum of mutation probabilities across all sites in a gene of size ×1.0
_BASE_PROB_SUM: dict[str, float] = {
    "synonymous":      5.00e-3,
    "missense":        1.00e-2,
    "stop_gained":     5.00e-4,
    "start_lost":      5.00e-5,
    "stop_lost":       5.00e-5,
    "splice_region":   5.00e-4,
    "splice_acceptor": 1.00e-4,
    "splice_donor":    1.00e-4,
}

_N_SITES = 10  # per-site rows per (gene, consequence)

###########
# Helpers #
###########


def _write_rates_tsv(path: Path, rows: list[dict]) -> None:
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def _make_full_rates_tsv(path: Path) -> None:
    """Write a rates TSV for all 10 regression genes (3 columns only: gene_id, consequence, prob)."""
    rows = [
        {"gene_id": gene_id, "consequence": csq,
         "prob": _BASE_PROB_SUM[csq] * mult / _N_SITES}
        for gene_id, mult in _GENE_MULTIPLIERS.items()
        for csq in CONSEQUENCES
        for _ in range(_N_SITES)
    ]
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


#############################
# Unit – _iter_rates        #
#############################


@pytest.mark.unit
def test_iter_rates_empty_file_yields_nothing(tmp_path):
    path = tmp_path / "rates.tsv"
    path.write_text("gene_id\tconsequence\tprob\n")

    assert list(_iter_rates(path)) == []


###################################################
# Unit – compute_rates_per_gene_per_consequence   #
###################################################


@pytest.mark.unit
def test_compute_rates_returns_correct_prob_sum(tmp_path):
    path = tmp_path / "rates.tsv"
    _write_rates_tsv(path, [
        {"gene_id": "G1", "consequence": "synonymous", "prob": 0.1},
        {"gene_id": "G1", "consequence": "synonymous", "prob": 0.2},
        {"gene_id": "G2", "consequence": "missense",   "prob": 0.5},
    ])

    prob, _ = compute_rates_per_gene_per_consequence(path)

    assert prob["G1"]["synonymous"] == pytest.approx(0.3)
    assert prob["G2"]["missense"]   == pytest.approx(0.5)


@pytest.mark.unit
def test_compute_rates_absent_consequence_is_zero(tmp_path):
    path = tmp_path / "rates.tsv"
    _write_rates_tsv(path, [
        {"gene_id": "G1", "consequence": "synonymous", "prob": 0.1},
    ])

    prob, _ = compute_rates_per_gene_per_consequence(path)

    assert prob["G1"]["missense"] == 0.0


@pytest.mark.unit
def test_compute_rates_count_matches_number_of_records(tmp_path):
    path = tmp_path / "rates.tsv"
    _write_rates_tsv(path, [
        {"gene_id": "G1", "consequence": "synonymous", "prob": 0.1},
        {"gene_id": "G1", "consequence": "synonymous", "prob": 0.2},
        {"gene_id": "G1", "consequence": "missense",   "prob": 0.3},
    ])

    _, count = compute_rates_per_gene_per_consequence(path)

    assert count["G1"]["synonymous"] == 2
    assert count["G1"]["missense"]   == 1


@pytest.mark.unit
def test_compute_rates_accumulates_non_adjacent_rows_of_same_gene(tmp_path):
    """Rows for the same gene that are not adjacent (cross-chunk boundary) must be summed."""
    path = tmp_path / "rates.tsv"
    _write_rates_tsv(path, [
        {"gene_id": "G1", "consequence": "synonymous", "prob": 0.1},
        {"gene_id": "G2", "consequence": "synonymous", "prob": 0.5},
        {"gene_id": "G1", "consequence": "synonymous", "prob": 0.2},  # G1 reappears
    ])

    prob, count = compute_rates_per_gene_per_consequence(path)

    assert prob["G1"]["synonymous"]  == pytest.approx(0.3)
    assert count["G1"]["synonymous"] == 2


@pytest.mark.unit
def test_compute_rates_all_consequences_present_in_output(tmp_path):
    path = tmp_path / "rates.tsv"
    _write_rates_tsv(path, [
        {"gene_id": "G1", "consequence": "synonymous", "prob": 0.1},
    ])

    prob, _ = compute_rates_per_gene_per_consequence(path)

    assert set(prob["G1"].keys()) == set(CONSEQUENCES)


##############################
# CLI – parameter validation #
##############################


@pytest.fixture
def runner():
    return CliRunner()


@pytest.mark.unit
def test_missing_rates_argument_exits_nonzero(runner):
    result = runner.invoke(main, [])

    assert result.exit_code != 0


@pytest.mark.unit
def test_nonexistent_rates_file_exits_nonzero(runner, tmp_path):
    result = runner.invoke(main, [str(tmp_path / "does_not_exist.tsv")])

    assert result.exit_code != 0


##############
# Regression #
##############


@pytest.fixture(scope="module")
def regression_rates_output(tmp_path_factory):
    """Run the CLI on the 10-gene synthetic dataset generated on-the-fly."""
    tmp = tmp_path_factory.mktemp("regression_rates")
    rates_path = tmp / "mutation_rates.tsv"
    prob_path  = tmp / "prob.tsv"
    count_path = tmp / "count.tsv"

    _make_full_rates_tsv(rates_path)

    result = CliRunner().invoke(main, [
        str(rates_path),
        "--output",    str(prob_path),
        "--output-nb", str(count_path),
    ])
    assert result.exit_code == 0, result.output

    return (
        pd.read_csv(prob_path,  sep="\t", index_col=0),
        pd.read_csv(count_path, sep="\t", index_col=0),
    )


@pytest.mark.regression
class TestRatesNonRegression:
    """Verify that aggregated prob sums and record counts match analytical expectations.

    Expected values are derived from the test-data constants (_BASE_PROB_SUM,
    _GENE_MULTIPLIERS, _N_SITES), not from stored reference files.  If the
    aggregation logic changes, these tests will fail.
    """

    def test_prob_sum_equals_base_times_multiplier(self, regression_rates_output):
        prob_df, _ = regression_rates_output
        for gene_id, mult in _GENE_MULTIPLIERS.items():
            for csq in CONSEQUENCES:
                expected = _BASE_PROB_SUM[csq] * mult
                assert prob_df.loc[gene_id, csq] == pytest.approx(expected, rel=1e-9), \
                    f"Unexpected prob_sum for {gene_id}.{csq}"

    def test_record_count_equals_n_sites_for_every_cell(self, regression_rates_output):
        _, count_df = regression_rates_output
        for gene_id in _GENE_MULTIPLIERS:
            for csq in CONSEQUENCES:
                assert count_df.loc[gene_id, csq] == _N_SITES, \
                    f"Unexpected record count for {gene_id}.{csq}"

    def test_all_genes_present_in_output(self, regression_rates_output):
        prob_df, _ = regression_rates_output
        assert set(prob_df.index) == set(_GENE_MULTIPLIERS.keys())

