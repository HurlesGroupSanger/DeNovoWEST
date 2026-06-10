"""Tests for denovowest.rates.get_rates_per_gene_per_consequence."""

from pathlib import Path

import pandas as pd
import pytest
from click.testing import CliRunner

from denovowest.rates.get_rates_per_gene_per_consequence import (
    CONSEQUENCES,
    _iter_rates,
    _write_matrix,
    compute_rates_per_gene_per_consequence,
    main,
)

###########
# Helpers #
###########


def _write_rates_tsv(path: Path, rows: list[dict]) -> None:
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


#############################
# Unit – _iter_rates        #
#############################


@pytest.mark.unit
def test_iter_rates_yields_one_group_per_gene(tmp_path):
    path = tmp_path / "rates.tsv"
    _write_rates_tsv(path, [
        {"gene_id": "G1", "consequence": "synonymous",  "prob": 0.1},
        {"gene_id": "G2", "consequence": "missense",    "prob": 0.2},
    ])

    gene_ids = [gid for gid, _ in _iter_rates(path)]

    assert gene_ids == ["G1", "G2"]


@pytest.mark.unit
def test_iter_rates_group_contains_correct_rows(tmp_path):
    path = tmp_path / "rates.tsv"
    _write_rates_tsv(path, [
        {"gene_id": "G1", "consequence": "synonymous",  "prob": 0.1},
        {"gene_id": "G1", "consequence": "missense",    "prob": 0.2},
        {"gene_id": "G2", "consequence": "stop_gained", "prob": 0.3},
    ])

    groups = dict(_iter_rates(path))

    assert len(groups["G1"]) == 2
    assert len(groups["G2"]) == 1


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


############################
# Unit – _write_matrix     #
############################


@pytest.mark.unit
def test_write_matrix_creates_file(tmp_path):
    output = tmp_path / "matrix.tsv"

    _write_matrix({"G1": {"synonymous": 0.1}}, output)

    assert output.exists()


@pytest.mark.unit
def test_write_matrix_is_readable_as_dataframe(tmp_path):
    output = tmp_path / "matrix.tsv"
    data = {"G1": {"synonymous": 0.1, "missense": 0.2}, "G2": {"synonymous": 0.3, "missense": 0.4}}
    _write_matrix(data, output)

    df = pd.read_csv(output, sep="\t", index_col=0)

    assert df.shape == (2, 2)
    assert set(df.index) == {"G1", "G2"}
    assert df.loc["G1", "synonymous"] == pytest.approx(0.1)


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


################################
# Functional – CLI end-to-end  #
################################


@pytest.fixture(scope="module")
def rates_cli_output(tmp_path_factory):
    """Run the CLI on minimal synthetic data and return (prob_df, count_df)."""
    tmp = tmp_path_factory.mktemp("rates_cli")
    rates_path = tmp / "rates.tsv"
    prob_path   = tmp / "prob.tsv"
    count_path  = tmp / "count.tsv"

    pd.DataFrame([
        {"gene_id": "G1", "consequence": "synonymous",  "prob": 0.1},
        {"gene_id": "G1", "consequence": "synonymous",  "prob": 0.2},
        {"gene_id": "G1", "consequence": "missense",    "prob": 0.05},
        {"gene_id": "G2", "consequence": "stop_gained", "prob": 0.3},
    ]).to_csv(rates_path, sep="\t", index=False)

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


@pytest.mark.functional
class TestRatesCLIOutput:

    def test_prob_matrix_has_all_consequence_columns(self, rates_cli_output):
        prob_df, _ = rates_cli_output
        assert set(CONSEQUENCES).issubset(prob_df.columns)

    def test_prob_matrix_has_expected_genes(self, rates_cli_output):
        prob_df, _ = rates_cli_output
        assert set(prob_df.index) == {"G1", "G2"}

    def test_prob_synonymous_sum_is_correct(self, rates_cli_output):
        prob_df, _ = rates_cli_output
        assert prob_df.loc["G1", "synonymous"] == pytest.approx(0.3)

    def test_prob_missense_is_correct(self, rates_cli_output):
        prob_df, _ = rates_cli_output
        assert prob_df.loc["G1", "missense"] == pytest.approx(0.05)

    def test_prob_stop_gained_is_correct(self, rates_cli_output):
        prob_df, _ = rates_cli_output
        assert prob_df.loc["G2", "stop_gained"] == pytest.approx(0.3)

    def test_count_synonymous_is_correct(self, rates_cli_output):
        _, count_df = rates_cli_output
        assert count_df.loc["G1", "synonymous"] == 2

    def test_count_missense_is_correct(self, rates_cli_output):
        _, count_df = rates_cli_output
        assert count_df.loc["G1", "missense"] == 1

    def test_count_stop_gained_is_correct(self, rates_cli_output):
        _, count_df = rates_cli_output
        assert count_df.loc["G2", "stop_gained"] == 1
