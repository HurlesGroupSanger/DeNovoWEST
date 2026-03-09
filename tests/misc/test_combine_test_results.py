import importlib.util
from pathlib import Path

import pandas as pd
import pytest
from click.testing import CliRunner

pytestmark = pytest.mark.misc


def _load_combine_test_results_command():
    script_path = (
        Path(__file__).resolve().parents[2] / "misc" / "denovonear" / "combine_test_results.py"
    )
    spec = importlib.util.spec_from_file_location("combine_test_results_script", script_path)
    module = importlib.util.module_from_spec(spec)
    assert spec is not None
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module.combine_test_results


def _run_combine_test_results(tmp_path, linear_df, clustering_3d_df):
    combine_test_results = _load_combine_test_results_command()

    linear_path = tmp_path / "linear.tsv"
    clustering_3d_path = tmp_path / "three_d.tsv"
    out_path = tmp_path / "combined.tsv"
    linear_df.to_csv(linear_path, sep="\t")
    clustering_3d_df.to_csv(clustering_3d_path, sep="\t")

    runner = CliRunner()
    result = runner.invoke(
        combine_test_results, [str(linear_path), str(clustering_3d_path), str(out_path)]
    )
    combined = pd.read_csv(out_path, sep="\t", index_col=0)
    return result, combined


def test_combine_test_results_prefers_3d_for_overlapping_missense_rows(tmp_path):
    linear_df = pd.DataFrame(
        [
            {
                "gene_id": "GENE1",
                "gene_symbol": "SYM1",
                "mutation_category": "missense",
                "events_n": 10,
                "dist": 0.11,
                "probability": 0.01,
            },
            {
                "gene_id": "GENE2",
                "gene_symbol": "SYM2",
                "mutation_category": "missense",
                "events_n": 20,
                "dist": 0.22,
                "probability": 0.02,
            },
        ]
    ).set_index("gene_id")

    clustering_3d_df = pd.DataFrame(
        [
            {
                "gene_id": "GENE1",
                "gene_symbol": "SYM1",
                "mutation_category": "missense",
                "events_n": 10,
                "dist": 1.11,
                "probability": 0.001,
            },
            {
                "gene_id": "GENE2",
                "gene_symbol": "SYM2",
                "mutation_category": "missense",
                "events_n": 20,
                "dist": 2.22,
                "probability": 0.002,
            },
        ]
    ).set_index("gene_id")

    result, combined = _run_combine_test_results(tmp_path, linear_df, clustering_3d_df)

    assert result.exit_code == 0
    assert "Number of genes with 3D clustering results: 2" in result.output
    assert "Number of genes with linear clustering results: 0" in result.output

    assert set(combined.index) == {"GENE1", "GENE2"}
    assert set(combined["mutation_category"]) == {"missense"}
    assert set(combined["mode"]) == {"3D"}

    assert combined.loc["GENE1", "mode"] == "3D"
    assert combined.loc["GENE1", "gene_symbol"] == "SYM1"
    assert combined.loc["GENE1", "events_n"] == 10
    assert combined.loc["GENE1", "dist"] == 1.11
    assert combined.loc["GENE1", "probability"] == 0.001

    assert combined.loc["GENE2", "mode"] == "3D"
    assert combined.loc["GENE2", "gene_symbol"] == "SYM2"
    assert combined.loc["GENE2", "events_n"] == 20
    assert combined.loc["GENE2", "dist"] == 2.22
    assert combined.loc["GENE2", "probability"] == 0.002


def test_combine_test_results_falls_back_when_3d_missing_or_probability_nan(tmp_path):
    linear_df = pd.DataFrame(
        [
            {
                "gene_id": "GENE_MISSING",
                "gene_symbol": "SYMM",
                "mutation_category": "missense",
                "events_n": 15,
                "dist": 0.15,
                "probability": 0.015,
            },
            {
                "gene_id": "GENE_NAN",
                "gene_symbol": "SYMN",
                "mutation_category": "missense",
                "events_n": 25,
                "dist": 0.25,
                "probability": 0.025,
            },
            {
                "gene_id": "GENE_3D",
                "gene_symbol": "SYM3",
                "mutation_category": "missense",
                "events_n": 35,
                "dist": 0.35,
                "probability": 0.035,
            },
        ]
    ).set_index("gene_id")

    clustering_3d_df = pd.DataFrame(
        [
            {
                "gene_id": "GENE_NAN",
                "gene_symbol": "SYMN",
                "mutation_category": "missense",
                "events_n": 25,
                "dist": 2.5,
                "probability": float("nan"),
            },
            {
                "gene_id": "GENE_3D",
                "gene_symbol": "SYM3",
                "mutation_category": "missense",
                "events_n": 35,
                "dist": 3.5,
                "probability": 0.0035,
            },
        ]
    ).set_index("gene_id")

    result, combined = _run_combine_test_results(tmp_path, linear_df, clustering_3d_df)

    assert result.exit_code == 0
    assert "Number of genes with 3D clustering results: 1" in result.output
    assert "Number of genes with linear clustering results: 2" in result.output

    assert set(combined.index) == {"GENE_MISSING", "GENE_NAN", "GENE_3D"}
    assert combined.loc["GENE_MISSING", "mode"] == "linear"
    assert combined.loc["GENE_MISSING", "dist"] == 0.15
    assert combined.loc["GENE_MISSING", "probability"] == 0.015

    assert combined.loc["GENE_NAN", "mode"] == "linear"
    assert combined.loc["GENE_NAN", "gene_symbol"] == "SYMN"
    assert combined.loc["GENE_NAN", "events_n"] == 25
    assert combined.loc["GENE_NAN", "dist"] == 0.25
    assert combined.loc["GENE_NAN", "probability"] == 0.025

    assert combined.loc["GENE_3D", "mode"] == "3D"
    assert combined.loc["GENE_3D", "gene_symbol"] == "SYM3"
    assert combined.loc["GENE_3D", "events_n"] == 35
    assert combined.loc["GENE_3D", "dist"] == 3.5
    assert combined.loc["GENE_3D", "probability"] == 0.0035
