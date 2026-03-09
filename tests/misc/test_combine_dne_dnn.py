import importlib.util
from pathlib import Path

import pandas as pd
import pytest
from click.testing import CliRunner
from scipy.stats import combine_pvalues

pytestmark = pytest.mark.misc


def _load_combine_dne_dnn_command():
    script_path = Path(__file__).resolve().parents[2] / "misc" / "combine_dne_dnn.py"
    spec = importlib.util.spec_from_file_location("combine_dne_dnn_script", script_path)
    module = importlib.util.module_from_spec(spec)
    assert spec is not None
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module.combine_dne_dnn


def _run_combine_dne_dnn(tmp_path, dne_all_df, dne_mis_df, dnn_df):
    combine_dne_dnn = _load_combine_dne_dnn_command()

    dne_all_path = tmp_path / "dne_all.tsv"
    dne_mis_path = tmp_path / "dne_mis.tsv"
    dnn_path = tmp_path / "dnn.tsv"
    out_path = tmp_path / "combined.tsv"
    dne_all_df.to_csv(dne_all_path, sep="\t", index=False)
    dne_mis_df.to_csv(dne_mis_path, sep="\t", index=False)
    dnn_df.to_csv(dnn_path, sep="\t", index=False)

    runner = CliRunner()
    result = runner.invoke(
        combine_dne_dnn, [str(dne_all_path), str(dne_mis_path), str(dnn_path), str(out_path)]
    )
    combined = pd.read_csv(out_path, sep="\t")
    return result, combined


def test_combine_dne_dnn_merges_sources_and_computes_min_pval_with_mode(tmp_path):
    dne_all_df = pd.DataFrame(
        [
            {"gene_id": "G1", "expected": 1.0, "observed": 2, "enrichment_pval": 0.2},
            {"gene_id": "G2", "expected": 1.5, "observed": 3, "enrichment_pval": 0.5},
            {"gene_id": "G3", "expected": 2.0, "observed": 1, "enrichment_pval": 0.05},
        ]
    )
    dne_mis_df = pd.DataFrame(
        [
            {"gene_id": "G1", "expected": 0.5, "observed": 1, "enrichment_pval": 0.1},
            {"gene_id": "G2", "expected": 0.7, "observed": 2, "enrichment_pval": 0.2},
            {"gene_id": "G4", "expected": 0.4, "observed": 1, "enrichment_pval": 0.3},
        ]
    )
    dnn_df = pd.DataFrame(
        [
            {
                "gene_id": "G1",
                "mutation_category": "missense",
                "probability": 0.01,
                "events_n": 11,
                "dist": 1.1,
                "mode": "3D",
            },
            {
                "gene_id": "G2",
                "mutation_category": "synonymous",
                "probability": 0.0001,
                "events_n": 99,
                "dist": 9.9,
                "mode": "linear",
            },
            {
                "gene_id": "G3",
                "mutation_category": "missense",
                "probability": 0.03,
                "events_n": 33,
                "dist": 3.3,
                "mode": "linear",
            },
            {
                "gene_id": "G5",
                "mutation_category": "missense",
                "probability": 0.02,
                "events_n": 55,
                "dist": 5.5,
                "mode": "3D",
            },
        ]
    )

    result, combined = _run_combine_dne_dnn(tmp_path, dne_all_df, dne_mis_df, dnn_df)

    assert result.exit_code == 0
    assert set(combined["gene_id"]) == {"G1", "G2", "G3", "G4"}
    assert "clustering_mode" in combined.columns

    by_gene = combined.set_index("gene_id")
    expected_g1_combined = combine_pvalues([0.01, 0.1], method="fisher")[1]

    assert by_gene.loc["G1", "clustering_mode"] == "3D"
    assert by_gene.loc["G1", "nb_mis_variants"] == 11
    assert by_gene.loc["G1", "clustering_dist"] == 1.1
    assert by_gene.loc["G1", "combined_mis_pval"] == pytest.approx(expected_g1_combined)
    assert by_gene.loc["G1", "min_pval"] == pytest.approx(min(expected_g1_combined, 0.2))

    assert pd.isna(by_gene.loc["G2", "clustering_pval"])
    assert by_gene.loc["G2", "min_pval"] == pytest.approx(0.5)

    assert pd.isna(by_gene.loc["G3", "enrichment_mis_pval"])
    assert by_gene.loc["G3", "clustering_mode"] == "linear"
    assert by_gene.loc["G3", "min_pval"] == pytest.approx(0.05)

    assert pd.isna(by_gene.loc["G4", "enrichment_all_pval"])
    assert pd.isna(by_gene.loc["G4", "combined_mis_pval"])
    assert by_gene.loc["G4", "min_pval"] == pytest.approx(1.0)


def test_combine_dne_dnn_omits_clustering_mode_when_not_present(tmp_path):
    dne_all_df = pd.DataFrame(
        [{"gene_id": "A1", "expected": 1.0, "observed": 2, "enrichment_pval": 0.4}]
    )
    dne_mis_df = pd.DataFrame(
        [{"gene_id": "A1", "expected": 0.5, "observed": 1, "enrichment_pval": 0.2}]
    )
    dnn_df = pd.DataFrame(
        [
            {
                "gene_id": "A1",
                "mutation_category": "missense",
                "probability": 0.05,
                "events_n": 4,
                "dist": 2.5,
            }
        ]
    )

    result, combined = _run_combine_dne_dnn(tmp_path, dne_all_df, dne_mis_df, dnn_df)

    assert result.exit_code == 0
    assert "clustering_mode" not in combined.columns
    expected_combined = combine_pvalues([0.05, 0.2], method="fisher")[1]
    assert combined.loc[0, "combined_mis_pval"] == pytest.approx(expected_combined)
    assert combined.loc[0, "min_pval"] == pytest.approx(min(expected_combined, 0.4))
