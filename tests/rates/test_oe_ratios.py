"""Tests for denovowest.rates.oe_ratios."""

import json
from pathlib import Path

import pandas as pd
import pytest
from click.testing import CliRunner

from denovowest.rates.oe_ratios import (
    _compute_expected_dnm_per_gene,
    _compute_observed_dnm_per_gene,
    _export_oe_ratios,
    _load_gene_infos,
    _oe_ratio_per_category,
    compute_oe_ratios,
    main,
)
from denovowest.utils.params import ConsequenceGroups

###########
# Fixtures #
###########


_MINIMAL_GFF = """\
##gff-version 3
chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=ENSG00000000001;hgnc_id=HGNC:1234
chrX\t.\tgene\t1000\t2000\t.\t+\t.\tID=ENSG00000000002
"""


@pytest.fixture(scope="module")
def minimal_gff_path(tmp_path_factory):
    """Write a two-gene (autosomal + X-linked) GFF to a temp directory."""
    gff_dir = tmp_path_factory.mktemp("gff")
    gff = gff_dir / "test.gff"
    gff.write_text(_MINIMAL_GFF)
    return gff


##########################################
# Unit – _load_gene_infos               #
##########################################


@pytest.mark.unit
def test_load_gene_infos_index_contains_all_genes(minimal_gff_path):
    df = _load_gene_infos(minimal_gff_path)

    assert "ENSG00000000001" in df.index
    assert "ENSG00000000002" in df.index


@pytest.mark.unit
def test_load_gene_infos_extracts_correct_chrom(minimal_gff_path):
    df = _load_gene_infos(minimal_gff_path)

    assert df.loc["ENSG00000000001", "chrom"] == "chr1"
    assert df.loc["ENSG00000000002", "chrom"] == "chrX"


@pytest.mark.unit
def test_load_gene_infos_extracts_hgnc_id(minimal_gff_path):
    df = _load_gene_infos(minimal_gff_path)

    assert df.loc["ENSG00000000001", "hgnc_id"] == "HGNC:1234"


@pytest.mark.unit
def test_load_gene_infos_missing_hgnc_id_is_empty_string(minimal_gff_path):
    df = _load_gene_infos(minimal_gff_path)

    assert df.loc["ENSG00000000002", "hgnc_id"] == ""


##################################################
# Unit – _compute_observed_dnm_per_gene          #
##################################################


@pytest.mark.unit
def test_observed_dnm_counts_per_gene_per_consequence():
    dnm = pd.DataFrame({
        "gene_id":     ["G1", "G1", "G1", "G2"],
        "consequence": ["synonymous", "synonymous", "missense", "synonymous"],
    })

    result = _compute_observed_dnm_per_gene(dnm)

    assert result.loc["G1", "synonymous"] == 2.0
    assert result.loc["G1", "missense"]   == 1.0
    assert result.loc["G2", "synonymous"] == 1.0


@pytest.mark.unit
def test_observed_dnm_missing_consequence_filled_with_zero():
    dnm = pd.DataFrame({
        "gene_id":     ["G1"],
        "consequence": ["synonymous"],
    })

    result = _compute_observed_dnm_per_gene(dnm)

    assert result.loc["G1", "missense"] == 0.0


@pytest.mark.unit
def test_observed_dnm_empty_input_returns_empty_dataframe():
    dnm = pd.DataFrame(columns=["gene_id", "consequence"])

    result = _compute_observed_dnm_per_gene(dnm)

    assert result.empty


##################################################
# Unit – _compute_expected_dnm_per_gene          #
##################################################


@pytest.mark.unit
def test_expected_dnm_autosomal_rate_scaled_by_cohort_size():
    rates = pd.DataFrame({"synonymous": [0.1]}, index=pd.Index(["G1"]))
    gene_infos = pd.DataFrame({"chrom": ["chr1"], "hgnc_id": [""]}, index=pd.Index(["G1"]))

    result = _compute_expected_dnm_per_gene(rates, gene_infos, nmales=50, nfemales=50)

    # rate × 2 × (nmales + nfemales) = 0.1 × 2 × 100 = 20
    assert result.loc["G1", "synonymous"] == pytest.approx(20.0)


@pytest.mark.unit
def test_expected_dnm_x_linked_rate_differs_from_autosomal():
    rates = pd.DataFrame(
        {"synonymous": [0.1, 0.1]},
        index=pd.Index(["autosomal_gene", "x_gene"]),
    )
    gene_infos = pd.DataFrame(
        {"chrom": ["chr1", "chrX"], "hgnc_id": ["", ""]},
        index=pd.Index(["autosomal_gene", "x_gene"]),
    )

    result = _compute_expected_dnm_per_gene(rates, gene_infos, nmales=50, nfemales=50)

    # The X-linked gene should receive a different correction than the autosomal one
    assert result.loc["autosomal_gene", "synonymous"] != pytest.approx(
        result.loc["x_gene", "synonymous"]
    )


#######################################
# Unit – compute_oe_ratios            #
#######################################


@pytest.mark.unit
def test_compute_oe_ratios_perfect_calibration_gives_one():
    observed = pd.DataFrame({"synonymous": [2.0, 4.0]}, index=["G1", "G2"])
    expected = pd.DataFrame({"synonymous": [2.0, 4.0]}, index=["G1", "G2"])

    oe_ratios, _ = compute_oe_ratios(observed, expected)

    assert oe_ratios["G1"] == pytest.approx(1.0)
    assert oe_ratios["G2"] == pytest.approx(1.0)


@pytest.mark.unit
def test_compute_oe_ratios_zero_observed_gives_zero():
    observed = pd.DataFrame({"synonymous": [0.0]}, index=["G1"])
    expected = pd.DataFrame({"synonymous": [2.0]}, index=["G1"])

    oe_ratios, _ = compute_oe_ratios(observed, expected)

    assert oe_ratios["G1"] == pytest.approx(0.0)


@pytest.mark.unit
def test_compute_oe_ratios_pads_genes_absent_from_observed():
    """Genes in the expected table but not the observed table are treated as 0 observed."""
    observed = pd.DataFrame({"synonymous": [2.0]}, index=["G1"])
    expected = pd.DataFrame({"synonymous": [2.0, 4.0]}, index=["G1", "G2"])

    oe_ratios, _ = compute_oe_ratios(observed, expected)

    assert oe_ratios["G2"] == pytest.approx(0.0)


@pytest.mark.unit
def test_compute_oe_ratios_per_consequence_independent():
    observed = pd.DataFrame({"synonymous": [2.0], "missense": [1.0]}, index=["G1"])
    expected = pd.DataFrame({"synonymous": [4.0], "missense": [1.0]}, index=["G1"])

    _, oe_per_csq = compute_oe_ratios(observed, expected)

    assert oe_per_csq.loc["G1", "synonymous"] == pytest.approx(0.5)
    assert oe_per_csq.loc["G1", "missense"]   == pytest.approx(1.0)


@pytest.mark.unit
def test_compute_oe_ratios_aggregate_sums_across_consequences():
    observed = pd.DataFrame({"synonymous": [1.0], "missense": [3.0]}, index=["G1"])
    expected = pd.DataFrame({"synonymous": [2.0], "missense": [2.0]}, index=["G1"])

    oe_ratios, _ = compute_oe_ratios(observed, expected)

    # (1 + 3) / (2 + 2) = 1.0
    assert oe_ratios["G1"] == pytest.approx(1.0)


#########################################
# Unit – _oe_ratio_per_category         #
#########################################


@pytest.fixture
def category_data():
    """Two categories: 'constrained' (G1, G2) and 'unconstrained' (G3)."""
    observed = pd.DataFrame({"synonymous": [2.0, 4.0, 6.0]}, index=["G1", "G2", "G3"])
    expected = pd.DataFrame({"synonymous": [4.0, 4.0, 6.0]}, index=["G1", "G2", "G3"])
    categories = pd.DataFrame(
        {"category": ["constrained", "constrained", "unconstrained"]},
        index=pd.Index(["G1", "G2", "G3"], name="gene_id"),
    )
    return observed, expected, categories


@pytest.mark.unit
def test_oe_ratio_per_category_correct_aggregation(category_data):
    observed, expected, categories = category_data

    df = _oe_ratio_per_category(observed, expected, categories)

    # constrained: (2 + 4) / (4 + 4) = 0.75
    assert df.loc["constrained",   "oe_ratio"] == pytest.approx(0.75)
    # unconstrained: 6 / 6 = 1.0
    assert df.loc["unconstrained", "oe_ratio"] == pytest.approx(1.0)


@pytest.mark.unit
def test_oe_ratio_per_category_nb_genes_is_correct(category_data):
    observed, expected, categories = category_data

    df = _oe_ratio_per_category(observed, expected, categories)

    assert df.loc["constrained",   "nb_genes_in_category"] == 2
    assert df.loc["unconstrained", "nb_genes_in_category"] == 1


@pytest.mark.unit
def test_oe_ratio_per_category_ci_lower_leq_oe(category_data):
    observed, expected, categories = category_data

    df = _oe_ratio_per_category(observed, expected, categories)

    assert (df["ci.lower"] <= df["oe_ratio"]).all()


@pytest.mark.unit
def test_oe_ratio_per_category_ci_upper_geq_oe(category_data):
    observed, expected, categories = category_data

    df = _oe_ratio_per_category(observed, expected, categories)

    assert (df["ci.upper"] >= df["oe_ratio"]).all()


@pytest.mark.unit
def test_oe_ratio_per_category_has_required_columns(category_data):
    observed, expected, categories = category_data

    df = _oe_ratio_per_category(observed, expected, categories)

    for col in ["oe_ratio", "nb_genes_in_category", "sum_observed", "sum_expected", "ci.lower", "ci.upper", "category"]:
        assert col in df.columns


######################################
# Unit – _export_oe_ratios           #
######################################


@pytest.fixture
def oe_data():
    oe_ratios = pd.Series({"G1": 0.8, "G2": 1.2})
    oe_per_csq = pd.DataFrame({
        "synonymous": {"G1": 0.9, "G2": 1.1},
        "missense":   {"G1": 0.7, "G2": 1.3},
    })
    return oe_ratios, oe_per_csq


@pytest.mark.unit
def test_export_oe_ratios_creates_json(tmp_path, oe_data):
    oe_ratios, oe_per_csq = oe_data

    _export_oe_ratios(oe_ratios, oe_per_csq, tmp_path)

    assert (tmp_path / "oe_ratios.json").exists()


@pytest.mark.unit
def test_export_oe_ratios_all_key_matches_aggregate(tmp_path, oe_data):
    oe_ratios, oe_per_csq = oe_data
    _export_oe_ratios(oe_ratios, oe_per_csq, tmp_path)

    data = json.loads((tmp_path / "oe_ratios.json").read_text())

    assert data["G1"]["all"] == pytest.approx(0.8)
    assert data["G2"]["all"] == pytest.approx(1.2)


@pytest.mark.unit
def test_export_oe_ratios_per_consequence_values_are_correct(tmp_path, oe_data):
    oe_ratios, oe_per_csq = oe_data
    _export_oe_ratios(oe_ratios, oe_per_csq, tmp_path)

    data = json.loads((tmp_path / "oe_ratios.json").read_text())

    assert data["G1"]["synonymous"] == pytest.approx(0.9)
    assert data["G1"]["missense"]   == pytest.approx(0.7)
    assert data["G2"]["synonymous"] == pytest.approx(1.1)
    assert data["G2"]["missense"]   == pytest.approx(1.3)


@pytest.mark.unit
def test_export_oe_ratios_all_genes_present(tmp_path, oe_data):
    oe_ratios, oe_per_csq = oe_data
    _export_oe_ratios(oe_ratios, oe_per_csq, tmp_path)

    data = json.loads((tmp_path / "oe_ratios.json").read_text())

    assert set(data.keys()) == {"G1", "G2"}


##############################
# CLI – parameter validation #
##############################


@pytest.fixture
def runner():
    return CliRunner()


@pytest.mark.unit
def test_missing_all_arguments_exits_nonzero(runner):
    result = runner.invoke(main, [])

    assert result.exit_code != 0


@pytest.mark.unit
def test_nonexistent_dnm_file_exits_nonzero(runner, tmp_path, minimal_gff_path):
    rates_path = tmp_path / "rates.tsv"
    pd.DataFrame(
        {csq: [0.01] for csq in ConsequenceGroups.OE_RATIO},
        index=pd.Index(["ENSG00000000001"], name="gene_id"),
    ).to_csv(rates_path, sep="\t")

    result = runner.invoke(main, [
        str(tmp_path / "missing_dnm.tsv"),
        str(rates_path),
        str(minimal_gff_path),
        "50", "50",
        str(tmp_path),
    ])

    assert result.exit_code != 0


################################
# Functional – CLI end-to-end  #
################################


@pytest.fixture(scope="module")
def oe_cli_output(tmp_path_factory, minimal_gff_path):
    """Run the CLI on minimal synthetic data and return the parsed oe_ratios.json."""
    tmp = tmp_path_factory.mktemp("oe_cli")

    dnm_path = tmp / "dnm.tsv"
    pd.DataFrame([
        {"gene_id": "ENSG00000000001", "chrom": "chr1", "pos": 1500, "ref": "A", "alt": "T", "consequence": "synonymous"},
        {"gene_id": "ENSG00000000001", "chrom": "chr1", "pos": 1600, "ref": "G", "alt": "C", "consequence": "missense"},
    ]).to_csv(dnm_path, sep="\t", index=False)

    rates_path = tmp / "rates.tsv"
    pd.DataFrame(
        {csq: [0.01] for csq in ConsequenceGroups.OE_RATIO},
        index=pd.Index(["ENSG00000000001"], name="gene_id"),
    ).to_csv(rates_path, sep="\t")

    outdir = tmp / "out"
    outdir.mkdir()

    result = CliRunner().invoke(main, [
        str(dnm_path),
        str(rates_path),
        str(minimal_gff_path),
        "50", "50",
        str(outdir),
    ])
    assert result.exit_code == 0, result.output

    return json.loads((outdir / "oe_ratios.json").read_text())


@pytest.mark.functional
class TestOeRatiosCLIOutput:

    def test_output_contains_gene(self, oe_cli_output):
        assert "ENSG00000000001" in oe_cli_output

    def test_output_has_all_key(self, oe_cli_output):
        assert "all" in oe_cli_output["ENSG00000000001"]

    def test_output_has_per_consequence_keys(self, oe_cli_output):
        gene = oe_cli_output["ENSG00000000001"]
        assert "synonymous" in gene
        assert "missense" in gene

    def test_aggregate_ratio_is_positive(self, oe_cli_output):
        assert oe_cli_output["ENSG00000000001"]["all"] > 0

    def test_synonymous_ratio_is_positive(self, oe_cli_output):
        assert oe_cli_output["ENSG00000000001"]["synonymous"] > 0
