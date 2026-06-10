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

#####################
# Regression data   #
#####################

_NMALES, _NFEMALES = 50, 50
_N_SITES = 10  # per-site rows per (gene, consequence) in the rates file

# 10 synthetic genes: (gene_id, chrom, category, hgnc_id, start, end, strand, size_multiplier)
_GENES_FULL = [
    ("ENSG00000000001", "chr1",  "constrained",   "HGNC:101",  1000, 11000, "+", 1.00),
    ("ENSG00000000002", "chr1",  "constrained",   "HGNC:102", 15000, 25000, "+", 1.20),
    ("ENSG00000000003", "chr1",  "constrained",   "HGNC:103", 30000, 40000, "-", 0.80),
    ("ENSG00000000004", "chr1",  "moderate",      "HGNC:104", 45000, 55000, "+", 1.50),
    ("ENSG00000000005", "chr1",  "moderate",      "HGNC:105", 60000, 70000, "-", 1.00),
    ("ENSG00000000006", "chr2",  "moderate",      "HGNC:106",  1000, 11000, "+", 1.10),
    ("ENSG00000000007", "chr2",  "moderate",      "HGNC:107", 15000, 25000, "+", 0.90),
    ("ENSG00000000008", "chr2",  "unconstrained", "HGNC:108", 30000, 40000, "-", 1.30),
    ("ENSG00000000009", "chr2",  "unconstrained", "HGNC:109", 45000, 55000, "+", 0.70),
    ("ENSG00000000010", "chrX",  "unconstrained", "HGNC:110",  1000, 11000, "+", 1.00),
]

_GENE_CHROM: dict[str, str]  = {g[0]: g[1] for g in _GENES_FULL}
_GENE_MULT:  dict[str, float] = {g[0]: g[7] for g in _GENES_FULL}

# Base prob_sum per consequence for a ×1.0 gene (spread across _N_SITES per-site rows)
_BASE_PROB_SUM: dict[str, float] = {
    "synonymous":      5.00e-3,
    "missense":        1.00e-2,
    "stop_gained":     5.00e-4,
    "stop_retained":   1.00e-4,
    "start_lost":      5.00e-5,
    "stop_lost":       5.00e-5,
    "splice_region":   5.00e-4,
    "splice_acceptor": 1.00e-4,
    "splice_donor":    1.00e-4,
}

# Observed de novo mutations: (gene_id, consequence, ref, alt, pos)
_DNM_SPEC = [
    # constrained gene 1: missense ×4, stop_gained ×1, splice_donor ×1
    ("ENSG00000000001", "missense",        "G", "A", 2100),
    ("ENSG00000000001", "missense",        "C", "T", 2200),
    ("ENSG00000000001", "missense",        "G", "C", 2300),
    ("ENSG00000000001", "missense",        "T", "A", 2400),
    ("ENSG00000000001", "stop_gained",     "C", "T", 3100),
    ("ENSG00000000001", "splice_donor",    "A", "G", 3200),
    # constrained gene 2: synonymous ×1, missense ×5, stop_gained ×1
    ("ENSG00000000002", "synonymous",      "T", "C", 16100),
    ("ENSG00000000002", "missense",        "G", "A", 16200),
    ("ENSG00000000002", "missense",        "C", "G", 16300),
    ("ENSG00000000002", "missense",        "A", "T", 16400),
    ("ENSG00000000002", "missense",        "G", "C", 16500),
    ("ENSG00000000002", "missense",        "T", "A", 16600),
    ("ENSG00000000002", "stop_gained",     "C", "T", 17100),
    # constrained gene 3: missense ×3, splice_acceptor ×1
    ("ENSG00000000003", "missense",        "G", "A", 31100),
    ("ENSG00000000003", "missense",        "C", "T", 31200),
    ("ENSG00000000003", "missense",        "T", "G", 31300),
    ("ENSG00000000003", "splice_acceptor", "A", "C", 32100),
    # moderate gene 4: synonymous ×2, missense ×3
    ("ENSG00000000004", "synonymous",      "T", "C", 46100),
    ("ENSG00000000004", "synonymous",      "A", "G", 46200),
    ("ENSG00000000004", "missense",        "G", "A", 46300),
    ("ENSG00000000004", "missense",        "C", "T", 46400),
    ("ENSG00000000004", "missense",        "G", "C", 46500),
    # moderate gene 5: synonymous ×1, missense ×2
    ("ENSG00000000005", "synonymous",      "A", "T", 61100),
    ("ENSG00000000005", "missense",        "G", "C", 61200),
    ("ENSG00000000005", "missense",        "C", "A", 61300),
    # moderate gene 6: synonymous ×1, missense ×2
    ("ENSG00000000006", "synonymous",      "T", "A", 2100),
    ("ENSG00000000006", "missense",        "G", "T", 2200),
    ("ENSG00000000006", "missense",        "A", "C", 2300),
    # moderate gene 7: synonymous ×1, missense ×2
    ("ENSG00000000007", "synonymous",      "C", "T", 16100),
    ("ENSG00000000007", "missense",        "G", "A", 16200),
    ("ENSG00000000007", "missense",        "T", "C", 16300),
    # unconstrained gene 8: synonymous ×2, missense ×2
    ("ENSG00000000008", "synonymous",      "A", "G", 31100),
    ("ENSG00000000008", "synonymous",      "T", "C", 31200),
    ("ENSG00000000008", "missense",        "G", "A", 31300),
    ("ENSG00000000008", "missense",        "C", "T", 31400),
    # unconstrained gene 9: synonymous ×1, missense ×1
    ("ENSG00000000009", "synonymous",      "A", "T", 46100),
    ("ENSG00000000009", "missense",        "C", "G", 46200),
    # unconstrained gene 10 (chrX): synonymous ×1, missense ×1
    ("ENSG00000000010", "synonymous",      "G", "C", 2100),
    ("ENSG00000000010", "missense",        "T", "A", 2200),
]

# ── Helpers for computing analytically-derived expected O/E values ─────────────


def _autosomal_factor() -> float:
    return 2.0 * (_NMALES + _NFEMALES)


def _x_factor() -> float:
    """X-chromosome scaling factor (same formula as compute_x_factor_correction)."""
    alpha = 3.4
    male_k = 2 / (1 + 1 / alpha)
    female_k = 2 / (1 + alpha)
    female_transmit = _NMALES + _NFEMALES
    male_transmit = _NFEMALES
    return (male_transmit * male_k + female_transmit * female_k) / _autosomal_factor()


def _expected_counts(gene_id: str) -> dict[str, float]:
    """Expected mutation count per consequence for one gene, given the cohort size."""
    factor = _autosomal_factor() * (_x_factor() if _GENE_CHROM[gene_id] == "chrX" else 1.0)
    return {csq: _BASE_PROB_SUM[csq] * _GENE_MULT[gene_id] * factor for csq in _BASE_PROB_SUM}


def _observed_counts(gene_id: str) -> dict[str, int]:
    """Observed mutation counts per consequence for one gene, read from _DNM_SPEC."""
    obs: dict[str, int] = {csq: 0 for csq in _BASE_PROB_SUM}
    for gid, csq, *_ in _DNM_SPEC:
        if gid == gene_id:
            obs[csq] += 1
    return obs


# ── Helpers for generating on-the-fly regression inputs ──────────────────────


def _write_full_gff(path: Path) -> None:
    lines = ["##gff-version 3"]
    for gene_id, chrom, _, hgnc, start, end, strand, _ in _GENES_FULL:
        lines.append(f"{chrom}\t.\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gene_id};hgnc_id={hgnc}")
    path.write_text("\n".join(lines) + "\n")


def _write_full_rates_matrix(path: Path) -> None:
    pd.DataFrame(
        {gene_id: {csq: _BASE_PROB_SUM[csq] * mult for csq in _BASE_PROB_SUM}
         for gene_id, mult in _GENE_MULT.items()}
    ).transpose().to_csv(path, sep="\t")


def _write_full_dnm(path: Path) -> None:
    pd.DataFrame([
        {"gene_id": gid, "chrom": _GENE_CHROM[gid], "pos": pos,
         "ref": ref, "alt": alt, "consequence": csq}
        for gid, csq, ref, alt, pos in _DNM_SPEC
    ]).to_csv(path, sep="\t", index=False)


def _write_full_categories(path: Path) -> None:
    pd.DataFrame(
        [{"gene_id": g[0], "category": g[2]} for g in _GENES_FULL]
    ).set_index("gene_id").to_csv(path, sep="\t")


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


##############
# Regression #
##############


@pytest.fixture(scope="module")
def regression_oe_output(tmp_path_factory):
    """Run the full CLI pipeline on the 10-gene dataset generated on-the-fly."""
    tmp = tmp_path_factory.mktemp("regression_oe")
    gff_path        = tmp / "genes.gff"
    rates_path      = tmp / "rates.tsv"
    dnm_path        = tmp / "dnm.tsv"
    categories_path = tmp / "gene_categories.tsv"
    outdir          = tmp / "out"
    outdir.mkdir()

    _write_full_gff(gff_path)
    _write_full_rates_matrix(rates_path)
    _write_full_dnm(dnm_path)
    _write_full_categories(categories_path)

    result = CliRunner().invoke(main, [
        str(dnm_path), str(rates_path), str(gff_path),
        str(_NMALES), str(_NFEMALES),
        str(outdir),
        "--gene-categories", str(categories_path),
    ])
    assert result.exit_code == 0, result.output
    return json.loads((outdir / "oe_ratios.json").read_text())


@pytest.mark.regression
class TestOeRatiosNonRegression:
    """Verify O/E ratios against values derived analytically from the test-data constants.

    Expected O/E = observed / expected, where:
      - observed is read from _DNM_SPEC
      - expected = BASE_PROB_SUM × multiplier × cohort_factor
        with cohort_factor = autosomal_factor × (x_factor if chrX else 1)

    If the pipeline's scaling or aggregation logic changes, these tests will fail.
    """

    def test_all_genes_present_in_output(self, regression_oe_output):
        assert set(regression_oe_output.keys()) == {g[0] for g in _GENES_FULL}

    def test_all_consequence_and_aggregate_keys_present(self, regression_oe_output):
        expected_keys = set(_BASE_PROB_SUM) | {"all"}
        for gene_id in regression_oe_output:
            assert set(regression_oe_output[gene_id].keys()) == expected_keys, gene_id

    def test_aggregate_oe_ratio_matches_analytical_value(self, regression_oe_output):
        for gene_id in _GENE_MULT:
            exp = _expected_counts(gene_id)
            obs = _observed_counts(gene_id)
            expected_oe = sum(obs.values()) / sum(exp.values())
            assert regression_oe_output[gene_id]["all"] == pytest.approx(expected_oe, rel=1e-6), \
                f"Aggregate O/E mismatch for {gene_id}"

    def test_per_consequence_oe_ratio_matches_analytical_value(self, regression_oe_output):
        for gene_id in _GENE_MULT:
            exp = _expected_counts(gene_id)
            obs = _observed_counts(gene_id)
            for csq in _BASE_PROB_SUM:
                expected_oe = obs[csq] / exp[csq]
                assert regression_oe_output[gene_id][csq] == pytest.approx(expected_oe, rel=1e-6), \
                    f"Per-consequence O/E mismatch for {gene_id}.{csq}"

    def test_constrained_genes_have_higher_aggregate_oe_than_unconstrained(self, regression_oe_output):
        constrained   = [g[0] for g in _GENES_FULL if g[2] == "constrained"]
        unconstrained = [g[0] for g in _GENES_FULL if g[2] == "unconstrained"]
        avg_constrained   = sum(regression_oe_output[g]["all"] for g in constrained)   / len(constrained)
        avg_unconstrained = sum(regression_oe_output[g]["all"] for g in unconstrained) / len(unconstrained)
        assert avg_constrained > avg_unconstrained
