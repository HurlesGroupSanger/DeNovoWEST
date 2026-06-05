"""Tests for denovowest.rates.create_rates_file."""

from pathlib import Path

import pandas as pd
import pytest
from click.testing import CliRunner

from denovowest.rates.create_rates_file import (
    get_alternates,
    load_gene_list,
    load_kmer_mutation_rate_model,
    main,
)

#############################
# Paths to fixtures         #
#############################

_ROOT = Path(__file__).parent.parent.parent

KMER_INPUT = _ROOT / "examples" / "rates" / "create_rates_kmer" / "input"
KMER_OUTPUT = _ROOT / "examples" / "rates" / "create_rates_kmer" / "output"
ROULETTE_INPUT = _ROOT / "examples" / "rates" / "create_rates_roulette" / "input"
ROULETTE_OUTPUT = _ROOT / "examples" / "rates" / "create_rates_roulette" / "output"

GFF = KMER_INPUT / "OR4F5.gff"
FASTA = KMER_INPUT / "OR4F5.fa"
KMER_MODEL = _ROOT / "resources" / "rates_model_dnn.tsv"
ROULETTE_VCF_DIR = ROULETTE_INPUT / "roulette_OR4F5"

ROULETTE_SCALING_FACTOR = "5.075e-08"  # 1.015e-7 / 2

#############################
# Unit tests                #
#############################


@pytest.mark.unit
def test_get_alternates_returns_three_kmers():
    alts = get_alternates("ACA", range_model=1)

    assert len(alts) == 3
    assert "ACA" not in alts


@pytest.mark.unit
def test_get_alternates_only_changes_central_nucleotide():
    alts = get_alternates("ACA", range_model=1)

    # Flanking bases must be unchanged; only the central "C" varies
    assert all(a[0] == "A" and a[2] == "A" for a in alts)
    assert {a[1] for a in alts} == {"A", "G", "T"}


@pytest.mark.unit
def test_load_kmer_mutation_rate_model_index_format():
    model = load_kmer_mutation_rate_model(KMER_MODEL)

    assert "mu_snp" in model.columns
    assert model.index[0] == "AAA_ACA"
    assert (model["mu_snp"] > 0).all()


@pytest.mark.unit
def test_load_gene_list_from_file(tmp_path):
    gene_file = tmp_path / "genes.txt"
    gene_file.write_text("ENSG00000001\nENSG00000002\n")

    result = load_gene_list(gene_file, gff_db=None)

    assert result == ["ENSG00000001", "ENSG00000002"]


@pytest.mark.unit
def test_load_gene_list_strips_trailing_newline(tmp_path):
    gene_file = tmp_path / "genes.txt"
    gene_file.write_text("ENSG00000001\n")

    result = load_gene_list(gene_file, gff_db=None)

    assert result == ["ENSG00000001"]


#############################
# CLI parameter validation  #
#############################


@pytest.fixture
def runner():
    return CliRunner()


@pytest.mark.unit
def test_missing_model_exits_nonzero(runner, tmp_path):
    result = runner.invoke(main, [
        "--rates_model_path", str(KMER_MODEL),
        "--gff", str(GFF),
        "--outdir", str(tmp_path),
    ])

    assert result.exit_code != 0


@pytest.mark.unit
def test_invalid_model_choice_exits_nonzero(runner, tmp_path):
    result = runner.invoke(main, [
        "--rates_model_path", str(KMER_MODEL),
        "--model", "not_a_model",
        "--gff", str(GFF),
        "--outdir", str(tmp_path),
    ])

    assert result.exit_code != 0
    assert "invalid value" in result.output.lower()


@pytest.mark.unit
def test_kmer_without_fasta_exits_with_usage_error(runner, tmp_path):
    result = runner.invoke(main, [
        "--rates_model_path", str(KMER_MODEL),
        "--model", "kmer",
        "--gff", str(GFF),
        "--outdir", str(tmp_path),
    ])

    assert result.exit_code != 0
    assert "--fasta" in result.output


@pytest.mark.unit
def test_roulette_without_scaling_factor_exits_with_usage_error(runner, tmp_path):
    result = runner.invoke(main, [
        "--rates_model_path", str(ROULETTE_VCF_DIR),
        "--model", "roulette",
        "--gff", str(GFF),
        "--outdir", str(tmp_path),
    ])

    assert result.exit_code != 0
    assert "--scaling_factor" in result.output


@pytest.mark.unit
def test_carlson_without_scaling_factor_exits_with_usage_error(runner, tmp_path):
    result = runner.invoke(main, [
        "--rates_model_path", str(ROULETTE_VCF_DIR),
        "--model", "carlson",
        "--gff", str(GFF),
        "--outdir", str(tmp_path),
    ])

    assert result.exit_code != 0
    assert "--scaling_factor" in result.output


@pytest.mark.unit
def test_nonexistent_gff_exits_nonzero(runner, tmp_path):
    result = runner.invoke(main, [
        "--rates_model_path", str(KMER_MODEL),
        "--model", "kmer",
        "--gff", str(tmp_path / "does_not_exist.gff"),
        "--fasta", str(FASTA),
        "--outdir", str(tmp_path),
    ])

    assert result.exit_code != 0


#############################
# Functional / non-regression
#############################


@pytest.fixture(scope="module")
def kmer_output(tmp_path_factory):
    """Run the kmer model once for the whole module and return the output DataFrame."""
    out = tmp_path_factory.mktemp("kmer")
    result = CliRunner().invoke(main, [
        "--rates_model_path", str(KMER_MODEL),
        "--model", "kmer",
        "--gff", str(GFF),
        "--fasta", str(FASTA),
        "--outdir", str(out),
    ])
    assert result.exit_code == 0, result.output
    return pd.read_csv(out / "mutation_rates.tsv", sep="\t")


@pytest.fixture(scope="module")
def roulette_vcf_dir(tmp_path_factory):
    """Return a directory whose VCF files follow the ``*all.vcf*.gz`` naming convention.

    The example file has a non-standard name; we expose it via symlinks so the
    production glob pattern is exercised without renaming the test data.
    """
    src_vcf = next(ROULETTE_VCF_DIR.glob("*.vcf.gz"))
    src_tbi = Path(str(src_vcf) + ".tbi")

    roulette_dir = tmp_path_factory.mktemp("roulette_vcf")
    (roulette_dir / "1_all.vcf.gz").symlink_to(src_vcf)
    (roulette_dir / "1_all.vcf.gz.tbi").symlink_to(src_tbi)

    return roulette_dir


@pytest.fixture(scope="module")
def roulette_output(tmp_path_factory, roulette_vcf_dir):
    """Run the roulette model once for the whole module and return the output DataFrame."""
    out = tmp_path_factory.mktemp("roulette")
    result = CliRunner().invoke(main, [
        "--rates_model_path", str(roulette_vcf_dir),
        "--model", "roulette",
        "--scaling_factor", ROULETTE_SCALING_FACTOR,
        "--gff", str(GFF),
        "--outdir", str(out),
    ])
    assert result.exit_code == 0, result.output
    return pd.read_csv(out / "mutation_rates.tsv", sep="\t")


@pytest.mark.functional
class TestKmerOutput:

    def test_exit_code_is_zero(self, kmer_output):
        # If the fixture didn't raise, the run succeeded
        assert kmer_output is not None

    def test_has_expected_columns(self, kmer_output):
        assert list(kmer_output.columns) == ["gene_id", "chrom", "pos", "ref", "alt", "prob"]

    def test_gene_is_OR4F5(self, kmer_output):
        assert set(kmer_output.gene_id) == {"ENSG00000186092.7"}

    def test_chromosome_is_chr1(self, kmer_output):
        assert (kmer_output.chrom == "chr1").all()

    def test_probabilities_are_positive(self, kmer_output):
        assert (kmer_output.prob > 0).all()

    def test_ref_and_alt_are_single_nucleotides(self, kmer_output):
        assert (kmer_output.ref.str.len() == 1).all()
        assert (kmer_output.alt.str.len() == 1).all()

    def test_no_duplicate_variants(self, kmer_output):
        duplicates = kmer_output[["gene_id", "chrom", "pos", "ref", "alt"]].duplicated()
        assert not duplicates.any()

    def test_non_regression(self, kmer_output):
        reference = pd.read_csv(KMER_OUTPUT / "mutation_rates.tsv", sep="\t")
        pd.testing.assert_frame_equal(kmer_output.reset_index(drop=True), reference.reset_index(drop=True))


@pytest.mark.functional
class TestRouletteOutput:

    def test_exit_code_is_zero(self, roulette_output):
        assert roulette_output is not None

    def test_has_expected_columns(self, roulette_output):
        assert list(roulette_output.columns) == ["gene_id", "chrom", "pos", "ref", "alt", "prob"]

    def test_gene_is_OR4F5(self, roulette_output):
        assert set(roulette_output.gene_id) == {"ENSG00000186092.7"}

    def test_probabilities_are_positive(self, roulette_output):
        assert (roulette_output.prob > 0).all()

    def test_probabilities_are_scaled(self, roulette_output):
        # All prob values should be in the ballpark of scaling_factor * raw_roulette_rates
        # Raw Roulette rates are ~1e-8 to 1e-6; after scaling by 5e-8 they should be <1e-9
        assert (roulette_output.prob < 1e-6).all()

    def test_no_duplicate_variants(self, roulette_output):
        duplicates = roulette_output[["gene_id", "chrom", "pos", "ref", "alt"]].duplicated()
        assert not duplicates.any()

    def test_non_regression(self, roulette_output):
        reference = pd.read_csv(ROULETTE_OUTPUT / "mutation_rates.tsv", sep="\t")
        pd.testing.assert_frame_equal(roulette_output.reset_index(drop=True), reference.reset_index(drop=True))
