import numpy as np
import pandas as pd
import pytest

from denovowest.simulation.config import Config
from denovowest.simulation.probabilities import calc_pn, get_pvalue, sim_score
from denovowest.simulation.simulation import (
    assign_meta_consequences,
    compute_expected_number_of_mutations,
    compute_x_factor_correction,
    extract_worst_consequence,
    filter_on_consequences,
)
from denovowest.utils.params import RunType


def test_extract_worst_consequence_picks_most_severe_term():
    assert extract_worst_consequence("synonymous&missense&splice_donor") == "splice_donor"


def test_filter_on_consequences_for_dnm_keeps_original_and_filters_unknown_terms():
    df = pd.DataFrame(
        {
            "gene_id": ["G1", "G1", "G2"],
            "consequence": ["missense&synonymous", "splice_donor", "intergenic"],
        }
    )

    filtered = filter_on_consequences(df.copy(), "dnm", Config(runtype=RunType.ALL_CODING))

    assert list(filtered["consequence"]) == ["missense", "splice_donor"]
    assert list(filtered["original_consequence"]) == ["missense&synonymous", "splice_donor"]


def test_assign_meta_consequences_maps_to_simulation_categories():
    df = pd.DataFrame({"consequence": ["stop_gained", "splice_donor", "synonymous"]})

    mapped = assign_meta_consequences(df.copy())

    assert list(mapped["consequence"]) == ["nonsense", "splice_lof", "synonymous"]


def test_compute_x_factor_correction_matches_expected_formula():
    correction = compute_x_factor_correction(10, 20)

    alpha = 3.4
    autosomal_factor = 2 * (10 + 20)
    female_transmit = 10 + 20
    male_transmit = 20
    male_k = 2 / (1 + (1 / alpha))
    female_k = 2 / (1 + alpha)
    expected = ((male_transmit * male_k) + (female_transmit * female_k)) / autosomal_factor

    assert correction == pytest.approx(expected)


def test_compute_expected_number_of_mutations_scales_autosomes_and_x():
    rates_df = pd.DataFrame(
        {
            "chrom": ["1", "X", "chrX"],
            "prob": [0.1, 0.1, 0.2],
        }
    )

    scaled = compute_expected_number_of_mutations(rates_df.copy(), nmales=10, nfemales=20)

    autosomal_expected = 0.1 * 2 * (10 + 20)
    x_factor = compute_x_factor_correction(10, 20)

    assert scaled.loc[0, "prob"] == pytest.approx(autosomal_expected)
    assert scaled.loc[1, "prob"] == pytest.approx(autosomal_expected * x_factor)
    assert scaled.loc[2, "prob"] == pytest.approx(0.2 * 2 * (10 + 20) * x_factor)


def test_get_pvalue_returns_one_when_observed_score_is_below_expected():
    generates = pd.DataFrame({"prob": [0.6, 0.4], "score": [2.0, 1.0]})

    pval, expected_score, logs = get_pvalue(
        generates=generates,
        obs_sum_scores=0.5,
        nb_observed_mutations=1,
        score_column="score",
        cfg=Config(),
        simulation_logs={},
    )

    assert pval == 1
    assert expected_score == pytest.approx(1.6)
    assert logs["nb_simulations"] == 0
    assert logs["info"] == "observed < expected, pvalue set at 1"


def test_get_pvalue_returns_one_when_observed_score_is_zero():
    generates = pd.DataFrame({"prob": [0.6, 0.4], "score": [0.0, 0.0]})

    pval, expected_score, logs = get_pvalue(
        generates=generates,
        obs_sum_scores=0,
        nb_observed_mutations=1,
        score_column="score",
        cfg=Config(),
        simulation_logs={},
    )

    assert pval == 1
    assert expected_score == pytest.approx(0.0)
    assert logs["nb_simulations"] == 0
    assert logs["info"] == "observed = 0, pvalue set at 1"


def test_sim_score_uses_explicit_nsim_for_chunking(monkeypatch):
    calls = []

    monkeypatch.setattr(
        "denovowest.simulation.probabilities.delayed",
        lambda func: (lambda *args, **kwargs: func(*args, **kwargs)),
    )

    class FakeParallel:
        def __init__(self, n_jobs):
            self.n_jobs = n_jobs

        def __call__(self, iterable):
            return list(iterable)

    monkeypatch.setattr("denovowest.simulation.probabilities.Parallel", FakeParallel)

    def fake_choice(scores, size, p):
        calls.append(size)
        return np.ones(size)

    monkeypatch.setattr("denovowest.simulation.probabilities.np.random.choice", fake_choice)

    rates = pd.DataFrame({"prob": [0.5, 0.5], "score": [1.0, 2.0]})
    cfg = Config(nsim=250, jobs=1)

    nb_more_extreme_scores = sim_score(
        mu=1.0,
        obs_sum_scores=0.0,
        rates=rates,
        nb_mutation_poisson=2,
        score_column="score",
        nsim=250,
        cfg=cfg,
    )

    assert nb_more_extreme_scores == 250
    assert calls == [(2, 100), (2, 100), (2, 50)]


def test_calc_pn_executes_the_adaptive_simulation_count(monkeypatch):
    simulated_nsim = []

    def fake_sim_score(mu, obs_sum_scores, rates, nb_mutation_poisson, score_column, nsim, cfg):
        simulated_nsim.append(nsim)
        return 0

    monkeypatch.setattr("denovowest.simulation.probabilities.sim_score", fake_sim_score)
    monkeypatch.setattr("denovowest.simulation.probabilities.stats.poisson.pmf", lambda n, mu: 0.5)

    rates = pd.DataFrame(
        {
            "prob": [1.0] * 20,
            "score": [0.0] * 10 + [1.0] * 10,
        }
    )
    cfg = Config(nsim=1_000, jobs=1)
    mu = rates["prob"].sum()

    pn, nsim, s = calc_pn(
        mu=mu,
        obs_sum_scores=5.5,
        rates=rates,
        nb_mutation_poisson=10,
        scores_sorted=np.sort(rates["score"].to_numpy()),
        score_column="score",
        ptot=0.0,
        cfg=cfg,
    )

    assert nsim == 10**7
    assert simulated_nsim == [nsim]
    assert s == 0
    assert pn >= 0
