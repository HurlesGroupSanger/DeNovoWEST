import numpy as np
import click

from scipy import stats
from joblib import Parallel, delayed
from denovowest.utils.params import DEFAULT_MAX_NB_MUTATIONS_SIM, DEFAULT_MIN_NB_SIM


def calc_p0(mu, obs_sum_scores):
    """
    Exact calculate P(S >= obs_sum_scores | N = 0)P(N = 0)

    Args:
        mu (float): poisson parameter correpsponding to the sum of the mutation rates
        obs_sum_scores (float): sum of observed DNM scores
    """

    if obs_sum_scores == 0:
        p0 = 1 * stats.poisson.pmf(0, mu)
    else:
        p0 = 0
    return p0


def calc_p1(mu, obs_sum_scores, rates, score_column):
    """
    Exact calculate P(S >= obs_sum_scores | N = 1)P(N = 1)

    Args:
        mu (float): poisson parameter correpsponding to the sum of the mutation rates
        obs_sum_scores (float): sum of observed DNM scores
        rates (pd.DataFrame): all possible mutations annotated
        score_column (str) : CEP scores

    """

    # Get the proportion of putative variants that have a score greater than the observed sum of scores
    p1c = rates["prob"][rates[score_column] >= obs_sum_scores].sum() / rates["prob"].sum()

    # Weight it by the probability of observing one mutation given the poisson rate (expected number of mutations)
    p1 = p1c * stats.poisson.pmf(1, mu)

    return p1


def calc_pn(mu, obs_sum_scores, rates, nb_mutation_poisson, nsim, scores_sorted, score_column):
    """
    Simulation to approximate  P(S >= s_obs | N = n)P(N = n)

    Args:
        mu (float): poisson parameter correpsponding to the sum of the mutation rates
        obs_sum_scores (float): sum of observed DNM scores
        rates (pd.DataFrame): all possible mutations annotated
        nb_mutation_poisson (int): number of mutations to draw
        nsim (int) : number of simulations to perform
        scores_sorted (list) : rates file variants sorted by score
        score_column (str) : CEP scores


    Returns :
        Tuple : (pn, nsim), pn represents the probability of observing a score greater than or equal to the observed score if we select randomly nb_mutation_poisson mutations
        nsim is the number of simulations performed
    """

    # Compute the probability of observing nb_mutation_poisson mutation given the poisson rate (expected number of mutations)
    pndnm = stats.poisson.pmf(nb_mutation_poisson, mu)

    # Scaling the number of simulations to be performed based on the probability of observing nb_mutation_poisson DNMs
    # TODO : see if there is a room for improvement here too
    nsim = max([int(round(nsim * pndnm)), DEFAULT_MIN_NB_SIM])

    # If pndm (see above) is really low there is no point in calculating a combined p-value as it will be epsilon
    if pndnm < 10 ** (-12):
        pscore = 0
        nsim = 0
    # If the top nb_mutation_poisson scores are not enough to reach the observed sum of scores, then there is no possible nb_mutation_poisson combination that
    # will achieve a higher score and the p-value is 0
    elif np.sum(scores_sorted[-nb_mutation_poisson:]) < obs_sum_scores:
        pscore = 0.0
        nsim = 0

    # On the contrary, if the lowest nb_mutation_poisson scores are enough to reach the observed sum of scores, then there is no possible nb_mutation_poisson combination that
    # will achieve a lower score and the p-value is 1
    elif np.sum(scores_sorted[0:nb_mutation_poisson] >= obs_sum_scores):
        pscore = 1.0
        nsim = 0

    # Otherwise, we simulate the cumulated scores for nb_mutation_poisson randomly picked mutations nsim times and calculate the proportion of simulations
    # for which we obtain a score greater than or equal to the observed score
    else:
        s = sim_score(mu, obs_sum_scores, rates, nb_mutation_poisson, nsim, score_column)
        pscore = float(s) / nsim

    # This probability is adjusted by the probability of actually observing nb_mutation_poisson mutations given the poisson rate (expected number of mutations)
    pn = pndnm * pscore

    return (pn, nsim)


def sim_score(mu, obs_sum_scores, rates, nb_mutation_poisson, nsim, score_column):
    """
    Draws nsim times nb_mutation_poisson mutations from set of all possible mutations according to their mutation rate.
    Count how many times their cumulated scores is greater than or equal to the observed sum of scores.

    Args:
        mu (float): poisson parameter correpsponding to the sum of the mutation rates
        obs_sum_scores (float): sum of observed DNM scores
        rates (pd.DataFrame): all possible mutations annotated
        nb_mutation_poisson (int): number of mutations to draw
        nsim (int) : number of simulations to perform
        score_column(str) : column to extract the scores from
    """

    # Precompute probabilities
    probabilities = rates["prob"] / mu
    scores = rates[score_column].values

    def simulate_chunk(chunk_size):
        """Simulate one chunk of scores."""
        simulated_scores = np.sum(np.random.choice(scores, (nb_mutation_poisson, chunk_size), p=probabilities), axis=0)
        return np.sum(simulated_scores >= obs_sum_scores)

    # Split simulations into chunks
    split_sim = 100000
    num_full_chunks = nsim // split_sim
    remaining_simulations = nsim % split_sim

    # Run full chunks in parallel
    ctx = click.get_current_context()
    full_chunk_results = Parallel(n_jobs=ctx.params["jobs"])(
        delayed(simulate_chunk)(split_sim) for _ in range(num_full_chunks)
    )

    # Run remaining simulations (if any)
    remaining_result = simulate_chunk(remaining_simulations) if remaining_simulations > 0 else 0

    # Sum up results
    nb_more_extreme_scores = sum(full_chunk_results) + remaining_result

    return nb_more_extreme_scores


def get_pvalue(rates, obs_sum_scores, nsim, pvalcap, nb_observed_mutations, score_column):
    """
    Calculate the p-value from the enrichment simulation test

    Args:
        rates (DataFrame): all possible SNVs annotated
        obs_sum_scores (float): sum of observed DNM scores
        nsim (int): Number of simulations to perform.
        pvalcap (float): P-value threshold to stop simulations.
        nb_observed_mutations (int) : Number of observed mutations in the current gene
        score_column (str) : CEP scores


    Returns:
        tuple: A tuple containing the p-value, simulation information, and expected score.
    """

    # The poisson rate is the sum of the adjusted mutation rates of all possible mutations in the gene
    mu = rates["prob"].sum()

    # The expected score is the sum of all possible mutations scores weighted by their mutation rates (already adjusted for cohort size)
    exp_sum_scores = np.sum(rates["prob"] * rates[score_column])

    # If the observed score is already lower than the expected one, no need to run the simulation
    if obs_sum_scores < exp_sum_scores:
        ptot = 1
        info = "0|0|True|observed < expected, pvalue set at 1"
        return ptot, info, exp_sum_scores

    # We sort the scores in order to use stopping rules that improve the speed of the simulations
    scores_sorted = np.sort(rates[score_column])

    # Calculate the probability of seeing a similar or more extreme observed score when selecting 0 or 1 mutation
    p0 = calc_p0(mu, obs_sum_scores)
    p1 = calc_p1(mu, obs_sum_scores, rates, score_column)

    # If the gene has a really high number of observed mutations, we increase the default threshold to two times the number of observed mutations
    nb_putative_mutations_sim = max(2 * nb_observed_mutations, DEFAULT_MAX_NB_MUTATIONS_SIM)

    # For long genes with a high number of expected mutations, it is a better strategy
    # to start with a poisson mutation rate close to the expected mutation rate,
    # that way, as most genes will not have an enrichment, we can discard them more quickly
    if mu > 50:
        range_mutation = list(diverging_range(int(mu), 2, nb_putative_mutations_sim))
        sequential = False
    else:
        range_mutation = range(2, nb_putative_mutations_sim)
        sequential = True

    # Running simulations
    ptot = p0 + p1
    nbsim_tot = 0
    info = ""
    for nb_mutation_poisson in range_mutation:

        # Calculate the probability of seeing a similar or more extreme observed gene score | nb_mutation_poisson mutations
        pi, nb_sim = calc_pn(mu, obs_sum_scores, rates, nb_mutation_poisson, nsim, scores_sorted, score_column)
        ptot = ptot + pi
        nbsim_tot = nbsim_tot + nb_sim

        # If we are increasing the number of mutation in sequential order we can
        # stop if the poisson p-value becomes too low from a certain point
        if sequential:
            # Probability of observing more than nb_mutation_poisson event given mu
            picdf = 1 - stats.poisson.cdf(nb_mutation_poisson, mu)

            # If this probability is too low, the resulting p-values will
            if (picdf < 10 ** (-12)) and (nb_mutation_poisson > mu):
                info = "probability of observing >= " + str(nb_mutation_poisson) + " mutations is too small"
                break

        # if p value is over threshold then stop going further
        if ptot > pvalcap:
            info = "pvalue > " + str(pvalcap) + ", stop simulations"
            break

    # Set min p-value to 10^-14
    if ptot == 0:
        ptot = 10 ** (-14)
        info = "p value was 0, set at 10^-14"

    # Add the number of simulations performed to the simulation information
    info = f"{nbsim_tot}|{nb_mutation_poisson}|{sequential}|{info}"

    return (ptot, info, exp_sum_scores)


def diverging_range(median, min_val, max_val):
    """
    Yields a sequence of integers starting from `median` and alternating outward
    in both directions (e.g., median+1, median-1, median+2, median-2, ...),
    within the bounds [min_val, max_val].

    Example:
        list(diverging_range(5, 2, 8))
        → [5, 6, 4, 7, 3, 8, 2]
    """

    yield median  # Start by yielding the median value
    offset = 1  # Initialize the offset from the median

    while True:
        lower = median - offset  # Step downward
        upper = median + offset  # Step upward
        did_yield = False  # Track if any value was yielded this round

        # If upper value is within bounds, yield it
        if upper <= max_val:
            yield upper
            did_yield = True

        # If lower value is within bounds, yield it
        if lower >= min_val:
            yield lower
            did_yield = True

        # Stop the loop if no values were yielded in this iteration
        if not did_yield:
            break

        # Increase the distance from the median for the next iteration
        offset += 1
