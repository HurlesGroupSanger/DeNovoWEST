#!/usr/bin/env python
# python enrichment.py ../../input/test/new_format/dnm_min.tsv ../../input/test/new_format/all_rates_min.tsv ../../input/weights_ppv_2020_01_17.tab --nmales 10000 --nfemales 10000
import logging
import os
import sys

import click
import pandas as pd
import logging
import numpy as np


from utils import init_log, log_configuration, CONSEQUENCES_MAPPING
from weights import assign_weights, get_indel_weights
from probabilities import get_pvalue


def load_weights(weightsfile: str):
    """
    Load weights file containing the positive predictive value (PPV) of each variant type.

    Args:
        weightsfile (str): weights file path
    """

    weights_df = pd.read_csv(weightsfile, sep="\t")
    return weights_df


def prepare_dnm(dnmfile: str, weights_df: pd.DataFrame):
    """
    Filter DNM file to remove functional consequence not handled.
    Assign a higher level consequence to each DNM.
    Assign weights to each DNM.

    Args:
        dnmfile (str): DNM file
        weights_df(pd.DataFrame): weights dataframe
    """

    dnm_df = pd.read_csv(dnmfile, sep="\t", dtype={"chrom": str, "pos": int, "score": float})

    dnm_df = filter_on_consequences(dnm_df)
    dnm_df = assign_meta_consequences(dnm_df)

    dnm_df = assign_weights(dnm_df, weights_df)

    return dnm_df


def filter_on_consequences(df: pd.DataFrame):
    """
    Filter all variants with a consequence not found in CONSEQUENCES_MAPPING

    Args:
        df (pd.DataFrame): variant table (rates or dnm) having a consequence column
    """

    logger = logging.getLogger("logger")

    # TODO : Replace with better consequence extraction
    # When bcftools report several consequences separated by the character "&", we extract the first one
    df.consequence = [csq.split("&")[0] if isinstance(csq, str) else csq for csq in list(df.consequence)]

    filt = df.consequence.isin(CONSEQUENCES_MAPPING.keys())
    kept_df = df.loc[filt]

    logger.info(f"Before consequence filtering : {df.shape[0]} records")

    discarded_df = df.loc[~filt]
    if not discarded_df.empty:
        count_discarded = discarded_df["consequence"].value_counts()
        logger.warning(f"{discarded_df.shape[0]}/{df.shape[0]} records were discarded")
    else:
        logger.info("All records have an acceptable consequence.")

    logger.info(f"After consequence filtering : {kept_df.shape[0]} records")

    return kept_df


def assign_meta_consequences(df: pd.DataFrame):
    """
    Assign a higher level consequence to each variant

    Args:
        df (pd.DataFrame): variant table (rates or dnm) having a consequence column
    """

    logger = logging.getLogger("logger")

    df.consequence = df.consequence.replace(CONSEQUENCES_MAPPING)

    consequence_counts = dict(df.consequence.value_counts())
    for consequence, count in consequence_counts.items():
        logger.info(f"{consequence} : {count}")

    return df


def prepare_rates(ratesfile: str, weights_df: pd.DataFrame, nmales: int, nfemales: int):
    """
    Load mutation rates file.

    Args:
        ratesfile (str): path to mutation rates file
    """

    rates_df = pd.read_csv(ratesfile, sep="\t", dtype={"chrom": str, "pos": int, "score": float})

    rates_df = filter_on_consequences(rates_df)
    rates_df = assign_meta_consequences(rates_df)

    rates_df = compute_expected_number_of_mutations(rates_df, nmales, nfemales)
    rates_df = assign_weights(rates_df, weights_df)

    return rates_df


def compute_expected_number_of_mutations(rates_df: pd.DataFrame, nmales: int, nfemales: int):
    """
    Adjust the expected number of mutations to the number of individuals in the cohort.

    Args:
        rates_df (pd.DataFrame): mutation rates dataframe
        nmales (int): number of males in the cohort
        nfemales (int): number of females in the cohort
    """

    # Compute the expected number of mutations as the product of the mutation rate and the number of individuals
    autosomal_factor = 2 * (nmales + nfemales)
    rates_df.prob = rates_df.prob * autosomal_factor

    # Apply an extra correction for the X chromosome
    x_factor_correction = compute_x_factor_correction(nmales, nfemales)
    rates_df.prob = np.where(rates_df.chrom.isin(["X", "chrX"]), x_factor_correction * rates_df.prob, rates_df.prob)

    return rates_df


def compute_x_factor_correction(nmales: int, nfemales: int):
    """
    Expected number of mutations on the X chromosome need to be adjusted to the number of male and female individuals.
    Scaling factors are computed using the alpha from the most recent SFHS (Scottish Family Health Study) phased de novo data.
    Correct the non-PAR chrX genes for fewer transmissions and lower rate (depends on alpha)

    Args:
        nmales (int): number of males in the cohort
        nfemales (int): number of females in the cohort
    """

    autosomal_factor = 2 * (nmales + nfemales)

    female_transmit = nmales + nfemales
    male_transmit = nfemales

    alpha = 3.4
    male_k = 2 / (1 + (1 / alpha))
    female_k = 2 / (1 + alpha)

    x_factor = ((male_transmit * male_k) + (female_transmit * female_k)) / autosomal_factor

    return x_factor


def run_simulations(
    dnm_df: pd.DataFrame, rates_df: pd.DataFrame, nsim: int, indel_weights: pd.DataFrame, pvalcap: float
):
    """
    For each gene in the DNM file, run nsim simulations and test whether or not this gene is significantly enriched in predictive DNM.

    Args:
        dnm_df (pd.DataFrame): DNM dataframe
        rates_df (pd.DataFrame): rates dataframe that contains all possible SNV
        nsim (int): number of simulations to run
        indel_weights (pd.DataFrame): indel weights
        pvalcap (float): stop simulations if cumulative p-value > pvalcap
    """

    logger = logging.getLogger("logger")

    genes = dnm_df.gene_id.unique()
    results = []
    cpt = 0
    for gene in genes:
        simulation_results = run_simulation(rates_df, dnm_df, gene, nsim, indel_weights, pvalcap)
        if simulation_results:
            results.append(simulation_results)

        cpt += 1
        if cpt % 10 == 0:
            logger.info(f"Processed {cpt}/{len(genes)} genes")

    return results


def run_simulation(rates_df, dnm_df, gene_id, nsim, indel_weights, pvalcap):
    """
    Run nsim simulations and test whether or not gene gene_id is significantly enriched in predictive DNM

    Args:
        rates_df (pd.DataFrame): rates dataframe that contains all possible SNV for the given gene
        dnm_df (pd.DataFrame): DNM dataframe that contains all observed DNM for the given gene
        gene_id (str) : gene identifier
        nsim (int): number of simulations to run
        indel_weights (pd.DataFrame): indels weights
        pvalcap (float): stop simulations if cumulative p-value > pvalcap
    """

    logger = logging.getLogger("logger")

    if gene_id not in rates_df.gene_id.unique():
        logger.debug("could not find " + str(gene_id))
        return

    logger.info(f"Testing {gene_id}")

    generates = rates_df.loc[rates_df.gene_id == gene_id]
    # Add the gene specific inframe and frameshift rates to the rates dataframe
    generates = get_indel_rates(generates, indel_weights)

    # Sum the PPV of all observed DNM in the gene
    obs_sum_ppv = dnm_df[dnm_df.gene_id == gene_id].ppv.sum()

    # Run nsim simulations
    pval, info, exp_sum_ppv = get_pvalue(generates, obs_sum_ppv, nsim, pvalcap)

    # Return the gene id, its expected and observed sum of PPV, the p-value from the enrichment simulation test and some informations about the simulation
    return (gene_id, exp_sum_ppv, obs_sum_ppv, pval, info)


def get_indel_rates(generates, indel_weights):
    """
    Infer the gene specific indel rates from the rate of missense and nonsense mutations.
    Rates file contain only SNP rates, so we need to infer the indel rates from the SNP rates.

    Args:
        generates (pd.DataFrame): rates for the current gene
        indel_weights (pd.DataFrame): indel weights

    Returns:
        pd.DataFrame: rates for the current gene updated with indel rates
    """

    logger = logging.getLogger("logger")

    # Get the overall probability of missense and nonsense mutation across the gene
    missense_rate = generates[generates.consequence == "missense"].prob.sum()
    nonsense_rate = generates[generates.consequence == "nonsense"].prob.sum()

    # TODO : Where does those factors come from ?
    # Set them up as parameters or constants
    inframe_rate = missense_rate * 0.03
    frameshift_rate = nonsense_rate * 1.3

    # Get the weights associated to frameshift and inframe indels
    shethigh = generates.shethigh.iloc[0]
    try:
        frameshift_weight = float(
            indel_weights.loc[(indel_weights.consequence == "frameshift") & (indel_weights.shethigh == shethigh)].ppv
        )
    except TypeError:
        logger.warning("No frameshift ppv found in weight file")
        frameshift_weight = np.NaN

    try:
        inframe_weight = float(
            indel_weights.loc[(indel_weights.consequence == "inframe") & (indel_weights.shethigh == shethigh)].ppv
        )
    except TypeError:
        logger.warning("No inframe ppv found in weight file")
        inframe_weight = np.NaN

    # Add the weights to the rates dataframe
    indelrates = pd.DataFrame(
        [
            ["inframe", inframe_rate, inframe_weight, False, shethigh],
            ["frameshift", frameshift_rate, frameshift_weight, False, shethigh],
        ],
        columns=["consequence", "prob", "ppv", "constrained", "shethigh"],
    )
    generates = pd.concat([generates, indelrates])
    return generates


def export_results(results: list, outdir: str):
    """
    Write enrichment results

    Args:
        results (list): list of per-gene enrichment simulation results
        outdir (str): output directory
    """

    os.makedirs(outdir, exist_ok=True)

    df = pd.DataFrame.from_records(results, columns=["symbol", "expected", "observed", "p-value", "info"])
    df.to_csv(f"{outdir}/enrichment_results.tsv", sep="\t", index=False)


def export_weighted_files(outdir: str, dnm_df: pd.DataFrame, rates_df: pd.DataFrame):
    """
    Export the DNM and rates file that have been annotated with the corresponding weights

    Args:
        outdir (str): output directory
        dnm_df (pd.DataFrame): DNM table with DNW weights
        rates_df (pd.DataFrame): Rates table with DNW weights

    """

    dnm_df.to_csv(f"{outdir}/weighted_DNM.tsv", sep="\t", index=False)
    rates_df.to_csv(f"{outdir}/weighted_rates.tsv", sep="\t", index=False)


@click.command()
@click.argument("dnm")
@click.argument("rates")
@click.argument("weights")
@click.option("--nmales", required=True, type=int, help="Number of males individual in your cohort")
@click.option("--nfemales", required=True, type=int, help="Number of females individual in your cohort")
@click.option(
    "--pvalcap", default=1.0, type=float, help="Stop simulations if cumulative p-value > pvalcap"
)  # TODO more details
@click.option("--nsim", default=10e9, type=int, help="Minimum number of simulations for each gene ")
@click.option("--outdir", default="./")
@click.option(
    "--export_weighted_dnmrates",
    is_flag=True,
    help="Export the DNM and rates file that have been annotated with the corresponding weights",
)
def main(dnm, rates, weights, nmales, nfemales, pvalcap, nsim, outdir, export_weighted_dnmrates):
    """
    DeNovoWEST is a simulation-based method to test for a statistically significant enrichment of damaging de novo mutations (DNMs) in individual genes.
    This method scores all classes of variants (e.g. nonsense, missense, splice site) on a unified severity scale based on the empirically-estimated positive predictive value of being pathogenic,
    and incorporates a gene-based weighting derived from the deficit of protein truncating variants in the general population.

    Args:
        dnm (str): DNM file
        rates (str): Per-generation mutation rates file
        weights (str): Weights file
        nmales (int): Number of males individual in your cohort
        nfemales (int): Number of females individual in your cohort
        pvalcap (float): Stop simulations if cumulative p-value > pvalcap
        nsim (int): Minimum number of simulations for each gene
        outdir (str): Output directory
        export_weighted_dnmrates (str) : Write the DNM and rates file annotated with weights
    """
    init_log()
    log_configuration(click.get_current_context().params)

    # Load weights
    weights_df = load_weights(weights)
    indel_weights = get_indel_weights(weights_df)

    # Assign weights to DNM and rates
    dnm_df = prepare_dnm(dnm, weights_df)
    rates_df = prepare_rates(rates, weights_df, nmales, nfemales)

    # Run simulations
    results = run_simulations(dnm_df, rates_df, nsim, indel_weights, pvalcap)

    # Export results
    export_results(results, outdir)

    # Export the weight annotated dnm and rates files
    if export_weighted_dnmrates:
        export_weighted_files(outdir, dnm_df, rates_df)


if __name__ == "__main__":
    main()
