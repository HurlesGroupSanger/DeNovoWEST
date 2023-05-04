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
from weights import assign_weights


def load_weights(weightsfile: str):
    """_summary_

    Args:
        weightsfile (str): _description_
    """

    weights_df = pd.read_csv(weightsfile, sep="\t")
    return weights_df


def prepare_dnm(dnmfile: str, weights_df: pd.DataFrame):
    """
    Load DNM file.

    Args:
        dnmfile (str): path to DNM file
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
    filt = df.consequence.isin(CONSEQUENCES_MAPPING.keys())
    kept_df = df.loc[filt]

    logger.info(f"Before consequence filtering : {df.shape[0]} DNMs")

    discarded_df = df.loc[~filt]
    if not discarded_df.empty:
        count_discarded = discarded_df["consequence"].value_counts()
        logger.warning(f"{discarded_df.shape[0]}/{df.shape[0]} DNMs were discarded ")
    else:
        logger.info("All DNM have an acceptable consequence.")

    logger.info(f"After consequence filtering : {kept_df.shape[0]} DNMs")

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
    rates_df = compute_expected_number_of_mutations(rates_df, nmales, nfemales)
    rates_df = assign_weights(rates_df, weights_df)

    return rates_df


def compute_expected_number_of_mutations(rates_df: pd.DataFrame, nmales: int, nfemales: int):
    """
    The rates file associates each possible variant with a mutation probability.
    Args:
        rates_df (pd.DataFrame): _description_
        nmales (int): _description_
        nfemales (int): _description_
    """

    # Compute the expected number of mutations as the product of the mutation rate and the number of individuals
    autosomal_factor = 2 * (nmales + nfemales)
    rates_df.prob = rates_df.prob * autosomal_factor

    # Correct for X chromosme
    x_factor_correction = compute_x_factor_correction(nmales, nfemales)
    rates_df.prob = np.where(rates_df.chrom.isin(["X", "chrX"]), x_factor_correction * rates_df.prob, rates_df.prob)

    return rates_df


def compute_x_factor_correction(nmales: int, nfemales: int):
    """# TODO : need some more explanation

    Args:
        nmales (int): _description_
        nfemales (int): _description_
    """

    autosomal_factor = 2 * (nmales + nfemales)

    female_transmit = nmales + nfemales
    male_transmit = nfemales

    # get scaling factors using the alpha from the most recent SFHS (Scottish
    # Family Health Study) phased de novo data.
    alpha = 3.4
    male_k = 2 / (1 + (1 / alpha))
    female_k = 2 / (1 + alpha)

    x_factor = ((male_transmit * male_k) + (female_transmit * female_k)) / autosomal_factor

    return x_factor


def correct_for_x_chrom(rates, male_n, female_n):
    """
    correct for X chromosome
    """
    autosomal = 2 * (male_n + female_n)
    female_transmit = male_n + female_n
    male_transmit = female_n

    # get scaling factors using the alpha from the most recent SFHS (Scottish
    # Family Health Study) phased de novo data.
    alpha = 3.4
    male_k = 2 / (1 + (1 / alpha))
    female_k = 2 / (1 + alpha)

    # correct the non-PAR chrX genes for fewer transmissions and lower rate
    # (dependent on alpha)
    chrX = rates["chrom"].isin(["X", "chrX"])
    x_factor = ((male_transmit * male_k) + (female_transmit * female_k)) / autosomal
    x_factor = pd.Series([x_factor] * len(chrX), index=rates.index)
    x_factor[~chrX] = 1

    rates["prob"] *= x_factor

    return rates


def run_simulation(dnm_df: pd.DataFrame, rates_df: pd.DataFrame):
    """_summary_

    Args:
        dnm_df (pd.DataFrame): _description_
        rates_df (pd.DataFrame): _description_
    """

    print("toto")


def export_results(results_df: pd.DataFrame, output: str):
    """_summary_

    Args:
        results_df (pd.DataFrame): _description_
        output (str): _description_
    """

    df = pd.DataFrame.from_records(results_df, columns=["symbol", "hgnc_id", "expected", "observed", "p-value", "info"])
    df.to_csv(output, sep="\t", index=False)


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
@click.option("--output", default="enrichment_results.tsv")
def main(dnm, rates, weights, nmales, nfemales, pvalcap, nsim, output):
    """#TODO

    Args:
        dnm (str): Path to denovo mutation file
        rates (str): Path to mutation rates file
        weights (str): Path to weights file
        nmales (int): Number of males individual in your cohort
        nfemales (int): Number of females individual in your cohort
        pvalcap (float): Stop simulations if cumulative p-value > pvalcap
        nsim (int): Minimum number of simulations for each gene
        output (str): _description_
    """

    init_log()
    log_configuration(click.get_current_context().params)

    weights_df = load_weights(weights)
    dnm_df = prepare_dnm(dnm, weights_df)
    rates_df = prepare_rates(rates, weights_df, nmales, nfemales)

    results_df = run_simulation(dnm_df, rates_df, nsim, pvalcap)

    export_results(results_df, output)


if __name__ == "__main__":
    main()
