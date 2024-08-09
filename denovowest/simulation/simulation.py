#!/usr/bin/env python
# python enrichment.py ../../input/test/new_format/dnm_min.tsv ../../input/test/new_format/all_rates_min.tsv ../../input/weights_ppv_2020_01_17.tab --nmales 10000 --nfemales 10000
import logging
import os
import sys

import click
import pandas as pd
import logging
import numpy as np
import utils


from utils import init_log, log_configuration, CONSEQUENCES_MAPPING
from weights import assign_weights
from probabilities import get_pvalue

RUNTYPE = None


def load_dnm_rates(dnm, rates, column):
    """
    Load DNM and rates files and limit the analysis to genes shared by both files.
    When parallelising DNW on HPC, the rates file is split per gene.

    Args:
        dnm (str): path to observed DNM
        rates (str): path to rates file
        column (str): column that stores scores to use as weights
    """

    logger = logging.getLogger("logger")

    dnm_df = pd.read_csv(dnm, sep="\t", dtype={"chrom": str, "pos": int, column: float}, na_values=[".", "NA"])
    rates_df = pd.read_csv(rates, sep="\t", dtype={"chrom": str, "pos": int, column: float}, na_values=[".", "NA"])

    shared_genes = set(dnm_df.gene_id.unique()) & set(rates_df.gene_id.unique())
    logger.info(
        f"{len(shared_genes)}/{dnm_df.gene_id.nunique()} genes with at least one observed DNM are found in the rates file"
    )

    dnm_df = dnm_df.loc[dnm_df.gene_id.isin(shared_genes)]
    rates_df = rates_df.loc[rates_df.gene_id.isin(shared_genes)]

    return dnm_df, rates_df


def prepare_dnm(dnm_df: pd.DataFrame):
    """
    Filter DNM file to remove functional consequence not handled.
    Assign a higher level consequence to each DNM.

    Args:
        dnmfile (str): DNM file
    """

    dnm_df = filter_on_consequences(dnm_df)
    dnm_df = assign_meta_consequences(dnm_df)

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

    # Filter variants depending on run type : non-synonymous or missense test
    if RUNTYPE == "ns":
        filt = df.consequence.isin(CONSEQUENCES_MAPPING.keys())
    else:
        filt = df.consequence.isin(["missense", "start_lost", "stop_lost"])

    kept_df = df.loc[filt]

    logger.info(f"Before consequence filtering : {df.shape[0]} records")

    discarded_df = df.loc[~filt]
    if not discarded_df.empty:
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

    df.loc[:, "consequence"] = df.consequence.replace(CONSEQUENCES_MAPPING)

    consequence_counts = dict(df.consequence.value_counts())
    for consequence, count in consequence_counts.items():
        logger.info(f"{consequence} : {count}")

    return df


def prepare_rates(rates_df: pd.DataFrame, nmales: int, nfemales: int):
    """
    Load mutation rates file.

    Args:
        ratesfile (str): path to mutation rates file
    """

    rates_df = filter_on_consequences(rates_df)
    rates_df = assign_meta_consequences(rates_df)

    rates_df = compute_expected_number_of_mutations(rates_df, nmales, nfemales)

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
    rates_df.loc[:, "prob"] = rates_df.prob * autosomal_factor

    # Apply an extra correction for the X chromosome
    x_factor_correction = compute_x_factor_correction(nmales, nfemales)
    rates_df.loc[:, "prob"] = np.where(
        rates_df.chrom.isin(["X", "chrX"]), x_factor_correction * rates_df.prob, rates_df.prob
    )

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


def assign_weights_dnm_rates(dnm_df, rates_df, score_column, indel_weights):
    """
    Assign weights to variants in the DNM and rates file.
    The score used for the weights is first min-max transformed

    Args:
        dnm_df (pd.DataFrame): DNM dataframe
        rates_df (pd.DataFrame): rates dataframe that contains all possible SNV
        score_column (str) : the column to use as weights
        indel_weights (pd.DataFrame) : Hardcoded indel weights depending on gene shet score

    """

    # We retrieve the minimum and maximum scores found across the DNM and rates file
    min_score = min(
        rates_df[score_column].min(),
        dnm_df[score_column].min(),
    )
    max_score = max(
        rates_df[score_column].max(),
        dnm_df[score_column].max(),
    )

    # Then we assign the weights
    rates_df = assign_weights(rates_df, score_column, min_score, max_score, RUNTYPE)
    dnm_df = assign_weights(dnm_df, score_column, min_score, max_score, RUNTYPE, indel_weights)

    return dnm_df, rates_df


def run_simulations(dnm_df: pd.DataFrame, rates_df: pd.DataFrame, nsim: int, pvalcap: float):
    """
    For each gene in the DNM file, run nsim simulations and test whether or not this gene is significantly enriched in predictive DNM.

    Args:
        dnm_df (pd.DataFrame): DNM dataframe
        rates_df (pd.DataFrame): rates dataframe that contains all possible SNV
        nsim (int): number of simulations to run
        pvalcap (float): stop simulations if cumulative p-value > pvalcap
    """

    logger = logging.getLogger("logger")

    genes = dnm_df.gene_id.unique()
    results = []
    cpt = 0
    for gene in genes:
        simulation_results = run_simulation(rates_df, dnm_df, gene, nsim, pvalcap)
        if simulation_results:
            results.append(simulation_results)

        cpt += 1
        if cpt % 10 == 0:
            logger.info(f"Processed {cpt}/{len(genes)} genes")

    return results


def run_simulation(rates_df, dnm_df, gene_id, nsim, pvalcap):
    """
    Run nsim simulations and test whether or not gene gene_id is significantly enriched in predictive DNM

    Args:
        rates_df (pd.DataFrame): rates dataframe that contains all possible SNV for the given gene
        dnm_df (pd.DataFrame): DNM dataframe that contains all observed DNM for the given gene
        gene_id (str) : gene identifier
        nsim (int): number of simulations to run
        pvalcap (float): stop simulations if cumulative p-value > pvalcap
    """

    logger = logging.getLogger("logger")

    if gene_id not in rates_df.gene_id.unique():
        logger.debug("could not find " + str(gene_id))
        return

    logger.info(f"Testing {gene_id}")

    # Subset rates file to current gene
    generates = rates_df.loc[rates_df.gene_id == gene_id]

    # Sum the PPV of all observed DNM in the gene
    nb_observed_mutations = dnm_df[dnm_df.gene_id == gene_id].shape[0]
    obs_sum_ppv = dnm_df[dnm_df.gene_id == gene_id].ppv.sum()

    # Run nsim simulations
    pval, info, exp_sum_ppv = get_pvalue(generates, obs_sum_ppv, nsim, pvalcap, nb_observed_mutations)

    # Return the gene id, its expected and observed sum of PPV, the p-value from the enrichment simulation test and some informations about the simulation
    return (gene_id, exp_sum_ppv, obs_sum_ppv, pval, info)


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
@click.argument("column")
@click.option("--nmales", required=True, type=int, help="Number of males individual in your cohort")
@click.option("--nfemales", required=True, type=int, help="Number of females individual in your cohort")
@click.option("--pvalcap", default=1.0, type=float, help="Stop simulations if cumulative p-value > pvalcap")
@click.option("--nsim", type=int, help="Minimum number of simulations for each gene", default=10e7)
@click.option(
    "--runtype",
    help="Run type is either missense test (mis) or non-synonymous (ns)",
    type=click.Choice(["ns", "mis"]),
    default="ns",
    show_default=True,
)
@click.option("--outdir", default="./")
@click.option(
    "--export_weighted_dnmrates",
    is_flag=True,
    help="Export the DNM and rates file that have been annotated with the corresponding weights",
)
def main(dnm, rates, column, nmales, nfemales, pvalcap, nsim, runtype, outdir, export_weighted_dnmrates):
    """
    DeNovoWEST is a simulation-based method to test for a statistically significant enrichment of damaging de novo mutations (DNMs) in individual genes.
    This method scores all classes of variants (e.g. nonsense, missense, splice site) on a unified severity scale based on the empirically-estimated positive predictive value of being pathogenic,
    and incorporates a gene-based weighting derived from the deficit of protein truncating variants in the general population.

    Args:
        dnm (str): DNM file
        rates (str): Per-generation mutation rates file
        column (str): Column to use for the simulation
        nmales (int): Number of males individual in your cohort
        nfemales (int): Number of females individual in your cohort
        pvalcap (float): Stop simulations if cumulative p-value > pvalcap
        nsim (int): Minimum number of simulations for each gene
        runtype (str): Run type is either missense test (mis) or non-synonymous (ns)
        outdir (str): Output directory
        export_weighted_dnmrates (str) : Write the DNM and rates file annotated with weights
    """
    init_log()
    log_configuration(click.get_current_context().params)

    global RUNTYPE
    RUNTYPE = runtype

    # Load DNM and rates files
    dnm_df, rates_df = load_dnm_rates(dnm, rates, column)

    # Prepare DNM and rates file for weight assignation
    dnm_df = prepare_dnm(dnm_df)
    rates_df = prepare_rates(rates_df, nmales, nfemales)

    # Assign weights to variants
    dnm_df, rates_df = assign_weights(dnm_df, rates_df, column)

    # Run simulations
    results = run_simulations(dnm_df, rates_df, nsim, pvalcap)

    # Export results
    export_results(results, outdir)

    # Export the weight annotated dnm and rates files
    if export_weighted_dnmrates:
        export_weighted_files(outdir, dnm_df, rates_df)


if __name__ == "__main__":
    main()
