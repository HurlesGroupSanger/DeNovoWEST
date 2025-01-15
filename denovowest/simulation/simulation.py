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
from scores import prepare_scores
from probabilities import get_pvalue

RUNTYPE = None


def load_dnm_rates(dnm, rates, column):
    """
    Load DNM and rates files and limit the analysis to genes shared by both files.
    When parallelising DNW on HPC, the rates file is split per gene.

    Args:
        dnm (str): path to observed DNM
        rates (str): path to rates file
        column (str): column that stores variant scores
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


def prepare_dnm_rates(
    dnm_df: pd.DataFrame, rates_df: pd.DataFrame, column: str, nmales: int, nfemales: int, impute_missing: bool
):
    """
    Prepare DNM and rates file for the simulation

    Args:
        dnm_df (pd.DataFrame): DNM dataframe
        rates_df (pd.DataFrame): mutation rates dataframe
        column (str): column that stores variant scores
        nmales (int): number of males in the cohort
        nfemales (int): number of females in the cohort
        impute_missing (bool): whether to impute the variant with missing scores or not
    """
    # Filter on variant consequence and calculate cohort based expected mutation rates
    dnm_df = prepare_dnm(dnm_df)
    rates_df = prepare_rates(rates_df, nmales, nfemales)

    # Prepare scores (indel scores, imputation missing scores...)
    dnm_df, rates_df = prepare_scores(dnm_df, rates_df, column, RUNTYPE, impute_missing)

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


def run_simulations(dnm_df: pd.DataFrame, rates_df: pd.DataFrame, nsim: int, pvalcap: float, score_column: str):
    """
    For each gene in the DNM file, run nsim simulations and test whether or not this gene is significantly enriched in predictive DNM.

    Args:
        dnm_df (pd.DataFrame): DNM dataframe
        rates_df (pd.DataFrame): rates dataframe that contains all possible SNV
        nsim (int): number of simulations to run
        pvalcap (float): stop simulations if cumulative p-value > pvalcap
        score_column (str) : CEP scores
    """

    logger = logging.getLogger("logger")

    genes = dnm_df.gene_id.unique()
    results = []
    cpt = 0
    for gene in genes:
        simulation_results = run_simulation(rates_df, dnm_df, gene, nsim, pvalcap, score_column)
        if simulation_results:
            results.append(simulation_results)

        cpt += 1
        if cpt % 10 == 0:
            logger.info(f"Processed {cpt}/{len(genes)} genes")

    return results


def run_simulation(rates_df, dnm_df, gene_id, nsim, pvalcap, score_column):
    """
    Run nsim simulations and test whether or not gene gene_id is significantly enriched in predictive DNM

    Args:
        rates_df (pd.DataFrame): rates dataframe that contains all possible SNV for the given gene
        dnm_df (pd.DataFrame): DNM dataframe that contains all observed DNM for the given gene
        gene_id (str) : gene identifier
        nsim (int): number of simulations to run
        pvalcap (float): stop simulations if cumulative p-value > pvalcap
        score_column (str) : CEP scores
    """

    logger = logging.getLogger("logger")

    if gene_id not in rates_df.gene_id.unique():
        logger.debug("could not find " + str(gene_id))
        return

    logger.info(f"Testing {gene_id}")

    # Subset rates file to current gene
    generates = rates_df.loc[rates_df.gene_id == gene_id]

    # Sum the scores of all observed DNM in the gene
    nb_observed_mutations = dnm_df[dnm_df.gene_id == gene_id].shape[0]
    obs_sum_scores = dnm_df[dnm_df.gene_id == gene_id][score_column].sum()

    # Run nsim simulations
    pval, info, exp_sum_scores = get_pvalue(
        generates, obs_sum_scores, nsim, pvalcap, nb_observed_mutations, score_column
    )

    # Return the gene id, its expected and observed sum of scores, the p-value from the enrichment simulation test and some informations about the simulation
    return (gene_id, exp_sum_scores, obs_sum_scores, pval, info)


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


@click.command()
@click.argument("dnm")
@click.argument("rates")
@click.argument("column")
@click.option("--nmales", required=True, type=int, help="Number of males individual in your cohort")
@click.option("--nfemales", required=True, type=int, help="Number of females individual in your cohort")
@click.option("--pvalcap", default=1, type=float, help="Stop simulations if cumulative p-value > pvalcap")
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
    "--impute-missing-scores",
    is_flag=True,
    help="Impute missing scores by taking the median score for similar variants (i.e. same consequence) in the gene",
)
def main(dnm, rates, column, nmales, nfemales, pvalcap, nsim, runtype, outdir, impute_missing_scores):
    """
    DeNovoWEST is a simulation-based method that tests for DNM enrichment in each gene separately.
    It uses computational effect predictors (CEP) scores to reflect the pathogenicity probability for each variant.


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
        impute_missing_scores (bool) : Impute missing scores by taking the median score for similar variants (i.e. same consequence) in the gene
    """
    init_log()
    log_configuration(click.get_current_context().params)

    global RUNTYPE
    RUNTYPE = runtype

    # Load DNM and rates files
    dnm_df, rates_df = load_dnm_rates(dnm, rates, column)

    # Prepare DNM and rates file for simulation
    dnm_df, rates_df = prepare_dnm_rates(dnm_df, rates_df, column, nmales, nfemales, impute_missing_scores)

    # Run simulations
    results = run_simulations(dnm_df, rates_df, nsim, pvalcap, column)

    # Export results
    export_results(results, outdir)


if __name__ == "__main__":
    main()
