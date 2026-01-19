#!/usr/bin/env python
# python enrichment.py ../../input/test/new_format/dnm_min.tsv ../../input/test/new_format/all_rates_min.tsv ../../input/weights_ppv_2020_01_17.tab --nmales 10000 --nfemales 10000
import json
import logging
import math
import os
import time
from dataclasses import asdict

import click
import numpy as np
import pandas as pd
from config import Config

from denovowest.simulation.probabilities import get_pvalue
from denovowest.simulation.scores import prepare_scores
from denovowest.utils.log import init_log, set_plain_log, set_regular_log
from denovowest.utils.params import CONSEQUENCES_MAPPING, CONSEQUENCES_SEVERITIES


def load_dnm_rates(dnm, rates, column, gene_list):
    """
    Load DNM and rates files and limit the analysis to genes shared by both files.
    When parallelising DNW on HPC, the rates file is split per gene.

    Args:
        dnm (str): path to observed DNM
        rates (str): path to rates file
        column (str): column that stores variant scores
        gene_list (str) : list of genes to consider
    """

    logger = logging.getLogger("logger")

    dnm_df = pd.read_csv(dnm, sep="\t", dtype={"chrom": str, "pos": int, column: float}, na_values=[".", "NA"])
    rates_df = pd.read_csv(rates, sep="\t", dtype={"chrom": str, "pos": int, column: float}, na_values=[".", "NA"])

    if gene_list:
        with open(gene_list, "r") as f:
            genes = [x.strip() for x in f.readlines()]

        dnm_df = dnm_df.loc[dnm_df.gene_id.isin(genes)]
        rates_df = rates_df.loc[rates_df.gene_id.isin(genes)]

    shared_genes = set(dnm_df.gene_id.unique()) & set(rates_df.gene_id.unique())

    set_plain_log()
    logger.info("=" * 50)
    logger.info("[DNM/RATES CONSISTENCY]")
    logger.info("=" * 50)
    set_regular_log()
    logger.info(
        f"{len(shared_genes)}/{rates_df.gene_id.nunique()} genes from the rates file have at least one observation in the DNM file"
    )

    dnm_df = dnm_df.loc[dnm_df.gene_id.isin(shared_genes)].copy()
    rates_df = rates_df.loc[rates_df.gene_id.isin(shared_genes)].copy()

    return dnm_df, rates_df


def prepare_dnm_rates(
    dnm_df: pd.DataFrame, rates_df: pd.DataFrame, score_column: str, nmales: int, nfemales: int, cfg: Config
):
    """
    Prepare DNM and rates file for the simulation

    Args:
        dnm_df (pd.DataFrame): DNM dataframe
        rates_df (pd.DataFrame): mutation rates dataframe
        score_column (str): column that stores variant scores
        nmales (int): number of males in the cohort
        nfemales (int): number of females in the cohort
        cfg (Config): configuration object that stores script parameters
    """
    # Filter on variant consequence and calculate cohort based expected mutation rates
    dnm_df = prepare_dnm(dnm_df, cfg)
    rates_df = prepare_rates(rates_df, nmales, nfemales, cfg)

    # Prepare scores (indel scores, imputation missing scores...)
    dnm_df, rates_df = prepare_scores(dnm_df, rates_df, score_column, cfg)

    return dnm_df, rates_df


def prepare_dnm(dnm_df: pd.DataFrame, cfg: Config):
    """
    Filter DNM file to remove functional consequence not handled.
    Assign a higher level consequence to each DNM.

    Args:
        dnm_df (pd.DataFrame): DNM dataframe
        cfg (Config): configuration object that stores script parameters
    """

    dnm_df = filter_on_consequences(dnm_df, "dnm", cfg)
    dnm_df = assign_meta_consequences(dnm_df)

    return dnm_df


def filter_on_consequences(df: pd.DataFrame, mode: str, cfg: Config):
    """
    Filter all variants with a consequence not found in CONSEQUENCES_MAPPING

    Args:
        df (pd.DataFrame): variant table (rates or dnm) having a consequence column
        mode(str) : dnm or rates
        cfg (Config): configuration object that stores script parameters
    """

    logger = logging.getLogger("logger")

    set_plain_log()
    logger.info("=" * 50)
    logger.info(f"[{mode.upper()}]")
    logger.info("=" * 50)
    set_regular_log()

    # Extract the worst consequence but also keep the original consequence called by bcftoolscsq
    if mode == "dnm":
        df["original_consequence"] = df.consequence
    df.consequence = [extract_worst_consequence(csq) if isinstance(csq, str) else csq for csq in list(df.consequence)]

    # Filter variants depending on run type : non-synonymous, missense or synonymous test
    if cfg.runtype == "ns":
        filt = df.consequence.isin(CONSEQUENCES_MAPPING.keys())
    elif cfg.runtype == "mis":
        filt = df.consequence.isin(["missense", "start_lost", "stop_lost"])
    else:  # syn
        filt = df.consequence.isin(["synonymous"])

    kept_df = df.loc[filt].copy()

    logger.info(f"{mode.upper()} - Before consequence filtering : {df.shape[0]} records")

    discarded_df = df.loc[~filt]
    if not discarded_df.empty:
        logger.warning(f"{discarded_df.shape[0]}/{df.shape[0]} records were discarded")
    else:
        logger.info("All records have an acceptable consequence.")

    logger.info(f"{mode.upper()} - After consequence filtering : {kept_df.shape[0]} records")

    return kept_df


def extract_worst_consequence(csq):
    """
    Bcftoolscsq sometimes return multiple consequences
    separated by an "&". We extract the worst one.

    Args:
        csq (str): consequences separated by "&"
    """

    list_csq = csq.split("&")
    worst_csq_value = 1
    worst_csq = ""
    for cur_csq in list_csq:
        if CONSEQUENCES_SEVERITIES[cur_csq] > worst_csq_value:
            worst_csq_value = CONSEQUENCES_SEVERITIES[cur_csq]
            worst_csq = cur_csq

    return worst_csq


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


def prepare_rates(rates_df: pd.DataFrame, nmales: int, nfemales: int, cfg: Config):
    """
    Prepare rates file for the simulation

    Args:
        rates_df (pd.DataFrame): mutation rates dataframe
        nmales (int): number of males in the cohort
        nfemales (int): number of females in the cohort
        cfg (Config): configuration object that stores script parameters

    """

    rates_df = filter_on_consequences(rates_df, "rates", cfg)
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
    Scaling factors are computed using the alpha (the male-to-female germline mutation rate ratio)
    from the most recent SFHS (Scottish Family Health Study) phased de novo data.
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


def run_simulations(dnm_df: pd.DataFrame, rates_df: pd.DataFrame, score_column: str, cfg: Config):
    """
    For each gene in the DNM file, run nsim simulations and test whether or not this gene is significantly enriched in predictive DNM.

    Args:
        dnm_df (pd.DataFrame): DNM dataframe
        rates_df (pd.DataFrame): rates dataframe that contains all possible SNV
        score_column (str) : CEP scores
        cfg (Config): configuration object that stores script parameters
    """

    logger = logging.getLogger("logger")

    set_plain_log()
    logger.info("=" * 50)
    logger.info("[SIMULATION]")
    logger.info("=" * 50)
    set_regular_log()

    genes = dnm_df.gene_id.unique()
    results = []
    logs = dict()
    cpt = 0
    for gene in genes:
        simulation_results, simulation_logs = run_simulation(rates_df, dnm_df, gene, score_column, cfg)
        if simulation_results:
            results.append(simulation_results)
            logs[gene] = simulation_logs

        cpt += 1
        if cpt % 10 == 0:
            logger.info(f"Processed {cpt}/{len(genes)} genes")

    return results, logs


def run_simulation(rates_df, dnm_df, gene_id, score_column, cfg):
    """
    Run nsim simulations and test whether or not gene gene_id is significantly enriched in predictive DNM

    Args:
        rates_df (pd.DataFrame): rates dataframe that contains all possible SNV for the given gene
        dnm_df (pd.DataFrame): DNM dataframe that contains all observed DNM for the given gene
        gene_id (str) : gene identifier
        score_column (str) : CEP scores
        cfg (Config): configuration object that stores script parameters
    """

    logger = logging.getLogger("logger")

    simulation_logs = dict()

    if gene_id not in rates_df.gene_id.unique():
        logger.debug("could not find " + str(gene_id))
        return

    logger.info(f"Testing {gene_id}")
    start_time = time.time()

    # Subset rates file to current gene
    generates = rates_df.loc[rates_df.gene_id == gene_id]

    # Sum the scores of all observed DNM in the gene
    nb_observed_mutations = dnm_df[dnm_df.gene_id == gene_id].shape[0]
    obs_sum_scores = dnm_df[dnm_df.gene_id == gene_id][score_column].sum()

    simulation_logs["nb_observed_mutations"] = nb_observed_mutations
    simulation_logs["obs_sum_scores"] = obs_sum_scores

    # Run nsim simulations
    results, simulation_logs = get_pvalue(
        generates, obs_sum_scores, nb_observed_mutations, score_column, cfg, simulation_logs
    )

    # Store how long the simulation took for each gene
    end_time = time.time()
    wall_time = end_time - start_time
    info = f"{results[1]}|{wall_time:.6f}"

    # Return the gene id, its expected and observed sum of scores, the p-value from the enrichment simulation test and some informations about the simulation
    simulation_results = (gene_id, results[2], obs_sum_scores, results[0], info)
    return simulation_results, simulation_logs


def export_results(results: list, outdir: str, outfile: str):
    """
    Write enrichment results

    Args:
        results (list): list of per-gene enrichment simulation results
        outdir (str): output directory
        outfile (str): enrichment results file
    """

    logger = logging.getLogger("logger")

    set_plain_log()
    logger.info("=" * 50)
    logger.info("[EXPORT]")
    logger.info("=" * 50)
    set_regular_log()

    os.makedirs(outdir, exist_ok=True)

    df = pd.DataFrame.from_records(results, columns=["symbol", "expected", "observed", "p-value", "info"])
    df["nb_sim"] = [x.split("|")[0] for x in df["info"]]
    df["nb_mutations_stop"] = [x.split("|")[1] for x in df["info"]]
    df["sequential_simulation"] = [x.split("|")[2] for x in df["info"]]
    df["log"] = [x.split("|")[3] for x in df["info"]]
    df["wall_time"] = [x.split("|")[4] for x in df["info"]]

    df.drop("info", axis=1, inplace=True)

    df.to_csv(f"{outdir}/{outfile}", sep="\t", index=False)

    logger.info(f"Simulation results exported to {outdir}/{outfile}")


def export_logs(logs, outdir):

    class NpEncoder(json.JSONEncoder):
        def default(self, obj):
            # Handle numpy integer
            if isinstance(obj, np.integer):
                return int(obj)
            # Handle numpy floating
            elif isinstance(obj, np.floating):
                if math.isnan(obj):  # NaN → null
                    return None
                return float(obj)
            # Handle numpy arrays
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            return super().default(obj)

    logger = logging.getLogger("logger")

    with open(f"{outdir}/simulation_logs.json", "w") as f:
        json.dump(logs, f, indent=4, cls=NpEncoder)

    logger.info(f"Simulation logs exported to {outdir}/simulation_logs.json")


def log_configuration(conf):
    """

    Log the configuration the simulation script was run with

    Args:
        conf (dict): script arguments
    """

    logger = logging.getLogger("logger")

    set_plain_log()
    logger.info("=" * 50)
    logger.info("[CONFIGURATION]")
    logger.info("=" * 50)
    set_regular_log()

    for key, value in sorted(conf.items()):
        logger.info(f"{key} : {value}")


@click.command()
@click.argument("dnm")
@click.argument("rates")
@click.argument("score_column")

# Required cohort parameters
@click.option("--nmales", required=True, type=int, help="Number of male individuals in the cohort")
@click.option("--nfemales", required=True, type=int, help="Number of female individuals in the cohort")
@click.option(
    "--gene-list", default=Config().gene_list, show_default=True, help="Restrict analysis to genes in the provided list"
)

# Variant processing
@click.option(
    "--impute-missing-scores",
    is_flag=True,
    default=Config().impute_missing_scores,
    show_default=True,
    help="Impute missing variant scores using the median of similar variants in the same gene",
)
@click.option(
    "--runtype",
    type=click.Choice(["ns", "mis", "syn"]),
    default=Config().runtype,
    show_default=True,
    help="Run type: 'mis' for missense, 'ns' for non-synonymous, 'syn' for synonymous",
)

# Inframe scoring
@click.option(
    "--inframe-missense-ratio",
    type=float,
    default=Config().inframe_missense_ratio,
    show_default=True,
    help="Observed inframe to missense ratio",
)
@click.option(
    "--inframe-score-quantile",
    type=float,
    default=Config().inframe_score_quantile,
    show_default=True,
    help="Inframes are assigned the score at this quantile of the missense score distribution in the gene (0.5 = median)",
)

# Frameshift scoring
@click.option(
    "--frameshift-nonsense-ratio",
    type=float,
    default=Config().frameshift_nonsense_ratio,
    show_default=True,
    help="Observed frameshift to nonsense ratio",
)
@click.option(
    "--frameshift-score-quantile",
    type=float,
    default=Config().frameshift_score_quantile,
    show_default=True,
    help="Frameshifts are assigned the score at this quantile of the nonsense score distribution in the gene (0.6 = 60th percentile)",
)

# Simulation
@click.option(
    "--nsim", type=int, default=Config().nsim, show_default=True, help="Minimum number of simulations per gene"
)
@click.option(
    "--pvalcap",
    type=float,
    default=Config().pvalcap,
    show_default=True,
    help="Stop simulations when cumulative p-value exceeds this threshold",
)
@click.option(
    "--jobs", type=int, default=Config().jobs, show_default=True, help="Number of cores to use during simulations"
)

# Output / debug
@click.option(
    "--debug",
    is_flag=True,
    default=Config().debug,
    show_default=True,
    help="Log detailed simulation information for each gene. Set a seed for reproducibility.",
)
@click.option("--outdir", default=Config().outdir, show_default=True, help="Output directory")
@click.option("--outfile", default=Config().outfile, show_default=True, help="Name of the enrichment results file")
def main(
    dnm,
    rates,
    score_column,
    nmales,
    nfemales,
    gene_list,
    impute_missing_scores,
    runtype,
    inframe_missense_ratio,
    inframe_score_quantile,
    frameshift_nonsense_ratio,
    frameshift_score_quantile,
    nsim,
    pvalcap,
    jobs,
    debug,
    outdir,
    outfile,
):
    """
    DeNovoWEST performs gene-level simulations to test for de novo mutation (DNM) enrichment, incorporating
    computational effect predictor (CEP) scores to account for predicted variant pathogenicity.
    """

    # Build config object to pass around
    cfg = Config(
        nmales=nmales,
        nfemales=nfemales,
        gene_list=gene_list,
        impute_missing_scores=impute_missing_scores,
        runtype=runtype,
        inframe_missense_ratio=inframe_missense_ratio,
        inframe_score_quantile=inframe_score_quantile,
        frameshift_nonsense_ratio=frameshift_nonsense_ratio,
        frameshift_score_quantile=frameshift_score_quantile,
        nsim=nsim,
        pvalcap=pvalcap,
        jobs=jobs,
        debug=debug,
        outdir=outdir,
        outfile=outfile,
    )

    # Set seed for reproducibility
    if debug:
        np.random.seed(42)

    # Initialize logger and log configuration
    init_log()
    log_configuration(asdict(cfg))

    # Load DNM and rates files
    dnm_df, rates_df = load_dnm_rates(dnm, rates, score_column, gene_list)

    # Prepare DNM and rates file for simulation
    dnm_df, rates_df = prepare_dnm_rates(dnm_df, rates_df, score_column, nmales, nfemales, cfg)

    # Run enrichment simulations
    results, logs = run_simulations(dnm_df, rates_df, score_column, cfg)

    # Export results
    export_results(results, outdir, outfile)

    # Export logs
    if debug:
        export_logs(logs, outdir)


if __name__ == "__main__":
    main()
