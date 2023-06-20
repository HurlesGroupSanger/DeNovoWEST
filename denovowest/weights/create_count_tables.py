#!/usr/bin/env python

import pandas as pd
import numpy as np
import itertools
import click
import utils


LOWER_BOUND_CADD_MISSENSE = 6
UPPER_BOUND_CADD_MISSENSE = 30
BINSIZE_CADD_MISSENSE = 6
LOWER_BOUND_CADD_NONSENSE_SHETLOW = 7.5
UPPER_BOUND_CADD_NONSENSE_SHETLOW = 45
BINSIZE_CADD_NONSENSE_SHETLOW = 7.5
LOWER_BOUND_CADD_NONSENSE_SHETHIGH = 15
UPPER_BOUND_CADD_NONSENSE_SHETHIGH = 45
BINSIZE_CADD_NONSENSE_SHETHIGH = 15

MIN_SCORE = 0
MAX_SCORE = 50

SHET_VALUES = [True, False]
CONSTRAINED_VALUES = [True, False]
CSQ_SYNONYMOUS_VALUES = ["synonymous", "splice_lof", "missense"]
CSQ_MISSENSE_VALUES = ["missense"]
CSQ_NONSENSE_VALUES = ["nonsense"]

NA = ["NA"]


@click.group()
def cli():
    pass


@cli.command()
@click.argument("rates_file")
@click.argument("nmales")
@click.argument("nfemales")
@click.argument("outfile")
def get_expected_counts(rates_file, nmales, nfemales, outfile):
    """
    Calculates the expected number of mutations in the cohort per category, given an annotated mutation rates file and
    the number of males and females individuals in the cohort. Mutations are classified in bins according
    to several variables : functional consequence, CADD score and whether they fall in a constrained region
    or a high/low shet gene.

    Args:
        rates_file (str): Path to the annotated mutation rates file
        nmales (int): Number of male individuals in the cohort
        nfemales (int): Number of female individuals in the cohort
        outfile (str): Path to the output table listing the expected number of mutations per category
    """

    utils.init_log()

    nmales = int(nmales)
    nfemales = int(nfemales)

    # Load the mutation rates file
    rates_df = pd.read_csv(
        rates_file,
        sep="\t",
        dtype={
            "chrom": str,
            "pos": int,
            "ref": str,
            "alt": str,
            "consequence": str,
            "score": float,
            "prob": float,
            "shethigh": bool,
            "constrained": bool,
        },
        usecols=["chrom", "pos", "ref", "alt", "consequence", "score", "prob", "shethigh", "constrained"],
    )

    # Only a subset of functional consequences are used
    rates_df = utils.filter_on_consequences(rates_df)
    rates_df = utils.assign_meta_consequences(rates_df)

    # Filter out missense|nonsense loci with missing CADD score
    rates_df = rates_df.loc[~(rates_df.consequence.isin(["missense", "nonsense"]) & rates_df.score.isna())]

    # Adjust mutation rates by chromosomal scaling factor
    autosomal_factor = 2 * (nmales + nfemales)
    x_factor = get_X_scaling_factor(nmales, nfemales)
    rates_df = rates_df.assign(
        exp=np.where(rates_df.chrom == "X", rates_df.prob * x_factor, rates_df.prob * autosomal_factor)
    )

    # Get expected number of mutations
    bins = define_bins()
    expected_df = create_expected_counts_table(rates_df, bins)
    expected_df.to_csv(outfile, sep="\t", index=False)


@cli.command()
@click.argument("dnm_file")
@click.argument("outfile")
def get_observed_counts(dnm_file, outfile):
    """
    Calculates the observed number of mutations in a DNM file

    Args:
        dnm_file (str): Path to a DNM file
        outfile (str): Output file
    """

    utils.init_log()

    # Load the DNM file
    dnm_df = pd.read_csv(
        dnm_file,
        sep="\t",
        dtype={
            "chrom": str,
            "pos": int,
            "ref": str,
            "alt": str,
            "consequence": str,
            "score": float,
            "shethigh": bool,
            "constrained": bool,
        },
        usecols=["chrom", "pos", "ref", "alt", "consequence", "score", "shethigh", "constrained"],
    )

    # Only a subset of functional consequences are used
    dnm_df = utils.filter_on_consequences(dnm_df)
    dnm_df = utils.assign_meta_consequences(dnm_df)

    # Filter out missense|nonsense with missing scores
    dnm_df = dnm_df.loc[~(dnm_df.consequence.isin(["missense", "nonsense"]) & dnm_df.score.isna())]

    # Get observed number of mutations
    bins = define_bins()
    observed_df = pd.concat([count_variants(dnm_df, condition, "obs") for condition in bins], axis=0, ignore_index=True)

    observed_df.to_csv(outfile, sep="\t", index=False)


def get_X_scaling_factor(male_n: int, female_n: int):
    """
    Scaling factors take into account the number of males and females individuals in the cohort,
    correct for the different inheritance pattern and the different mutation rates in the male and female germline for the X chromosome.
    They are further used to calculate the true expected number of mutations given those parameters.

    More informations about the scaling factors can be found in https://www.nature.com/articles/s41467-020-20852-3

    Args:
        male_n (int):  Number of males individuals in the cohort
        female_n (int):  Number of females individuals in the cohort

    Returns:
        float: scaling factor for X chromosome
    """

    ALPHA = 3.4  # ratio of the mutation rate in fathers to mothers in DDD
    male_factor = 2 / (1 + (1 / ALPHA))
    female_factor = 2 / (1 + ALPHA)

    male_transmissions = female_n
    female_transmissions = male_n + female_n
    x_factor = (male_transmissions * male_factor) + (female_transmissions * female_factor)

    return x_factor


def define_bins() -> list:
    """
    Builds the bins used to categorize mutations.

    Returns:
        list: list of bins
    """

    synonymous_splice_inframe_bins = define_bins_synonymous()
    missense_bins = define_bins_missense()
    nonsense_bins = define_bins_nonsense()

    return itertools.chain(synonymous_splice_inframe_bins, missense_bins, nonsense_bins)


def define_bins_synonymous():
    """
    Builds the bins containing variants having a synonymous functional consequence.
    CADD score is not considered here.

    Returns:
        list: list of bins
    """

    synonymous_splice_inframe = list(itertools.product(CSQ_SYNONYMOUS_VALUES, NA, NA, SHET_VALUES))
    return synonymous_splice_inframe


def define_bins_missense():
    """
    Builds the bins containing variants having a missense functional consequence

    Returns:
        list: list of bins
    """

    cadd_missense = generate_cadd_bins(LOWER_BOUND_CADD_MISSENSE, UPPER_BOUND_CADD_MISSENSE, BINSIZE_CADD_MISSENSE)

    missense_bins = list(itertools.product(CSQ_MISSENSE_VALUES, cadd_missense, CONSTRAINED_VALUES, SHET_VALUES))

    return missense_bins


def define_bins_nonsense():
    """
    Builds the bins containing variants having a nonsense functional consequence
    """

    cadd_nonsense_low = generate_cadd_bins(
        LOWER_BOUND_CADD_NONSENSE_SHETLOW, UPPER_BOUND_CADD_NONSENSE_SHETLOW, BINSIZE_CADD_NONSENSE_SHETLOW
    )
    cadd_nonsense_high = generate_cadd_bins(
        LOWER_BOUND_CADD_NONSENSE_SHETHIGH, UPPER_BOUND_CADD_NONSENSE_SHETHIGH, BINSIZE_CADD_NONSENSE_SHETHIGH
    )

    nonsense_bins_shethigh = list(itertools.product(CSQ_NONSENSE_VALUES, cadd_nonsense_high, NA, [SHET_VALUES[0]]))
    nonsense_bins_shetlow = list(itertools.product(CSQ_NONSENSE_VALUES, cadd_nonsense_low, NA, [SHET_VALUES[1]]))

    return itertools.chain(nonsense_bins_shethigh, nonsense_bins_shetlow)


def generate_cadd_bins(lower_bound: float, upper_bound: float, bin_size: float) -> list:
    """
    Create bins defined by CADD scores.
    One bin is created to capture all loci with a CADD score < lower_bound
    One bin is created to capture all loci with a CADD score > upper_bound
    Several intermediary bins of length bin_size are created between lower_bound and upper_bound CADD scores.

    Args:
        lower_bound (float): CADD score defining the upper bound of the lower bin
        upper_bound (float): CADD score defining the lower bound of the upper bin
        bin_size (float): Size of the intermediary bins

    Returns:
        list: list of CADD scores bins
    """

    lower_bin = [[MIN_SCORE, lower_bound]]
    upper_bin = [[upper_bound, MAX_SCORE]]

    intermediate_bin = [[step, step + bin_size] for step in np.arange(lower_bound, upper_bound, bin_size)]

    return itertools.chain(lower_bin, intermediate_bin, upper_bin)


def create_expected_counts_table(rates_df: pd.DataFrame, bins: list):
    """
    Calculate and assign expected number of mutations per bin.

    Args:
        rates_df (pd.DataFrame): mutatation rates data frame
        bins (list): list of bins
    """

    expected_df = pd.concat([count_variants(rates_df, bin, "exp") for bin in bins], axis=0, ignore_index=True)
    return expected_df


def count_variants(df: pd.DataFrame, bin: list, mode: str) -> pd.DataFrame:
    """
    Calculate the expected or count the observed number of mutations in mutation rates or DNM file

    Args:
        df (pd.DataFrame): Mutation rates or DNM data frame
        bin (list): List of bins
        mode (str): exp for expected, obs for observed

    Returns:
        count_df (pd.DataFrame): Data frame with each row corresponding to a bin with a value corresponding to
        the observed or expected number of mutations falling in it
    """

    # TODO : replace bin list per dict to make it more explicit
    filter_consequence = df.consequence == bin[0]

    # Add score filter for bins using the score information (e.g. missense)
    try:
        filter_score = (df.score >= bin[1][0]) & (df.score < bin[1][1])
        score_range = f"{bin[1][0]}-{bin[1][1]}"
    except TypeError as e:
        filter_score = True
        score_range = "NA"

    # Some bins (e.g. synonymous variants) do not use the constraint regions information
    if bin[2] == "NA":
        filter_constrained = True
    else:
        filter_constrained = df.constrained == bin[2]

    # All bins use the shet information
    filter_shet = df.shethigh == bin[3]

    # Get only variants falling in the current bin
    min_df = df.loc[filter_consequence & filter_score & filter_constrained & filter_shet]

    # Count number of expected or observed variants in the bin
    if mode == "exp":
        res = min_df["exp"].sum()
    else:
        res = len(min_df)

    # Store this result in a data frame
    count_df = pd.DataFrame(
        {"consequence": bin[0], "score": score_range, "constrained": bin[2], "shethigh": bin[3], mode: res}, index=[0]
    )

    return count_df


@cli.command()
@click.argument("expected_file")
@click.argument("observed_file")
@click.argument("outfile")
def merge_expected_observed(expected_file: str, observed_file: str, outfile: str):
    """
    Merge observed and expected number of mutations per bin, and compute their ratio.

    Args:
        expected_file (str): Path to the dataframe listing expected number of mutations per bin
        observed_file (str): Path to the dataframe listing observed number of mutations per bin
    """

    df_expected = pd.read_csv(expected_file, sep="\t")
    df_observed = pd.read_csv(observed_file, sep="\t")

    df_merged = df_expected.merge(df_observed, on=["consequence", "score", "constrained", "shethigh"], how="outer")

    df_merged["obs_exp"] = df_merged["obs"] / df_merged["exp"]

    df_merged.to_csv(outfile, sep="\t", index=False)


if __name__ == "__main__":
    cli()
