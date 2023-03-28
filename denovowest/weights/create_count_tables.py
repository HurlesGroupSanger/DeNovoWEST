import pandas as pd
import numpy as np
import itertools
import os
import click

SHET_BINS = ["shethigh == True", "shethigh == False"]
MCR_BINS = ["constrained == True", "constrained == False"]

LOWER_BOUND_CADD_MISSENSE = 6
UPPER_BOUND_CADD_MISSENSE = 30
BINSIZE_CADD_MISSENSE = 6
LOWER_BOUND_CADD_NONSENSE_SHETLOW = 7.5
UPPER_BOUND_CADD_NONSENSE_SHETLOW = 45
BINSIZE_CADD_NONSENSE_SHETLOW = 7.5
LOWER_BOUND_CADD_NONSENSE_SHETHIGH = 15
UPPER_BOUND_CADD_NONSENSE_SHETHIGH = 45
BINSIZE_CADD_NONSENSE_SHETHIGH = 15


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
    Calculates the expected number of mutations per gene given a categorized mutability rate file and
    the number of males and females individual in the cohort.

    Args:
        rates_file (str): Path to a mutation rates file
        nmales (int): Number of males individuals in the cohort
        nfemales (int): Number of females individuals in the cohort
        outfile (str): Output file
    """

    nmales = int(nmales)
    nfemales = int(nfemales)

    # Load the mutation rates file
    rates_df = pd.read_csv(
        rates_file,
        sep="\t",
        dtype={"score": float, "prob": float, "chrom": str, "consequence": str},
        usecols=["chrom", "ref", "alt", "consequence", "score", "prob", "shethigh", "constrained"],
    )

    # Filter out missense|nonsense with missing scores
    rates_df = rates_df.loc[
        ~(rates_df.consequence.isin(["missense", "nonsense", "stop_gained"]) & rates_df.score.isna())
    ]

    # Adjust mutation rates by chromosomal scaling factor
    autosomal_factor = 2 * (nmales + nfemales)
    x_factor = get_X_scaling_factor(nmales, nfemales)
    rates_df = rates_df.assign(
        exp=np.where(rates_df.chrom == "X", rates_df.prob * x_factor, rates_df.prob * autosomal_factor)
    )

    # Get expected number of mutations
    expected_df = pd.concat(
        [count_variants(rates_df, condition, "exp") for condition in define_bins("exp")], axis=0, ignore_index=True
    )

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

    # Load the DNM file
    dnm_df = pd.read_csv(
        dnm_file,
        sep="\t",
        dtype={"score": float, "chrom": str, "consequence": str},
        usecols=["chrom", "ref", "alt", "consequence", "score", "shethigh", "constrained"],
    )

    # Filter out missense|nonsense with missing scores
    dnm_df = dnm_df.loc[~(dnm_df.consequence.isin(["missense", "nonsense", "stop_gained"]) & dnm_df.score.isna())]

    # Get observed number of mutations
    observed_df = pd.concat(
        [count_variants(dnm_df, condition, "obs") for condition in define_bins("obs")], axis=0, ignore_index=True
    )

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


def define_bins(mode: str) -> list:
    """
    Builds a query combining different attributes (CADD score, shet, consequence etc...) to create bins into each
    variant/loci will be categorized.

    Args:
        mode (str): exp for expected, obs for observed

    Returns:
        list: list of bins (defined by their query)
    """

    synonymous_splice_inframe_bins = define_bins_synonymous(mode)
    missense_bins = define_bins_missense()
    nonsense_bins = define_bins_nonsense()

    return [" and ".join(list(x)) for x in synonymous_splice_inframe_bins + missense_bins + nonsense_bins]


def define_bins_synonymous(mode):
    """
    Builds a query combining different attributes (shet, consequence) to create meta bins for synonymous variants
    """

    if mode == "obs":

        synonymous_splice_inframe = list(
            itertools.product(
                SHET_BINS,
                [
                    f'consequence.str.contains("{cq}")'
                    for cq in ["synonymous", "splice_donor|splice_acceptor|splice_lof", "missense"]
                ],
                ["(ref.str.len() - alt.str.len()) == 0"],
            )
        )
    else:
        synonymous_splice_inframe = list(
            itertools.product(
                SHET_BINS,
                [
                    f'consequence.str.contains("{cq}")'
                    for cq in ["synonymous", "splice_donor|splice_acceptor|splice_lof", "missense"]
                ],
            )
        )

    return synonymous_splice_inframe


def define_bins_missense():
    """
    Builds a query combining different attributes (CADD, shet, consequence etc..) to create bins for missense variants
    """

    cadd_missense = generate_cadd_bins(LOWER_BOUND_CADD_MISSENSE, UPPER_BOUND_CADD_MISSENSE, BINSIZE_CADD_MISSENSE)

    missense_bins = list(
        itertools.product(cadd_missense, SHET_BINS, MCR_BINS, ['consequence.str.contains("missense")'])
    )

    return missense_bins


def define_bins_nonsense():
    """
    Builds a query combining different attributes (CADD, shet, consequence etc..) to create bins for nonsense variants
    """

    cadd_nonsense_low = generate_cadd_bins(
        LOWER_BOUND_CADD_NONSENSE_SHETLOW, UPPER_BOUND_CADD_NONSENSE_SHETLOW, BINSIZE_CADD_NONSENSE_SHETLOW
    )
    cadd_nonsense_high = generate_cadd_bins(
        LOWER_BOUND_CADD_NONSENSE_SHETHIGH, UPPER_BOUND_CADD_NONSENSE_SHETHIGH, BINSIZE_CADD_NONSENSE_SHETHIGH
    )

    nonsense_bins = list(
        itertools.product(
            cadd_nonsense_high, ["shethigh == True"], ['consequence.str.contains("nonsense|stop_gained")']
        )
    ) + list(
        itertools.product(
            cadd_nonsense_low, ["shethigh == False"], ['consequence.str.contains("nonsense|stop_gained")']
        )
    )

    return nonsense_bins


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

    lower_bin = [f'score.between(0, {lower_bound}, inclusive="left")']
    upper_bin = [f"score >= {upper_bound}"]

    intermediate_bin = [
        f'score.between({step}, {min(step + lower_bound, upper_bound)}, inclusive="left")'
        for step in np.arange(lower_bound, upper_bound, bin_size)
    ]

    return itertools.chain(lower_bin, intermediate_bin, upper_bin)


def count_variants(df: pd.DataFrame, condition: str, mode: str) -> pd.DataFrame:
    """
    Count the expected or observed number of mutations in mutation rates or DNM file

    Args:
        df (pd.DataFrame): Observed DNM or rates file
        condition (str): Formal description of the bin
        mode (str): exp for expected, obs for observed

    Returns:
        count_df (pd.DataFrame): Data frame with each row corresponding to a bin with a value corresponding to
        the observed or expected number of mutations falling in it
    """

    if mode == "exp":
        count_df = pd.DataFrame({"bin": condition, "exp": df.query(condition)["exp"].sum()}, index=[0])
    else:
        count_df = pd.DataFrame({"bin": condition, "obs": len(df.query(condition))}, index=[0])

    return count_df


if __name__ == "__main__":
    cli()
