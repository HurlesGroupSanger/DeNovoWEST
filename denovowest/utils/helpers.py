import numpy as np
import pandas as pd


def compute_expected_number_of_mutations_stats(
    gene_rates: pd.Series, gene_infos: pd.DataFrame, nmales: int, nfemales: int
):

    # Compute the expected number of mutations as the product of the mutation rate and the number of individuals
    autosomal_factor = 2 * (nmales + nfemales)
    gene_expected = gene_rates * autosomal_factor

    if gene_infos.loc[gene_rates.name, "chrom"] == "chrX":
        x_factor_correction = compute_x_factor_correction(nmales, nfemales)
        gene_expected = x_factor_correction * gene_expected

    return gene_expected


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
