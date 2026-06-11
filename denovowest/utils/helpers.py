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


def ensembl_ids_map_version(source_ids: list[str], target_ids: list[str]) -> dict[str, str]:
    """Map Ensembl IDs from a source list to the best-matching IDs in a target list.

    Handles all version-mismatch combinations:

    - source unversioned, target versioned  (``ENSG…`` → ``ENSG….N``)
    - source versioned, target unversioned  (``ENSG….N`` → ``ENSG…``)
    - source and target share the same version (identity)
    - source and target carry different versions of the same gene

    Exact matches take priority; base-ID matching (strip the ``.N`` suffix) is used
    as the fallback. For the base-ID lookup, the first encountered target ID wins
    if multiple target entries share the same base (version conflicts are rare).

    Args:
        source_ids: Ensembl gene IDs to map from (versioned or bare).
        target_ids: Reference Ensembl gene IDs to map to (versioned or bare).

    Returns:
        dict mapping each source ID to its best-matching target ID.  Source IDs
        with no match in ``target_ids`` are absent from the result.
    """
    exact_target: set[str] = set(target_ids)

    # Build a base-ID lookup: strip version suffix from every target ID
    base_to_target: dict[str, str] = {}
    for tid in target_ids:
        base = tid.split(".")[0]
        if base not in base_to_target:
            base_to_target[base] = tid

    mapping: dict[str, str] = {}
    for src in source_ids:
        if src in exact_target:
            mapping[src] = src
        else:
            base = src.split(".")[0]
            if base in base_to_target:
                mapping[src] = base_to_target[base]

    return mapping
