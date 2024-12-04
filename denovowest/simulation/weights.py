import pandas as pd
import numpy as np
import utils


def assign_weights(
    dnm_df,
    rates_df,
    score_column,
):
    """
    Assign weights to variants in the DNM and rates file.
    The score used for the weights is first min-max transformed (accounts for CEP with negative scores)

    Args:
        dnm_df (pd.DataFrame): DNM dataframe
        rates_df (pd.DataFrame): rates dataframe that contains all possible SNV
        score_column (str) : the column to use as weights
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

    # Apply min-max transform
    rates_df = min_max_transformation(rates_df, score_column, min_score, max_score)
    dnm_df = min_max_transformation(dnm_df, score_column, min_score, max_score)

    # Infer indel weights and rates
    indel_rates_df = infer_indel_weights_and_rates(rates_df)
    dnm_df = assign_dnm_indel_weights(dnm_df, indel_rates_df)

    # Impute weights for variants with missing scores
    dnm_df, rates_df = impute_missing_weights(dnm_df, rates_df)

    # Consolidate the rates df by adding the indel rates
    rates_df = pd.concat([rates_df, indel_rates_df])

    return dnm_df, rates_df


def min_max_transformation(df, column, min_score, max_score):
    """
    Transform scores between 0 and 1

    Args:
        df (pd.DataFrame): DNM or rates dataframe
        column (str): Column containing the CEP score to use
        min_score (float): Minimum score across rates and DNM for this CEP
        max_score (float): Maximum score across rates and DNM for this CEP

    Returns:
        pd.DataFrame: DNM or rates with normalized scores
    """

    df["ppv"] = df[column].apply(lambda x: (x - min_score) / (max_score - min_score))

    return df


def infer_indel_weights_and_rates(rates_df):
    """
    Infer expected inframe and frameshift weights and mutation rates per gene based
    on the missense and nonsense mutations respectively.

    Args:
        rates_df (pd.DataFrame): rates datafrane

    Returns:
        pd.DataFrame: per-gene inframe and frameshift mutation rates and associated weight
    """

    indel_rates_list = list()
    for gene_id, gene_df in rates_df.groupby("gene_id"):

        # Get the inframe mutation rate based on the cumulative missense mutation rates, and assign a weight from gene-based median missense weight
        gene_inframe_rate = gene_df.loc[gene_df.consequence == "missense", "prob"].sum() * utils.INFRAME_MISSENSE_RATIO
        gene_inframe_weights = gene_df.loc[gene_df.consequence == "missense", "ppv"].median()

        # Get the frameshift mutation rate based on the cumulative nonsense mutation rates, and assign a weight from gene-based median nonsense weight
        gene_frameshift_rate = (
            gene_df.loc[gene_df.consequence == "nonsense", "prob"].sum() * utils.FRAMESHIFT_NONSENSE_RATIO
        )
        gene_frameshift_weights = gene_df.loc[gene_df.consequence == "nonsense", "ppv"].median()

        inframe_row = {
            "gene_id": gene_id,
            "consequence": "inframe",
            "prob": gene_inframe_rate,
            "ppv": gene_inframe_weights,
        }

        frameshift_row = {
            "gene_id": gene_id,
            "consequence": "frameshift",
            "prob": gene_frameshift_rate,
            "ppv": gene_frameshift_weights,
        }

        indel_rates_list.append(inframe_row)
        indel_rates_list.append(frameshift_row)

    indel_rates_df = pd.DataFrame(indel_rates_list)

    return indel_rates_df


def assign_dnm_indel_weights(dnm_df, indel_rates_df):
    """
    Assign gene-based indel weights taken from the rates file to observed inframe
    and franeshift variants

    Args:
        dnm_df (pd.DataFrame): observed DNM
        rates_df (pd.DataFrame): expected mutations

    Returns:
        pd.DataFrame: observed DNM with weights associated to indels
    """

    list_weights = list()
    for _, dnm in dnm_df.iterrows():

        # If the DNM is not an indel we keep the current weight
        if not (dnm.consequence in ["inframe", "frameshift"]):
            list_weights.append(dnm.ppv)
        # Otherwise we take the corresponding gene-based indel weight from the indel rates dataframe
        else:
            list_weights.append(
                indel_rates_df.loc[
                    (indel_rates_df.gene_id == dnm.gene_id) & (indel_rates_df.consequence == dnm.consequence), "ppv"
                ].iloc[0]
            )

    dnm_df["ppv"] = list_weights
    return dnm_df


def impute_missing_weights(dnm_df, rates_df):
    """
    Some CEPs do not assign a score to each variant.
    Here we impute the missing weights looking at the median weight for each type
    of consequence per gene.

    Args:
        dnm_df (pd.DataFrame): observed DNM
        rates_df (pd.DataFrame): expected mutations

    Returns:
        tuple(pd.DataFrame, pd.DataFrame): DNM and rates dataframes with imputed missing weights
    """

    # Get the median values per gene and per consequence type
    median_values_dict = dict()
    for gene_id, generates_df in rates_df.groupby("gene_id"):
        median_values_dict[gene_id] = dict()
        for consequence, generates_cq_df in generates_df.groupby("consequence"):
            median_value = generates_cq_df["ppv"].median()
            median_values_dict[gene_id][consequence] = median_value

    # Impute missing weights on rates dataframe
    imputed_rates_weights = list()
    for _, variant in rates_df.iterrows():
        if np.isnan(variant.ppv):
            imputed_rates_weights.append(median_values_dict[variant.gene_id][variant.consequence])
        else:
            imputed_rates_weights.append(variant.ppv)
    rates_df["ppv_before_imputation"] = rates_df["ppv"]
    rates_df["ppv"] = imputed_rates_weights

    # Impute missing weights on DNM dataframe
    imputed_dnm_weights = list()
    for _, variant in dnm_df.iterrows():
        if np.isnan(variant.ppv):
            imputed_dnm_weights.append(median_values_dict[variant.gene_id][variant.consequence])
        else:
            imputed_dnm_weights.append(variant.ppv)

    dnm_df["ppv_before_imputation"] = dnm_df["ppv"]
    dnm_df["ppv"] = imputed_dnm_weights

    return dnm_df, rates_df
