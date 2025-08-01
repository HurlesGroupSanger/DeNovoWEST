import pandas as pd
import numpy as np
import logging


from denovowest.utils.params import INFRAME_MISSENSE_RATIO, FRAMESHIFT_NONSENSE_RATIO
from denovowest.utils.log import set_plain_log, set_regular_log


def prepare_scores(dnm_df, rates_df, score_column, runtype, impute_missing=False):
    """
    Assign scores to variants in the DNM and rates file.
    The score used for the scores is first min-max transformed (accounts for CEP with negative scores)

    Args:
        dnm_df (pd.DataFrame): DNM dataframe
        rates_df (pd.DataFrame): rates dataframe that contains all possible SNV
        score_column (str) : CEP scores
    """

    # Infer indel scores and mutation rates
    if runtype == "ns":
        indel_rates_df = infer_indel_scores_and_rates(rates_df, score_column)
        dnm_df = assign_dnm_indel_scores(dnm_df, indel_rates_df, rates_df, score_column)

    # Impute scores for variants with missing scores
    if impute_missing:
        dnm_df, rates_df = impute_missing_scores(dnm_df, rates_df, score_column)
    else:
        dnm_df, rates_df = remove_missing_scores(dnm_df, rates_df, score_column)

    # Consolidate the rates df by adding the indel rates
    if runtype == "ns":
        rates_df = pd.concat([rates_df, indel_rates_df])

    return dnm_df, rates_df


def infer_indel_scores_and_rates(rates_df, score_column):
    """
    Infer expected inframe and frameshift scores and mutation rates per gene based
    on the missense and nonsense mutations respectively.

    Args:
        rates_df (pd.DataFrame): rates dataframe
        score_column (str) : CEP scores


    Returns:
        pd.DataFrame: per-gene inframe and frameshift mutation rates and associated scores
    """

    indel_rates_list = list()
    for gene_id, gene_df in rates_df.groupby("gene_id"):

        # Get the inframe mutation rate based on the cumulative missense mutation rates, and assign a score from gene-based median missense score
        gene_inframe_rate = gene_df.loc[gene_df.consequence == "missense", "prob"].sum() * INFRAME_MISSENSE_RATIO
        gene_inframe_scores = gene_df.loc[gene_df.consequence == "missense", score_column].median()

        # Get the frameshift mutation rate based on the cumulative nonsense mutation rates, and assign a score from gene-based median nonsense score
        gene_frameshift_rate = gene_df.loc[gene_df.consequence == "nonsense", "prob"].sum() * FRAMESHIFT_NONSENSE_RATIO
        gene_frameshift_scores = gene_df.loc[gene_df.consequence == "nonsense", score_column].median()

        inframe_row = {
            "gene_id": gene_id,
            "consequence": "inframe",
            "prob": gene_inframe_rate,
            score_column: gene_inframe_scores,
        }

        frameshift_row = {
            "gene_id": gene_id,
            "consequence": "frameshift",
            "prob": gene_frameshift_rate,
            score_column: gene_frameshift_scores,
        }

        indel_rates_list.append(inframe_row)
        indel_rates_list.append(frameshift_row)

    indel_rates_df = pd.DataFrame(indel_rates_list)

    return indel_rates_df


def assign_dnm_indel_scores(dnm_df, indel_rates_df, rates_df, score_column):
    """
    Assign gene-based indel scores taken from the rates file to observed inframe
    and frameshift variants

    Args:
        dnm_df (pd.DataFrame): observed DNM
        indel_rates_df (pd.DataFrame) : per-gene score for frameshit and inframe annotated indels
        rates_df (pd.DataFrame): rates dataframe
        score_column (str) : CEP scores

    Returns:
        pd.DataFrame: observed DNM with scores associated to indels
    """

    list_scores = list()
    for _, dnm in dnm_df.iterrows():

        # If the DNM is not an indel we keep the current score
        if (len(dnm.alt) - len(dnm.ref)) == 0:
            list_scores.append(dnm[score_column])
            continue

        # If it is an indel but not annotated as inframe or frameshift, we get the
        # gene median corresponding score
        if not (dnm.consequence in ["inframe", "frameshift"]):

            generates_df = rates_df.loc[rates_df.gene_id == dnm.gene_id]
            score = generates_df.loc[generates_df.consequence == dnm.consequence, score_column].median()
            list_scores.append(score)

        # Otherwise we get the gene inframe/frameshift score
        else:
            list_scores.append(
                indel_rates_df.loc[
                    (indel_rates_df.gene_id == dnm.gene_id) & (indel_rates_df.consequence == dnm.consequence),
                    score_column,
                ].iloc[0]
            )

    dnm_df.loc[:, score_column] = list_scores
    return dnm_df


def impute_missing_scores(dnm_df, rates_df, score_column):
    """
    Some CEPs do not assign a score to each variant.
    Here we impute the missing scores looking at the median score for each type
    of consequence per gene.

    Args:
        dnm_df (pd.DataFrame): observed DNM
        rates_df (pd.DataFrame): expected mutations
        score_column (str) : score to use for the simulation


    Returns:
        tuple(pd.DataFrame, pd.DataFrame): DNM and rates dataframes with imputed missing scores
    """

    logger = logging.getLogger("logger")

    set_plain_log()
    logger.info("=" * 50)
    logger.info("[MISSING SCORES]")
    logger.info("=" * 50)
    set_regular_log()

    logger.info("Imputing missing scores")

    # Get the median values per gene and per consequence type
    median_values_dict = dict()
    for gene_id, generates_df in rates_df.groupby("gene_id"):
        median_values_dict[gene_id] = dict()
        for consequence, generates_cq_df in generates_df.groupby("consequence"):
            median_value = generates_cq_df[score_column].median()
            median_values_dict[gene_id][consequence] = median_value

    # Impute missing scores on rates dataframe
    imputed_rates_scores = list()
    for _, variant in rates_df.iterrows():
        if np.isnan(variant[score_column]):
            imputed_rates_scores.append(median_values_dict[variant.gene_id][variant.consequence])
        else:
            imputed_rates_scores.append(variant[score_column])
    rates_df["score_before_imputation"] = rates_df[score_column]
    rates_df[score_column] = imputed_rates_scores

    # Impute missing scores on DNM dataframe
    imputed_dnm_scores = list()
    for _, variant in dnm_df.iterrows():
        if np.isnan(variant[score_column]):
            imputed_dnm_scores.append(median_values_dict[variant.gene_id][variant.consequence])
        else:
            imputed_dnm_scores.append(variant[score_column])

    dnm_df["score_before_imputation"] = dnm_df[score_column]
    dnm_df[score_column] = imputed_dnm_scores

    return dnm_df, rates_df


def remove_missing_scores(dnm_df, rates_df, score_column):
    """
    Some CEPs do not assign a score to each variant. Or depending on the source (dbNSFP) some
    annotations can be missing.
    Here we remove all records that do not have a score assigned.

    Args:
        dnm_df (pd.DataFrame): observed DNM
        rates_df (pd.DataFrame): expected mutations
        score_column (str) : score to use for the simulation

    Returns:
        tuple(pd.DataFrame, pd.DataFrame): DNM and rates dataframes with records with missing scores removed
    """

    logger = logging.getLogger("logger")

    set_plain_log()
    logger.info("=" * 50)
    logger.info("[MISSING SCORES]")
    logger.info("=" * 50)
    set_regular_log()

    nb_records_dnm_before = dnm_df.shape[0]
    nb_records_rates_before = rates_df.shape[0]

    dnm_df = dnm_df.loc[~dnm_df[score_column].isna()]
    rates_df = rates_df.loc[~rates_df[score_column].isna()]

    nb_records_dnm_after = dnm_df.shape[0]
    nb_records_rates_after = rates_df.shape[0]

    if nb_records_dnm_after - nb_records_dnm_before:
        logger.info(
            f"{nb_records_dnm_before - nb_records_dnm_after } observed DNMs were removed as they do not have a {score_column} score"
        )

    if nb_records_rates_after - nb_records_rates_before:
        logger.info(
            f"{nb_records_rates_before - nb_records_rates_after} variants were removed from the rates file as they do not have a {score_column} score"
        )

    return dnm_df, rates_df
