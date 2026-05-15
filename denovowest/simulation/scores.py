import logging

import numpy as np
import pandas as pd

from denovowest.utils.log import set_plain_log, set_regular_log
from denovowest.utils.params import RATIO_SPLICE_REGION_INDEL_SNV, RATIO_SPLICE_SITE_INDEL_SNV, RunType


def prepare_scores(dnm_df, rates_df, score_column, prep_logs, cfg):
    """
    Assign scores to variants in the DNM and rates file.
    The score used for the scores is first min-max transformed (accounts for CEP with negative scores)

    Args:
        dnm_df (pd.DataFrame): DNM dataframe
        rates_df (pd.DataFrame): rates dataframe that contains all possible SNV
        score_column (str) : CEP scores
        prep_logs (dict): preparation logs structure
        cfg (Config): configuration object that stores script parameters
    """

    # Infer indel scores and mutation rates
    if cfg.runtype in [RunType.ALL_CODING, RunType.PROTEIN_ALTERING, RunType.PTV]:
        indel_rates_df = infer_indel_scores_and_rates(rates_df, score_column, cfg)
        dnm_df = assign_dnm_indel_scores(dnm_df, indel_rates_df, rates_df, score_column)

    # Impute scores for variants with missing scores
    if cfg.impute_missing_scores:
        prep_logs, dnm_df, rates_df = impute_missing_scores(dnm_df, rates_df, score_column, prep_logs)
    else:
        prep_logs, dnm_df, rates_df = remove_missing_scores(dnm_df, rates_df, score_column, prep_logs)

    # Consolidate the rates df by adding the indel rates
    if cfg.runtype in [RunType.ALL_CODING, RunType.PROTEIN_ALTERING, RunType.PTV]:
        rates_df = pd.concat([rates_df, indel_rates_df])

    return prep_logs, dnm_df, rates_df


def infer_indel_scores_and_rates(rates_df, score_column, cfg):
    """
    Infer expected inframe and frameshift scores and mutation rates per gene based
    on the missense and nonsense mutations respectively.

    Args:
        rates_df (pd.DataFrame): rates dataframe
        score_column (str) : CEP scores
        cfg (Config): configuration object that stores script parameters


    Returns:
        pd.DataFrame: per-gene inframe and frameshift mutation rates and associated scores
    """

    indel_rates_list = list()
    for gene_id, gene_df in rates_df.groupby("gene_id"):

        ##### INFRAME #####

        # Get the inframe mutation rate based on the cumulative missense mutation rates, and assign a score from gene-based median missense score
        gene_inframe_rate = gene_df.loc[gene_df.consequence == "missense", "prob"].sum() * cfg.inframe_missense_ratio
        gene_inframe_scores = gene_df.loc[gene_df.consequence == "missense", score_column].quantile(
            cfg.inframe_score_quantile
        )

        inframe_row = {
            "gene_id": gene_id,
            "consequence": "inframe",
            "prob": gene_inframe_rate,
            score_column: gene_inframe_scores,
        }

        indel_rates_list.append(inframe_row)

        ##### FRAMESHIFT #####

        # Get the frameshift mutation rate based on the cumulative nonsense mutation rates, and assign a score from gene-based median nonsense score
        gene_frameshift_rate = (
            gene_df.loc[gene_df.consequence == "nonsense", "prob"].sum() * cfg.frameshift_nonsense_ratio
        )
        gene_frameshift_scores = gene_df.loc[gene_df.consequence == "nonsense", score_column].quantile(
            cfg.frameshift_score_quantile
        )

        frameshift_row = {
            "gene_id": gene_id,
            "consequence": "frameshift",
            "prob": gene_frameshift_rate,
            score_column: gene_frameshift_scores,
        }

        indel_rates_list.append(frameshift_row)

        ##### SPLICE_LOF #####

        # Get the indel splice_lof mutation rate based on the cumulative SNV splice_lof rates, and assign a ascore from gene-based median splice_lof score
        gene_splicesite_indel_rate = (
            gene_df.loc[gene_df.consequence == "splice_lof", "prob"].sum() * RATIO_SPLICE_SITE_INDEL_SNV
        )
        gene_splicesite_indel_scores = gene_df.loc[gene_df.consequence == "splice_lof", score_column].median()

        # Some genes have a single exon, introducing NaN lead to calculation errors later
        if gene_splicesite_indel_rate != 0:

            splice_lof_row = {
                "gene_id": gene_id,
                "consequence": "splice_lof",
                "prob": gene_splicesite_indel_rate,
                score_column: gene_splicesite_indel_scores,
            }

            indel_rates_list.append(splice_lof_row)

        ##### SPLICE_REGION #####

        # Get the indel splice_region mutation rate based on the cumulative SNV splice_region rates, and assign a ascore from gene-based median splice_region score
        if cfg.runtype == RunType.ALL_CODING:
            gene_spliceregion_indel_rate = (
                gene_df.loc[gene_df.consequence == "splice_region", "prob"].sum() * RATIO_SPLICE_REGION_INDEL_SNV
            )
            gene_spliceregion_indel_scores = gene_df.loc[gene_df.consequence == "splice_region", score_column].median()

            # Some genes have a single exon, introducing NaN lead to calculation errors later
            if gene_spliceregion_indel_rate != 0:

                splice_region_row = {
                    "gene_id": gene_id,
                    "consequence": "splice_region",
                    "prob": gene_spliceregion_indel_rate,
                    score_column: gene_spliceregion_indel_scores,
                }

                indel_rates_list.append(splice_region_row)

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

        # If the DNM is an inframe or frameshift we get the corresponding score in the gene
        if dnm.consequence in ["inframe", "frameshift", "splice_lof", "splice_region"]:

            list_scores.append(
                indel_rates_df.loc[
                    (indel_rates_df.gene_id == dnm.gene_id) & (indel_rates_df.consequence == dnm.consequence),
                    score_column,
                ].iloc[0]
            )

        # A handful of indels in CDS will be annotated as a category for which we do not directly compute an indel mutation rate (e.g. stop_gained, start_gained)
        # In these rare cases we retrieve the median score for this category in the gene.
        else:

            mask = (rates_df.gene_id == dnm.gene_id) & (rates_df.consequence == dnm.consequence)
            score = rates_df.loc[mask, score_column].median()
            list_scores.append(score)

    dnm_df.loc[:, score_column] = list_scores
    return dnm_df


def impute_missing_scores(dnm_df, rates_df, score_column, prep_logs):
    """
    Some CEPs do not assign a score to each variant.
    Here we impute the missing scores looking at the median score for each type
    of consequence per gene.

    Args:
        dnm_df (pd.DataFrame): observed DNM
        rates_df (pd.DataFrame): expected mutations
        score_column (str) : score to use for the simulation
        prep_logs (dict): preparation logs structure


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

    # Get the median scores per gene and per consequence type
    median_values_dict = dict()
    for gene_id, generates_df in rates_df.groupby("gene_id"):
        median_values_dict[gene_id] = dict()
        for consequence, generates_cq_df in generates_df.groupby("consequence"):
            median_value = generates_cq_df[score_column].median()
            median_values_dict[gene_id][consequence] = median_value

    # Identify genes with all missing scores
    list_genes_all_scores_missing = list()
    for gene_id, gene_dict in median_values_dict.items():
        median_values = gene_dict.values()
        if all(np.isnan(value) for value in median_values):
            logger.warning(
                f"All scores are missing for gene {gene_id}. Reverting to a classic burden test for this gene."
            )
            list_genes_all_scores_missing.append(gene_id)
            dnm_df.loc[dnm_df.gene_id == gene_id, score_column] = np.nan
            prep_logs[gene_id]["all_scores_missing"] = True

    ##### RATES #####

    # Impute missing scores for rates dataframe based on median values per gene and consequence type
    imputed_rates_scores = list()
    for _, variant in rates_df.iterrows():
        if np.isnan(variant[score_column]):
            # If all scores are missing for the gene, assign a score of 1 to all variants which revert to a classic burden test
            if variant.gene_id in list_genes_all_scores_missing:
                imputed_rates_scores.append(1)
            else:
                imputed_rates_scores.append(median_values_dict[variant.gene_id][variant.consequence])
        else:
            imputed_rates_scores.append(variant[score_column])
    rates_df.loc[:, "score_before_imputation"] = rates_df[score_column]
    rates_df.loc[:, score_column] = imputed_rates_scores

    # Log imputing results on rates
    for consequence, consequence_rates_df in rates_df.groupby("consequence"):
        nb_imputed = sum(
            (~consequence_rates_df[score_column].isna()) & (consequence_rates_df["score_before_imputation"].isna())
        )
        nb_still_no_score = sum(consequence_rates_df[score_column].isna())
        logger.info(f"RATES - {consequence} : {nb_imputed} / {consequence_rates_df.shape[0]} imputed")

        if nb_still_no_score:
            logger.warning(f"RATES - {consequence} : {nb_still_no_score} were not imputed")

    for gene_id, generates_df in rates_df.groupby("gene_id"):
        nb_imputed = sum((~generates_df[score_column].isna()) & (generates_df["score_before_imputation"].isna()))
        prep_logs[gene_id]["nb_rates_variants_imputed"] = nb_imputed

        nb_still_no_score = sum(generates_df[score_column].isna())
        if nb_still_no_score:
            prep_logs[gene_id]["nb_rates_variants_still_no_score"] = nb_still_no_score

    ##### DNM #####

    # Impute missing scores for DNM dataframe based on median values per gene and consequence type
    imputed_dnm_scores = list()
    for _, variant in dnm_df.iterrows():
        if np.isnan(variant[score_column]):
            # If all scores are missing for the gene, assign a score of 1 to all variants which revert to a classic burden test
            if variant.gene_id in list_genes_all_scores_missing:
                imputed_dnm_scores.append(1)
            else:
                imputed_dnm_scores.append(median_values_dict[variant.gene_id][variant.consequence])
        else:
            imputed_dnm_scores.append(variant[score_column])

    dnm_df.loc[:, "score_before_imputation"] = dnm_df[score_column]
    dnm_df.loc[:, score_column] = imputed_dnm_scores

    # Log inmputing results on DNM
    for consequence, consequence_dnm_df in dnm_df.groupby("consequence"):
        nb_imputed = sum(
            (~consequence_dnm_df[score_column].isna()) & (consequence_dnm_df["score_before_imputation"].isna())
        )
        nb_still_no_score = sum(consequence_dnm_df[score_column].isna())
        logger.info(f"DNM - {consequence} : {nb_imputed} / {consequence_dnm_df.shape[0]} imputed")
        if nb_still_no_score:
            logger.warning(f"DNM - {consequence} : {nb_still_no_score} were not imputed")

    for gene_id, gene_dnm_df in dnm_df.groupby("gene_id"):
        nb_imputed = sum((~gene_dnm_df[score_column].isna()) & (gene_dnm_df["score_before_imputation"].isna()))
        prep_logs[gene_id]["nb_observed_variants_imputed"] = nb_imputed

        nb_still_no_score = sum(gene_dnm_df[score_column].isna())
        if nb_still_no_score:
            prep_logs[gene_id]["nb_observed_variants_still_no_score"] = nb_still_no_score

    return prep_logs, dnm_df, rates_df


def remove_missing_scores(dnm_df, rates_df, score_column, prep_logs):
    """
    Some CEPs do not assign a score to each variant. Or depending on the source (dbNSFP) some
    annotations can be missing.
    Here we remove all records that do not have a score assigned.

    Args:
        dnm_df (pd.DataFrame): observed DNM
        rates_df (pd.DataFrame): expected mutations
        score_column (str) : score to use for the simulation
        prep_logs (dict): preparation logs structure

    Returns:
        tuple(pd.DataFrame, pd.DataFrame): DNM and rates dataframes with records with missing scores removed
    """

    logger = logging.getLogger("logger")

    set_plain_log()
    logger.info("=" * 50)
    logger.info("[MISSING SCORES]")
    logger.info("=" * 50)
    set_regular_log()

    # Remove locus without scores in both DNM and rates dataframes
    dnm_kept_df = dnm_df.loc[~dnm_df[score_column].isna()].copy()
    rates_kept_df = rates_df.loc[~rates_df[score_column].isna()].copy()

    if dnm_df.shape[0] != dnm_kept_df.shape[0]:
        logger.info(
            f"{dnm_df.shape[0] - dnm_kept_df.shape[0] } observed DNMs were removed as they do not have a {score_column} score"
        )

    if rates_df.shape[0] != rates_kept_df.shape[0]:
        logger.info(
            f"{rates_df.shape[0]  - rates_kept_df.shape[0]} variants were removed from the rates file as they do not have a {score_column} score"
        )

    # Log number of removed variants per gene
    dnm_discarded_df = dnm_df.loc[dnm_df[score_column].isna()].copy()
    rates_discarded_df = rates_df.loc[rates_df[score_column].isna()].copy()

    for gene_id, gene_dnm_df in dnm_discarded_df.groupby("gene_id"):
        prep_logs[gene_id]["nb_observed_dnms_discarded"] = gene_dnm_df.shape[0]

    for gene_id, gene_rates_df in rates_discarded_df.groupby("gene_id"):
        prep_logs[gene_id]["nb_rates_variants_discarded"] = gene_rates_df.shape[0]

    return prep_logs, dnm_kept_df, rates_kept_df
