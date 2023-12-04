"""
@queenjobo @ksamocha

16/08/2019

Functions to assign weights to different variant classes

"""

import pandas as pd
import logging as lg
import numpy as np
import itertools
from scipy import stats
import sys


def assign_weights(df, weights_df):
    """
    Assign weights to each variant in the DNM or rates dataframe

    Args:
        df (pd.DataFrame): DNM or rates dataframe
        weights_df (pd.DataFrame): weights dataframe
    """

    weighted_df = assign_standard_weights(df, weights_df)
    weighted_df = assign_outlier_weights(weighted_df, weights_df)

    return weighted_df


def assign_standard_weights(df, weights_df):
    """
    Assign standard weights to each variant in the DNM or rates dataframe.
    By standard, we mean variants that can be directly assigned a weight from the weights dataframe.

    Args:
        df (pd.DataFrame): DNM or rates dataframe
        weights_df (pd.DataFrame): weights dataframe
    """

    # Variant scores have several decimals, but the weights bin are limited to 3 decimals
    df["score_rounded"] = df.score.round(3)

    # Assign weights to each variant depending on the variant type
    missense_weighted_df = assign_weights_missense(df, weights_df)
    nonsense_weighted_df = assign_weights_nonsense(df, weights_df)
    synonymous_weighted_df = assign_weights_synonymous(df, weights_df)
    splicelof_weighted_df = assign_weights_splicelof(df, weights_df)
    inframe_weighted_df = assign_weights_inframe(df, weights_df)
    frameshift_weighted_df = assign_weights_frameshift(df, weights_df)

    weighted_df = pd.concat(
        [
            missense_weighted_df,
            nonsense_weighted_df,
            synonymous_weighted_df,
            splicelof_weighted_df,
            inframe_weighted_df,
            frameshift_weighted_df,
        ]
    )

    # TODO : Further check here (splice_region not handled)
    assert weighted_df.shape[0] <= df.shape[0]

    # Clean merging residuals
    weighted_df.drop(["score_y", "score_rounded"], axis=1, inplace=True)
    weighted_df.rename({"score_x": "score"}, axis=1, inplace=True)

    return weighted_df


def assign_weights_missense(df, weights_df):
    """
    Assign weights to each standard missense variant in the DNM or rates dataframe.
    By standard, we mean variants having a score corresponding to a bin in the weights dataframe.

    Args:
        df (pd.DataFrame): DNM or rates dataframe
        weights_df (pd.DataFrame): weights dataframe
    """

    min_df = df.loc[df.consequence == "missense"]
    weighted_df = min_df.merge(
        weights_df,
        left_on=["consequence", "constrained", "shethigh", "score_rounded"],
        right_on=["consequence", "constrained", "shethigh", "score"],
        how="left",
    )
    assert weighted_df.shape[0] == min_df.shape[0]

    return weighted_df


def assign_weights_nonsense(df, weights_df):
    """
    Assign weights to each standard nonsense variant in the DNM or rates dataframe.
    By standard, we mean variants having a score corresponding to a bin in the weights dataframe.

    Args:
        df (pd.DataFrame): DNM or rates dataframe
        weights_df (pd.DataFrame): weights dataframe
    """

    min_df = df.loc[df.consequence == "nonsense"]

    # TODO : During the annotation steps, all nonsense variants are assigned to False for the constrained column
    # When rewriting the weights module I set the constrained column to NA for nonsense variants
    # This is because the constrained column is not used to compute enrichments for nonsense variants
    # In order to assign weights to nonsense variants, we need to set constrained column to False
    # This would need to be better handled in the future (using NA for example)
    weights_df_copy = weights_df.copy()
    weights_df_copy.constrained = False
    weighted_df = min_df.merge(
        weights_df_copy,
        left_on=["consequence", "constrained", "shethigh", "score_rounded"],
        right_on=["consequence", "constrained", "shethigh", "score"],
        how="left",
    )
    assert weighted_df.shape[0] == min_df.shape[0]

    return weighted_df


def assign_weights_synonymous(df, weights_df):
    """
    Assign weights to each synonymous variant in the DNM or rates dataframe.

    Args:
        df (pd.DataFrame): DNM or rates dataframe
        weights_df (pd.DataFrame): weights dataframe
    """

    min_df = df.loc[df.consequence == "synonymous"]

    # TODO : see assign_weights_nonsense
    weights_df_copy = weights_df.copy()
    weights_df_copy.constrained = False
    weighted_df = min_df.merge(
        weights_df_copy,
        left_on=["consequence", "constrained", "shethigh"],
        right_on=["consequence", "constrained", "shethigh"],
        how="left",
    )
    assert weighted_df.shape[0] == min_df.shape[0]

    return weighted_df


def assign_weights_splicelof(df, weights_df):
    """
    Assign weights to each splice_lof variant in the DNM or rates dataframe.

    Args:
        df (pd.DataFrame): DNM or rates dataframe
        weights_df (pd.DataFrame): weights dataframe
    """

    min_df = df.loc[df.consequence == "splice_lof"]

    # TODO : see assign_weights_nonsense
    weights_df_copy = weights_df.copy()
    weights_df_copy.constrained = False
    weighted_df = min_df.merge(
        weights_df_copy,
        left_on=["consequence", "constrained", "shethigh"],
        right_on=["consequence", "constrained", "shethigh"],
        how="left",
    )
    assert weighted_df.shape[0] == min_df.shape[0]

    return weighted_df


def assign_weights_inframe(df, weights_df):
    """
    Assign weights to each inframe variant in the DNM or rates dataframe.

    Args:
        df (pd.DataFrame): DNM or rates dataframe
        weights_df (pd.DataFrame): weights dataframe
    """

    min_df = df.loc[df.consequence == "inframe"]

    # TODO : see assign_weights_nonsense
    weights_df_copy = weights_df.copy()
    weights_df_copy.constrained = False
    weighted_df = min_df.merge(
        weights_df_copy,
        left_on=["consequence", "constrained", "shethigh"],
        right_on=["consequence", "constrained", "shethigh"],
        how="left",
    )
    assert weighted_df.shape[0] == min_df.shape[0]

    return weighted_df


def assign_weights_frameshift(df, weights_df):
    """
    Assign weights to each frameshift variant in the DNM or rates dataframe.

    Args:
        df (pd.DataFrame): DNM or rates dataframe
        weights_df (pd.DataFrame): weights dataframe
    """

    min_df = df.loc[df.consequence == "frameshift"]

    # TODO : see assign_weights_nonsense
    weights_df_copy = weights_df.copy()
    weights_df_copy.constrained = False
    weighted_df = min_df.merge(
        weights_df_copy,
        left_on=["consequence", "constrained", "shethigh"],
        right_on=["consequence", "constrained", "shethigh"],
        how="left",
    )
    assert weighted_df.shape[0] == min_df.shape[0]

    return weighted_df


def assign_outlier_weights(df, weights_df):
    """
    Assign weights to each missense or nonsense outlier variant in the DNM or rates dataframe.
    By outlier, we mean variants that cannot be directly assigned a weight from the weights dataframe, that are variant
    with a CADD score higher than the maximum score in the weights dataframe.

    Args:
        df (pd.DataFrame): DNM or rates dataframe
        weights_df (pd.DataFrame): weights dataframe
    """

    # Get the maximum CADD score in the weights dataframe
    max_score_bin = weights_df.score.max()

    consequence = ["missense", "nonsense"]
    score = [max_score_bin]
    constrained = [True, False]
    shethigh = [True, False]

    combinations = itertools.product(consequence, score, constrained, shethigh)

    for combination in combinations:
        # There is no such thing as a nonsense variant that is constrained (see assign_weights_nonsense)
        if combination[0] == "nonsense" and combination[2] == True:
            continue

        # Retrieve the maximum PPV for the given combination of variant characteristics
        max_score = get_max_score(combination, df)

        # Assign the maximum PPV to the outlier variants (variants with a CADD score higher than the maximum CADD score)
        df["ppv"] = np.where(
            (df.consequence == combination[0])
            & (df.score > combination[1])
            & (df.constrained == combination[2])
            & (df.shethigh == combination[3]),
            max_score,
            df["ppv"],
        )

    return df


def get_max_score(combination, df):
    """
    For a given category of variants, retrieve the maximum PPV in the DNM or rates dataframe

    Args:
        combination (list): list of variant characteristics defining a category of variants
        df (pd.DataFrame): DNM or rates dataframe
    """

    min_df = df.loc[
        (df.consequence == combination[0]) & (df.constrained == combination[2]) & (df.shethigh == combination[3])
    ]

    max_score = min_df.ppv.max()

    return max_score


def get_indel_weights(weights_df):
    """
    Returns a dataframe containing only indel weights

    Args:
        weights_df (pd.DataFrame): weights dataframe
    """
    indel_weights = weights_df.loc[weights_df.consequence.str.contains("frame")]

    return indel_weights
