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
    """#TODO

    Args:
        df (pd.DataFrame): DNM or mutation df table
        weights_df (pd.DataFrame): Weights
    """

    weighted_df = assign_standard_weights(df, weights_df)
    weighted_df = assign_outlier_weights(weighted_df, weights_df)

    return weighted_df


def assign_standard_weights(df, weights_df):
    """#TODO

    Args:
        df (pd.DataFrame): DNM or mutation df table
        weights_df (pd.DataFrame): Weights
    """

    # Variant scores have several decimals, but the weights bin are limited to 3 decimals
    df["score_rounded"] = df.score.round(3)

    # Merge on consequence, score, constraints and shethigh in order to retrieve the weights to assign to each variant
    # weighted_df = df.merge(
    #     weights_df,
    #     left_on=["consequence", "#", "constrained", "shethigh"],
    #     right_on=["consequence", "score", "constrained", "shethigh"],
    #     how="left",
    # )
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
    """_summary_

    Args:
        df (_type_): _description_
        weights_df (_type_): _description_
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
    """_summary_

    Args:
        df (_type_): _description_
        weights_df (_type_): _description_
    """

    min_df = df.loc[df.consequence == "nonsense"]
    # TODO : Explain why we are not setting the constraint information to nan for nonsense mutations
    # min_df.constrained = np.NaN
    # min_df.constrained = min_df.constrained.astype(object)
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
    """_summary_

    Args:
        df (_type_): _description_
        weights_df (_type_): _description_
    """

    min_df = df.loc[df.consequence == "synonymous"]
    weighted_df = min_df.merge(
        weights_df,
        left_on=["consequence", "constrained", "shethigh"],
        right_on=["consequence", "constrained", "shethigh"],
        how="left",
    )
    assert weighted_df.shape[0] == min_df.shape[0]

    return weighted_df


def assign_weights_splicelof(df, weights_df):
    """_summary_

    Args:
        df (_type_): _description_
        weights_df (_type_): _description_
    """

    min_df = df.loc[df.consequence == "splice_lof"]
    weighted_df = min_df.merge(
        weights_df,
        left_on=["consequence", "constrained", "shethigh"],
        right_on=["consequence", "constrained", "shethigh"],
        how="left",
    )
    assert weighted_df.shape[0] == min_df.shape[0]

    return weighted_df


def assign_weights_inframe(df, weights_df):
    """_summary_

    Args:
        df (_type_): _description_
        weights_df (_type_): _description_
    """

    min_df = df.loc[df.consequence == "inframe"]
    weighted_df = min_df.merge(
        weights_df,
        left_on=["consequence", "constrained", "shethigh"],
        right_on=["consequence", "constrained", "shethigh"],
        how="left",
    )
    assert weighted_df.shape[0] == min_df.shape[0]

    return weighted_df


def assign_weights_frameshift(df, weights_df):
    """_summary_

    Args:
        df (_type_): _description_
        weights_df (_type_): _description_
    """

    min_df = df.loc[df.consequence == "frameshift"]
    weighted_df = min_df.merge(
        weights_df,
        left_on=["consequence", "constrained", "shethigh"],
        right_on=["consequence", "constrained", "shethigh"],
        how="left",
    )
    assert weighted_df.shape[0] == min_df.shape[0]

    return weighted_df


def assign_outlier_weights(df, weights_df):
    """#TODO

    Args:
        df (pd.DataFrame): DNM or mutation df table
        weights_df (pd.DataFrame): Weights
    """
    max_score_bin = weights_df.score.max()

    consequence = ["missense", "nonsense"]
    score = [max_score_bin]
    constrained = [True, False]
    shethigh = [True, False]

    combinations = itertools.product(consequence, score, constrained, shethigh)

    for combination in combinations:
        if combination[0] == "nonsense" and combination[2] == True:
            continue
        max_score = get_max_score(combination, df)

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
    """_summary_

    Args:
        combination (_type_): _description_
        df (_type_): _description_
    """

    min_df = df.loc[
        (df.consequence == combination[0]) & (df.constrained == combination[2]) & (df.shethigh == combination[3])
    ]

    max_score = min_df.ppv.max()

    return max_score


def get_indel_weights(weights_df):
    """_summary_

    Args:
        weights_df (_type_): _description_
    """
    indel_weights = weights_df.loc[weights_df.consequence.str.contains("frame")]

    return indel_weights
