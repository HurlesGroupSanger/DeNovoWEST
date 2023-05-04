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
    weighted_df = df.merge(
        weights_df,
        left_on=["consequence", "score_rounded", "constrained", "shethigh"],
        right_on=["consequence", "score", "constrained", "shethigh"],
        how="left",
    )
    assert df.shape[0] == weighted_df.shape[0]

    # Clean merging residuals
    weighted_df.drop(["score_y", "score_rounded"], axis=1, inplace=True)
    weighted_df.rename({"score_x": "score"}, axis=1, inplace=True)

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
        max_score = get_max_score(combination, df)
        df["ppv"] = np.where(
            (df.consequence == combination[0])
            & (df.score > combination[1])
            & (df.constrained == combination[2])
            & (df.shethigh == combination[3]),
            max_score,
            df["ppv"],
        )


def get_max_score(combination, df):
    """_summary_

    Args:
        combination (_type_): _description_
        df (_type_): _description_
    """

    min_df = df.loc[
        (df.consequence == combination[0])
        & (df.score > combination[1])
        & (df.constrained == combination[2])
        & (df.shethigh == combination[3])
    ]

    max_score = min_df.score.max()

    return max_score


def get_indel_weights(weights_df):
    """_summary_

    Args:
        weights_df (_type_): _description_
    """
    indel_weights = weights_df.loc[weights_df.consequence.str.contains("frame")]

    return indel_weights
