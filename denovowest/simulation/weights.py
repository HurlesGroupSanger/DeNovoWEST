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


def assign_weights(df, column, min_score, max_score, runtype, indel_weights=pd.DataFrame()):
    """
    Assign weights to each variant in the DNM or rates dataframe

    Args:
        df (pd.DataFrame): DNM or rates dataframe
        weights_df (pd.DataFrame): weights dataframe
    """

    # Map the scores between 0 and 1
    weighted_df = min_max_transformation(df, column, min_score, max_score)

    if runtype == "ns":
        # Some splice_lof variants are not scored by most CEPs, we assign them a weight based on previous runs in a similar fashion as for indels
        weighted_df = fix_splice_lof(df, column)

        # We do not assign any indel weights for variants in the rates file as there are no indels in it
        if not indel_weights.empty:
            weighted_df = assign_indel_weights(weighted_df, indel_weights)
        return weighted_df


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


def fix_splice_lof(df, column):
    """
    Add a score to splice_lof variants missing one

    Args:
        df (pd.DataFrame): DNM or rates dataframe
        column (str): Column containing the CEP score to use

    Returns:
        pd.DataFrame: DNM or rates with imputed missing splice_lof scores

    """

    SHET_HIGH_SPLICELOF_WEIGHT = 0.779862
    SHET_LOW_SPLICELOF_WEIGHT = 0.309483

    df.loc[(df["shethigh"] is True) & (df.consequence == "splice_lof"), column] = df.loc[
        (df["shethigh"] is True) & (df.consequence == "splice_lof"), column
    ].fillna(SHET_HIGH_SPLICELOF_WEIGHT)
    df.loc[(df["shethigh"] is False) & (df.consequence == "splice_lof"), column] = df.loc[
        (df["shethigh"] is False) & (df.consequence == "splice_lof"), column
    ].fillna(SHET_LOW_SPLICELOF_WEIGHT)

    return df


def assign_indel_weights(df, indel_weights):
    """
    Assigns weights to indels in a DataFrame based on the provided indel_weights.

    Args:
        df (pandas.DataFrame): The DataFrame containing indel data.
        indel_weights (pandas.DataFrame): The DataFrame containing indel weights.

    Returns:
        pandas.DataFrame: The DataFrame with assigned indel weights.
    """

    min_df = df.loc[df.consequence.isin(["frameshift", "inframe"])]
    df2 = min_df.drop("ppv", axis=1).merge(indel_weights, on=["consequence", "shethigh"])
    df2.index = min_df.index
    df.update(df2)

    return df


def get_indel_weights():
    """
    Returns a dataframe containing only indel weights

    """

    # Here I picked weights obtained with my b38 BayesDel run
    indel_weights = pd.DataFrame(
        [
            ["inframe", False, 0.159436],
            ["inframe", True, 0.506993],
            ["frameshift", False, 0.530357],
            ["frameshift", True, 0.907829],
        ],
        columns=["consequence", "shethigh", "ppv"],
    )

    return indel_weights
