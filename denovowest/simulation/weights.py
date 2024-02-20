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


def assign_weights(df, column, indel_weights=pd.DataFrame()):
    """
    Assign weights to each variant in the DNM or rates dataframe

    Args:
        df (pd.DataFrame): DNM or rates dataframe
        weights_df (pd.DataFrame): weights dataframe
    """

    weighted_df = min_max_transformation(df, column)
    if not indel_weights.empty:
        weighted_df = assign_indel_weights(weighted_df, indel_weights)
    return weighted_df


def min_max_transformation(df, column):
    max_val = df[column].max()
    min_val = df[column].min()

    max_val = 100
    min_val = 0

    df["ppv"] = df[column].apply(lambda x: (x - min_val) / (max_val - min_val))

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

    Args:
        weights_df (pd.DataFrame): weights dataframe
    """

    # Here I picked weights obtained with classic approach on 31k cohort
    indel_weights = pd.DataFrame(
        [
            ["inframe", False, 0.1594364909515928],
            ["inframe", True, 0.5069925869911193],
            ["frameshift", False, 0.35831160882134394],
            ["frameshift", True, 0.8825922897781597],
        ],
        columns=["consequence", "shethigh", "ppv"],
    )

    return indel_weights
