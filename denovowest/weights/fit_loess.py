#!/usr/bin/env python

import pandas as pd
import numpy as np
import itertools
import os
import click
import re
import utils

from rpy2.robjects import Formula, FloatVector
from rpy2.robjects.packages import importr
import rpy2.robjects as ro


SHET_VALUES = [True, False]
CONSTRAINED_VALUES = [True, False]
CSQ_MISSENSE_VALUES = ["missense"]
CSQ_NONSENSE_VALUES = ["nonsense"]

CADD_MIN = 0
CADD_MAX = 52.5
CADD_STEP = 0.001
# CADD_STEP = 0.1


@click.command()
@click.argument("obs_exp_table_path")
@click.argument("outfile")
@click.option("--outdir", default="enrichment_plots")
def main(obs_exp_table_path: str, outfile: str, outdir: str):
    """
    Fits a loess regression using CADD score bins midpoint to infer expected and observed counts for every intermediary CADD score.

    Args:
        obs_exp_table_path (str): Path to a table listing for each bin of variants the ratio observed/expected
        outfile (str): Path to the enrichment results
    """

    obs_exp_table = pd.read_csv(obs_exp_table_path, sep="\t")

    # Builds the meta categories of mutations (not considering CADD score)
    categories = define_categories()

    # New CADD scores we want to infer the observed/expected ratio on
    # TODO deport all code dedicated to CADD loess in separate function
    new_cadd_scores = np.arange(CADD_MIN, CADD_MAX + CADD_STEP, CADD_STEP)

    # For each category we want to run loess on
    list_df_lowess = list()
    for category in categories["lowess"]:
        # Extract bins corresponding to each category
        min_obs_exp_table = prepare_for_loess(obs_exp_table, category)
        if min_obs_exp_table.empty:
            print("No observed data for category ")
            continue

        # Run the loess regression
        # TODO : find out about the different use of span parameter and its sensitivity
        if category["consequence"] == "nonsense":
            span = 1
        elif category["consequence"] == "missense":
            span = 0.99

        loess_prediction = fit_loess(
            min_obs_exp_table, x="midpoint", y="obs_exp", w="obs", new_cadd_scores=new_cadd_scores, span=span
        )

        obs_exp_cadd_df = get_obs_exp_ratio_per_score(min_obs_exp_table, loess_prediction, new_cadd_scores)
        list_df_lowess.append(obs_exp_cadd_df)

    # Concatenate dataframes computed on each category
    df_lowess = pd.concat(list_df_lowess, axis=0, ignore_index=True)

    df_lowess.constrained = df_lowess.constrained.astype("boolean")

    # Compute frameshift and inframe variants observed/expected ratio per category, without any CADD score consideration
    # df_meta = get_obs_exp_ratio(obs_exp_table, categories["inframe"], categories["frameshift"])
    df_frameshift = get_obs_exp_ratio_frameshift(df_lowess, categories["frameshift"])
    df_inframe = get_obs_exp_ratio_inframe(obs_exp_table, categories["inframe"])

    # Compute for other variants (synonymous, splice regions)
    df_other = get_obs_exp_ratio_other(obs_exp_table)

    # Concatenate all indepents dataframes
    res_df = pd.concat([df_lowess, df_frameshift, df_inframe, df_other], axis=0, ignore_index=True)
    # Plot the enrichment curves before setting the minimum observed/expected ratio to 1
    utils.plot_enrichment(obs_exp_table, mode="obs_exp", loess=True, df_loess=res_df, outdir=outdir)
    df_lowess.loc[df_lowess.obs_exp < 1, "obs_exp"] = 1  # TODO : Why ?
    res_df = pd.concat([df_lowess, df_frameshift, df_inframe, df_other], axis=0, ignore_index=True)

    # Turn ratio into PPVs
    res_df = positive_predictive_value(res_df)

    # Plot enrichment both at the observed/expected ratio level and ppv level
    utils.plot_enrichment(obs_exp_table, mode="ppv", loess=True, df_loess=res_df, outdir=outdir)

    # Ad there are some missing values in constrained columns, pandas turn boolean values to float automatically. We force boolean use here.
    # res_df.constrained = res_df.constrained.astype("boolean")

    # Export results
    res_df.drop(["obs", "exp"], axis=1, inplace=True)
    if "." in str(CADD_STEP):
        res_df.score = res_df.score.round(len(str(CADD_STEP)) - str(CADD_STEP).index(".") - 1)
    res_df.to_csv(outfile, sep="\t", index=False, na_rep="NA")


def prepare_for_loess(obs_exp_table, category):
    """
    Retrieve all the CADD range bins corresponding to that category

    Args:
        obs_exp_table (pd.DataFrame) : Table listing for each bin of variants the ratio observed/expected
        category (dict): Description of the meta category we want to gather bins in

    Returns:
        pd.DataFrame : Subset of obs_exp_table containing all bins falling in the input category

    """

    # Retrieve all the CADD range bins corresponding to that category
    min_obs_exp_table = extract_bins_in_category(obs_exp_table, category)
    min_obs_exp_table = min_obs_exp_table.loc[~min_obs_exp_table.score.isna()]

    # Get the midpoint score of each CADD range bin
    min_obs_exp_table["midpoint"] = min_obs_exp_table["score"].apply(
        lambda x: get_CADD_midpoint(x, category["consequence"])
    )

    # Sort table by ascending CADD midpoint, and keep only bins for which we observed at least one mutation
    min_obs_exp_table = min_obs_exp_table.sort_values("midpoint").query("obs_exp > 0")

    return min_obs_exp_table


def extract_bins_in_category(obs_exp_table, category):
    """
    Retrieve all the CADD range bins corresponding to that category

    Args:
        obs_exp_table (pd.DataFrame) : Table listing for each bin of variants the ratio observed/expected
        category (dict): Description of the meta category we want to gather bins in

    Returns:
        pd.DataFrame : Subset of obs_exp_table containing all bins falling in the input category

    """

    min_obs_exp_table = obs_exp_table
    for key, value in category.items():
        min_obs_exp_table = min_obs_exp_table.loc[min_obs_exp_table[key] == value]

    return min_obs_exp_table


def define_categories() -> list:
    """
    Builds the meta categories (not including CADD scores) used to categorize mutations.

    Returns:
        list: list of categories
    """

    categories = dict()

    categories["inframe"] = define_inframe_categories()
    categories["frameshift"] = define_frameshift_categories()
    categories["lowess"] = define_lowess_categories()

    return categories


def define_inframe_categories():
    """
    Builds the categories containing inframe variants.

    Returns:
        list: list of categories
    """
    inframe_categories = list(itertools.product(SHET_VALUES, CSQ_MISSENSE_VALUES))

    return categories_tuple_to_dict(inframe_categories, ["shethigh", "consequence"])


def define_frameshift_categories():
    """
    Builds the categories containing frameshift variants.

    Returns:
        list: list of categories
    """
    frameshift_categories = list(itertools.product(SHET_VALUES, CSQ_NONSENSE_VALUES))
    return categories_tuple_to_dict(frameshift_categories, ["shethigh", "consequence"])


def define_lowess_categories():
    """
    Builds the categories containing inframe variants.

    Returns:
        list: list of categories
    """
    lowess_nonsenses = list(itertools.product(SHET_VALUES, CSQ_NONSENSE_VALUES))
    lowess_missense = list(itertools.product(SHET_VALUES, CONSTRAINED_VALUES, CSQ_MISSENSE_VALUES))

    return categories_tuple_to_dict(lowess_nonsenses, ["shethigh", "consequence"]) + categories_tuple_to_dict(
        lowess_missense, ["shethigh", "constrained", "consequence"]
    )


def categories_tuple_to_dict(categories, keys):
    categories_dict = [dict(zip(keys, category)) for category in categories]
    return categories_dict


def get_CADD_midpoint(cadd_range, consequence):
    """Returns the midpoint for each CADD bin depending on its category

    Args:
        cadd_range (str): lower and upper bounds of the CADD bin
        consequence (str) : functional consequence

    Returns:
        float: midpoint of the CADD bin
    """

    # TODO : handle hardcoded upper bins properly
    if "+" in cadd_range:
        cadd_lb = float(cadd_range.replace("+", ""))
        if consequence == "missense":
            midpoint = cadd_lb + 3
        if consequence == "nonsense":
            midpoint = cadd_lb + 3.75

    else:
        cadd_lb = float(cadd_range.split("-")[0])
        if consequence == "missense":
            midpoint = cadd_lb + 3
        if consequence == "nonsense":
            if cadd_lb == 0:
                midpoint = 15
            else:
                midpoint = cadd_lb + 3.75

    return midpoint


def fit_loess(obs_exp_table, x, y, w, new_cadd_scores, span=0.99):
    """
    Fit a loess regression between CADD score bins midpoints and the observed to expected mutations ratio.

    Args:
        obs_exp_table (pd.DataFrame): Table listing for each category of variants the ratio observed/expected
        x (list): CADD bins midpoints scores (e.g. 7.5 if CADD bin represents variants with a CADD score between 5 and 10])
        y (list): Observed/expected mutation ratio for each bin
        w (list): Weights used in the regression model, number of observed mutations for each bin
        new_cadd_scores (list): New "continuous" CADD scores for which we will infer the observed/expected ratio thanks to loess
        span (int, optional): TODO. Defaults to 1.


    Returns:
        pd.DataFrame: Table with observed/expected ratio for each continuous CADD score
    """
    rstats = importr("stats")

    # Prepare model
    x, y, w, vals = map(FloatVector, [obs_exp_table[x], obs_exp_table[y], obs_exp_table[w], new_cadd_scores])

    fmla = Formula("y ~ x")
    env = fmla.environment
    env["x"] = x
    env["y"] = y
    env["w"] = w
    env["vals"] = vals

    # Run loess
    model = rstats.loess(fmla, weights=w, span=span)
    # model_dict = dict(zip(model.names, list(model)))
    # print(x)
    # print(y)
    # print(w)
    # print(model_dict.keys())
    # print(model_dict["fitted"])

    prediction = ro.r["predict"](model, vals)

    return prediction


def get_obs_exp_ratio_per_score(obs_exp_table, loess_prediction, new_cadd_scores):
    """_summary_

    Args:
        obs_exp_table (_type_): _description_
        loess_prediction (_type_): _description_
        new_cadd_scores (_type_): _description_

    Returns:
        _type_: _description_
    """
    # Build the observed/expected ratio table for each intermediate CADD score
    res_df = pd.DataFrame(
        {
            "consequence": [obs_exp_table["consequence"].iloc[0]] * len(loess_prediction),
            "score": new_cadd_scores,
            "constrained": [obs_exp_table["constrained"].iloc[0]] * len(loess_prediction),
            "shethigh": [obs_exp_table["shethigh"].iloc[0]] * len(loess_prediction),
            "obs_exp": loess_prediction,
        }
    )

    # Get the minimum and maximum observed/expected ratio across all bins
    min_predicted, max_predicted = res_df.obs_exp.min(), res_df.obs_exp.max()

    # Get the CADD score that corresponds
    min_score, max_score = (
        res_df[res_df.obs_exp == min_predicted].score.min(),
        res_df[res_df.obs_exp == max_predicted].score.max(),
    )
    res_df["obs_exp"] = np.select([res_df.score < min_score], [min_predicted], default=res_df.obs_exp)
    res_df["obs_exp"] = np.select([res_df.score > max_score], [max_predicted], default=res_df.obs_exp)

    return res_df


def get_obs_exp_ratio_frameshift(obs_exp_table, frameshift_categories):
    """TODO"""

    df = pd.DataFrame(columns=obs_exp_table.columns)
    for category in frameshift_categories:
        min_obs_exp_table = extract_bins_in_category(obs_exp_table, category)
        max_obs_exp = min_obs_exp_table.obs_exp.max()

        df = pd.concat(
            [
                df,
                pd.DataFrame(
                    [
                        {
                            "consequence": "frameshift",
                            "score": np.nan,
                            "constrained": np.nan,
                            "shethigh": min_obs_exp_table["shethigh"].iloc[0],
                            "obs_exp": max_obs_exp,
                        }
                    ]
                ),
            ]
        )

    return df


def get_obs_exp_ratio_inframe(obs_exp_table, inframe_categories):
    """TODO"""

    df = pd.DataFrame(columns=obs_exp_table.columns)
    for category in inframe_categories:
        min_obs_exp_table = extract_bins_in_category(obs_exp_table, category)
        min_obs_exp_table = min_obs_exp_table.loc[min_obs_exp_table.score.isna()]

        max_obs_exp = min_obs_exp_table.obs_exp.max()

        df = pd.concat(
            [
                df,
                pd.DataFrame(
                    [
                        {
                            "consequence": "inframe",
                            "score": np.nan,
                            "constrained": np.nan,
                            "shethigh": min_obs_exp_table["shethigh"].iloc[0],
                            "obs_exp": max_obs_exp,
                        }
                    ]
                ),
            ]
        )

    return df


def get_obs_exp_ratio_other(obs_exp_table):
    """ """

    df = obs_exp_table.loc[(~obs_exp_table["consequence"].str.contains("nonsense|missense"))]
    return df


def positive_predictive_value(df) -> pd.DataFrame:
    """TODO

    Args:
        df (_type_): _description_

    Returns:
        pd.DataFrame: _description_
    """
    odds_ratio = df.obs_exp - df.obs_exp.min() + 1
    df["ppv"] = np.select(
        [df["consequence"].str.contains("synonymous")], [0.001], default=(odds_ratio - 1) / odds_ratio
    )

    return df


if __name__ == "__main__":
    main()
