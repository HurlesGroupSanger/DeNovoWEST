import pandas as pd
import numpy as np
import itertools
import os
import click
import re
from rpy2.robjects import Formula, FloatVector
from rpy2.robjects.packages import importr
import rpy2.robjects as ro

SHET_BINS = ["shethigh == True", "shethigh == False"]
MCR_BINS = ["constrained == True", "constrained == False"]


@click.command()
@click.argument("count_table")
@click.argument("outfile")
def main(count_table: str, outfile: str):
    """

    Args:
        count_table (str): Path to a table listing for each bin of variants the ratio observed/expected
        outfile (str): Path to the enrichment results
    """

    count_table = pd.read_csv(count_table, sep="\t")
    enrichment_df = adjust_enrichment(count_table)
    enrichment_df.to_csv(outfile, sep="\t")


def adjust_enrichment(count_table: pd.DataFrame) -> pd.DataFrame:

    bins = define_bins()

    # fit lowess to missense and nonsense
    df_lowess = pd.concat(
        [fit_loess_missense_nonsense(count_table, condition) for condition in bins["lowess"]],
        axis=0,
        ignore_index=True,
    )

    # set values below 1 to 1
    df_lowess.loc[df_lowess.obs_exp < 1, "obs_exp"] = 1

    # add frameshifts and inframe
    df_lowess = add_frameshifts_inframe(df_lowess, count_table, bins["inframe"], bins["frameshift"])

    # concatenate lowess with count table
    results = pd.concat(
        [count_table[~count_table.bin.str.contains("nonsense|missense")], df_lowess], axis=0, ignore_index=True
    )[["obs_exp", "score", "bin"]]
    results.score = results.score if results.bin.str.contains("nonsense|missense").any() else np.nan

    # get ppv
    return positive_predictive_value(results)


def define_bins() -> list:

    bins = dict()
    bins["inframe"] = define_inframe_bins()
    bins["frameshift"] = define_frameshift_bins()
    bins["lowess"] = define_lowess_bins()

    return bins


def define_inframe_bins():

    inframe_bins = [
        " and ".join(list(x)) for x in itertools.product(SHET_BINS, ['consequence.str.contains("missense")'])
    ]

    return inframe_bins


def define_frameshift_bins():

    frameshift_bins = list(itertools.product(SHET_BINS, ['consequence.str.contains("nonsense")']))
    return frameshift_bins


def define_lowess_bins():

    lowess_nonsenses = list(itertools.product(SHET_BINS, ["nonsense"]))
    lowess_missense = list(itertools.product(SHET_BINS, MCR_BINS, ["missense"]))

    return lowess_nonsenses + lowess_missense


def fit_loess_missense_nonsense(df, condition) -> pd.DataFrame:

    try:
        if "nonsense" in condition[1]:
            df = df[df["bin"].str.contains(condition[0]) & df["bin"].str.contains("nonsense")]
            condition = " and ".join(condition).replace("nonsense", 'cq.str.contains("nonsense")')
        else:
            df = df[
                df["bin"].str.contains(condition[0])
                & df["bin"].str.contains(condition[1])
                & df["bin"].str.contains("missense")
            ]
            condition = " and ".join(condition).replace("missense", 'cq.str.contains("missense")')

        df["midpoint"] = df["bin"].apply(lambda x: get_CADD_midpoint(x))

        vals = np.arange(0, 52.5, 0.001)
        df = fit_loess(df, x="midpoint", y="obs_exp", w="obs", new_vals=vals, condition=condition)

        return df

    except:
        pass


def get_CADD_midpoint(bin_name):
    """_summary_

    Args:
        bin_name (_type_): _description_

    Returns:
        _type_: _description_
    """
    minimum = 3 if "missense" in bin_name else 3.75
    cadd_lower_bound = [float(x) for x in re.findall(r"\d+(?:\.\d+)?", bin_name)][0]
    return cadd_lower_bound + minimum


def fit_loess(df, x, y, w, new_vals, span=1, condition=str):

    rstats = importr("stats")

    df = df.sort_values(x).query("obs_exp > 0")

    x, y, w, vals = map(FloatVector, [df[x], df[y], df[w], new_vals])

    fmla = Formula("y ~ x")
    env = fmla.environment
    env["x"] = x
    env["y"] = y
    env["w"] = w
    env["vals"] = vals

    model = rstats.loess(fmla, weights=w, span=span)
    # fitted = np.array(model.rx2('fitted'))

    prediction = ro.r["predict"](model, vals)

    new_df = pd.DataFrame({"obs_exp": prediction, "bin": [condition] * len(prediction), "score": new_vals.round(3)})
    new_df["bin"] = new_df["bin"] + " and score == " + new_df["score"].astype(str)

    min_predicted, max_predicted = new_df.obs_exp.min(), new_df.obs_exp.max()
    min_score, max_score = (
        new_df[new_df.obs_exp == min_predicted].score.min(),
        new_df[new_df.obs_exp == max_predicted].score.max(),
    )
    new_df["obs_exp"] = np.select(
        [new_df.score < min_score, new_df.score > max_score], [min_predicted, max_predicted], default=new_df.obs_exp
    )

    return new_df


def add_frameshifts_inframe(df, count_table, inframe_bins, frameshift_bins) -> pd.DataFrame:

    df_original = df.copy()

    # frameshifts
    for condition in frameshift_bins:
        df_temp = df[(df["bin"].str.contains(condition[0])) & (df["bin"].str.contains("nonsense"))]
        max_obs_exp = df_temp.obs_exp.max()
        condition = " and ".join(condition).replace("nonsense", "frameshift")
        df = df.append({"obs_exp": max_obs_exp, "bin": condition, "score": np.nan}, ignore_index=True)

    # inframe
    df_temp = count_table[count_table["bin"].isin(inframe_bins)]
    df_temp["bin"] = df_temp["bin"].str.replace("missense", "inframe")

    df = df.append(df_temp, ignore_index=True).query('bin.str.contains("frameshift") | bin.str.contains("inframe")')

    return pd.concat([df_original, df], axis=0, ignore_index=True)


def positive_predictive_value(df) -> pd.DataFrame:

    odds_ratio = df.obs_exp - df.obs_exp.min() + 1
    df["ppv"] = np.select([df.bin.str.contains("synonymous")], [0.001], default=(odds_ratio - 1) / odds_ratio)

    return df


if __name__ == "__main__":
    main()
