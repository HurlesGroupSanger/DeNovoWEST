import logging.config
import logging
import pandas as pd
import seaborn as sns
import os
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import poisson


CONSEQUENCES_MAPPING = {
    "frameshift_variant": "frameshift",
    "frameshift": "frameshift",
    "inframe_insertion": "inframe",
    "inframe_deletion": "inframe",
    "missense_variant": "missense",
    "missense": "missense",
    "stop_gained": "nonsense",
    "synonymous_variant": "synonymous",
    "splice_acceptor_variant": "splice_lof",
    "splice_donor_variant": "splice_lof",
    "splice_acceptor": "splice_lof",
    "splice_donor": "splice_lof",
    "splice_region_variant": "splice_region",
    "splice_region": "splice_region",
    "conserved_exon_terminus_variant": "splice_lof",
    "start_lost": "missense",
    "stop_lost": "missense",
    "stop_retained": "synonymous",
    "synonymous": "synonymous",
    "nonsense": "nonsense",
    "splice_lof": "splice_lof",
    "inframe": "inframe",
}


def filter_on_consequences(df: pd.DataFrame):
    """
    Filter all variants with a consequence not found in CONSEQUENCES_MAPPING

    Args:
        df (pd.DataFrame): variant table (rates or dnm) having a consequence column
    """

    logger = logging.getLogger("logger")

    # TODO : Replace with better consequence extraction
    # When bcftools report several consequences separated by the character "&", we extract the first one
    df.consequence = [csq.split("&")[0] if isinstance(csq, str) else csq for csq in list(df.consequence)]

    filt = df.consequence.isin(CONSEQUENCES_MAPPING.keys())
    kept_df = df.loc[filt]

    logger.info(f"Before consequence filtering : {df.shape[0]} records")

    discarded_df = df.loc[~filt]
    if not discarded_df.empty:
        count_discarded = discarded_df["consequence"].value_counts()
        logger.warning(f"{discarded_df.shape[0]}/{df.shape[0]} records were discarded")
        for consequence, count in count_discarded.items():
            logger.warning(f"{consequence} : {count}")
    else:
        logger.info("All records have an acceptable consequence.")

    logger.info(f"After consequence filtering : {kept_df.shape[0]} records")

    return kept_df


def assign_meta_consequences(df: pd.DataFrame):
    """
    Assign a higher level consequence to each variant

    Args:
        df (pd.DataFrame): variant table (rates or dnm) having a consequence column
    """

    logger = logging.getLogger("logger")

    df.consequence = df.consequence.replace(CONSEQUENCES_MAPPING)

    consequence_counts = dict(df.consequence.value_counts())
    for consequence, count in consequence_counts.items():
        logger.info(f"{consequence} : {count}")

    return df


def init_log():
    """Initialize logging configuration"""

    MY_LOGGING_CONFIG = {
        "version": 1,
        "disable_existing_loggers": False,
        "formatters": {
            "default_formatter": {"format": "[%(levelname)s] %(message)s"},
        },
        "handlers": {
            "stream_handler": {
                "class": "logging.StreamHandler",
                "formatter": "default_formatter",
            },
        },
        "loggers": {
            "logger": {
                "handlers": ["stream_handler"],
                "level": "INFO",
                "propagate": True,
            }
        },
    }

    logging.config.dictConfig(MY_LOGGING_CONFIG)


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
            midpoint = cadd_lb + 3.75

    return midpoint


def plot_enrichment(df, mode="obs_exp", loess=False, df_loess=pd.DataFrame(), outdir="enrichment_plots"):
    """
    Plot enrichment based on observed/expected ratio or positive predictive value (ppv) for all bins


    Args:
        df (pd.DataFrame): contains obs/exp ratio or ppv for each bin
        mode (str, optional): observed/expected ratio or ppv. Defaults to "obs_exp".
        loess (bool, optional): plot loess results if True. Defaults to False.
        df_loess (pd.DataFrame, optional): contains loess bins obs/exp and ppv . Defaults empty.
        outdir (str, optional): output directory. Defaults to "enrichment_plots".
    """

    os.makedirs(outdir, exist_ok=True)

    # Compute confidence intervals around the observed number of DNM for each bin
    if mode == "obs_exp":
        df["ci_lower"] = df["obs"].apply(lambda x: poisson.interval(0.95, x)[0]) / df["exp"]
        df["ci_upper"] = df["obs"].apply(lambda x: poisson.interval(0.95, x)[1]) / df["exp"]

    plot_enrichment_missense(df, mode, loess, df_loess, outdir)
    plot_enrichment_nonsense(df, mode, loess, df_loess, outdir)

    # TODO : right now the inframe, and frameshift categories are added in fit_loess.py, which means that we can generate the plots only after this stage for
    # those categories, this could be done beforehand
    if not df_loess.empty:
        plot_enrichment_other(df_loess, mode, outdir)


def plot_enrichment_missense(df, mode, loess, df_loess, outdir):
    """
    Plot enrichment based on observed/expected ratio or positive predictive value (ppv) for the missense bins

    Args:
        df (pd.DataFrame): contains obs/exp ratio or ppv for each bin
        mode (str): observed/expected ratio or ppv.
        loess (bool): plot loess results if True.
        df_loess (pd.DataFrame): contains loess bins obs/exp and ppv, empty if mode obs_exp.
        outdir (str): output directory.
    """

    df = df.loc[(df.consequence == "missense") & ~(df.score.isna())].copy()

    # Technical adjustments for plotting
    plt.figure(figsize=(12, 9))
    plt.xlabel("CADD score")
    plt.ylabel(f"{mode} missense")
    df["shet"] = df["shethigh"].map({False: "low", True: "high"})
    df["region"] = df["constrained"].map({False: "unconstrained", True: "constrained"})
    df["color"] = df.shet.map({"high": "blue", "low": "orange"})

    # In observed/expected ratio mode we want to show on the plot the ratio values for the CADD range bins
    if mode == "obs_exp":
        df["x_label"] = df["score"].map({item: index for index, item in enumerate(list(df.score.unique()))})
        df["midpoint"] = df["score"].apply(lambda x: get_CADD_midpoint(x, "missense"))
        plt.axhline(y=1, color="grey", linestyle=":")  # Obs = exp line

    # Show loess values as the line plot
    if loess:
        df_loess = df_loess.loc[(df_loess.consequence == "missense") & ~(df_loess.score.isna())].copy()
        df_loess["shet"] = df_loess["shethigh"].map({False: "low", True: "high"})
        df_loess["region"] = df_loess["constrained"].map({False: "unconstrained", True: "constrained"})
        df_loess["color"] = df_loess.shet.map({"high": "blue", "low": "orange"})

    # Show loess values as the line plot
    if loess:
        sns.lineplot(df_loess, x="score", y=mode, hue="shet", color=["blue", "orange"], style="region", markers=False)
    # Show CADD range bin values as the line plot
    else:
        sns.lineplot(df, x="midpoint", y=mode, hue="shet", color=["blue", "orange"], style="region", markers=False)

    # Plot error bars for CADD range bins (not when plotting ppv)
    if mode == "obs_exp":
        for cat, min_df in df.groupby("shet"):
            plt.errorbar(
                min_df.midpoint,
                min_df.obs_exp,
                yerr=(min_df.obs_exp - min_df.ci_lower, min_df.ci_upper - min_df.obs_exp),
                fmt="o",
                color=list(min_df["color"])[0],
            )

    # Export plot
    if loess:
        plt.savefig(f"{outdir}/enrichment_missense_{mode}_loess.png", format="png")
    else:
        plt.savefig(f"{outdir}/enrichment_missense_{mode}.png", format="png")


def plot_enrichment_nonsense(df, mode, loess, df_loess, outdir):
    """
    Plot enrichment based on observed/expected ratio or positive predictive value (ppv) for the nonsense bins

    Args:
        df (pd.DataFrame): contains obs/exp ratio or ppv for each bin
        mode (str): observed/expected ratio or ppv.
        loess (bool): plot loess results if True.
        df_loess (pd.DataFrame): contains loess bins obs/exp and ppv, empty if mode obs_exp.
        outdir (str): output directory.
    """

    df = df.loc[(df.consequence == "nonsense")].copy()

    # Get midpoint that will serve as the X axis
    if mode == "obs_exp":
        df["midpoint"] = df["score"].apply(lambda x: get_CADD_midpoint(x, "nonsense"))

    # Technical adjustments for plotting
    showMarkers = True if mode == "obs_exp" else False
    plt.figure(figsize=(12, 9))
    df["shet"] = df["shethigh"].map({False: "low", True: "high"})
    plt.xlabel("CADD score")
    plt.ylabel(f"{mode} nonsense")
    plt.legend(loc="upper left", title="sHet")
    if mode == "obs_exp":
        plt.axhline(y=1, color="grey", linestyle=":")  # Obs = exp line
        df["color"] = df.shet.map({"high": "blue", "low": "orange"})

    # Show loess values as the line plot
    if loess:
        df_loess = df_loess.loc[(df_loess.consequence == "nonsense")].copy()
        df_loess["shet"] = df_loess["shethigh"].map({False: "low", True: "high"})

    # Show loess values as the line plot
    if loess:
        sns.lineplot(df_loess, x="score", y=mode, hue="shet", markers=False)
    # Show CADD range bin values as the line plot
    else:
        sns.lineplot(df, x="midpoint", y=mode, hue="shet", markers=False)

    # Plot error bars for CADD range bins (not when plotting ppv)
    if mode == "obs_exp":
        for cat, min_df in df.groupby("shet"):
            plt.errorbar(
                min_df.midpoint,
                min_df.obs_exp,
                yerr=(min_df.obs_exp - min_df.ci_lower, min_df.ci_upper - min_df.obs_exp),
                fmt="D",
                color=list(min_df["color"])[0],
            )

    # Export plot
    if loess:
        plt.savefig(f"{outdir}/enrichment_nonsense_{mode}_loess.png", format="png")
    else:
        plt.savefig(f"{outdir}/enrichment_nonsense_{mode}.png", format="png")


def plot_enrichment_other(df, mode, outdir):
    """
    Plot enrichment based on observed/expected ratio or positive predictive value (ppv) for the other (not nonsense, not missense) bins

    Args:
        df (pd.DataFrame): contains obs/exp ratio or ppv for each bin
        mode (str): observed/expected ratio or ppv.
        outdir (str): output directory.
    """

    # Keep all categories but nonsense or missense
    df = df.loc[(~df["consequence"].str.contains("nonsense|missense"))].copy()

    if mode == "obs_exp":
        df["ci_lower"] = df["obs"].apply(lambda x: poisson.interval(0.95, x)[0]) / df["exp"]
        df["ci_upper"] = df["obs"].apply(lambda x: poisson.interval(0.95, x)[1]) / df["exp"]

    # Create main plot
    plt.figure(figsize=(12, 9))
    ax = sns.pointplot(data=df, x="consequence", y=mode, hue="shethigh", dodge=True)

    # In obs_exp mode we add the confidence interval on the plot
    if mode == "obs_exp":
        # Find the x,y coordinates for each point
        x_coords = []
        y_coords = []
        for point_pair in ax.collections:
            for x, y in point_pair.get_offsets():
                x_coords.append(x)
                y_coords.append(y)

        # Order coordinates according to data frame
        y_coords_first = y_coords[: int(len(y_coords) / 2)]
        y_coords_last = y_coords[int(len(y_coords) / 2) :]
        x_coords_first = x_coords[: int(len(x_coords) / 2)]
        x_coords_last = x_coords[int(len(x_coords) / 2) :]
        y_new_coords = list()
        x_new_coords = list()
        for i, x in enumerate(x_coords_first):
            y_new_coords.append(y_coords_first[i])
            y_new_coords.append(y_coords_last[i])
            x_new_coords.append(x_coords_first[i])
            x_new_coords.append(x_coords_last[i])

        # Add error bars
        colors = ["steelblue", "coral"] * int(len(x_coords) / 2)
        ax.errorbar(
            x_new_coords,
            y_new_coords,
            yerr=(df["obs_exp"] - df["ci_lower"], df["ci_upper"] - df["obs_exp"]),
            ecolor=colors,
            fmt=" ",
            zorder=-1,
        )

        # Add obs=exp line
        plt.axhline(y=1, color="grey", linestyle=":")

    # Export figure
    plt.savefig(f"{outdir}/enrichment_others_{mode}.png", format="png")
