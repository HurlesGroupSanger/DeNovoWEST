import pandas as pd
import numpy as np
import click
import logging
import logging.config
from jinja2 import Environment, FileSystemLoader
from jinja2.exceptions import TemplateNotFound
import os
import pysam
import subprocess


def compare_enrichment_results(enrichment_results_run1_df, enrichment_results_run2_df, bonferonni_threshold):

    # Check how many genes were analysed in both runs
    shared_genes = get_shared_genes(enrichment_results_run1_df, enrichment_results_run2_df)

    # Check how many genes were found significant in both runs
    significant_genes_run1 = get_significant_genes(enrichment_results_run1_df, bonferonni_threshold, "run 1")
    significant_genes_run2 = get_significant_genes(enrichment_results_run2_df, bonferonni_threshold, "run 2")
    shared_significant_genes = get_shared_significant_genes(significant_genes_run1, significant_genes_run2)

    # Find the distribution of observed and expected ratios between the two runs
    ratios = find_observed_and_expected_ratios_between_runs(
        enrichment_results_run1_df.loc[enrichment_results_run1_df.symbol.isin(shared_genes)].sort_values("symbol"),
        enrichment_results_run2_df.loc[enrichment_results_run2_df.symbol.isin(shared_genes)].sort_values("symbol"),
    )

    return shared_genes, shared_significant_genes, significant_genes_run1, significant_genes_run2, ratios


def get_shared_genes(enrichment_results_run1_df, enrichment_results_run2_df):

    logger = logging.getLogger("logger")

    genes_run1 = set(enrichment_results_run1_df.symbol)
    genes_run2 = set(enrichment_results_run2_df.symbol)

    shared_genes = genes_run1 & genes_run2

    logger.info(f"{len(shared_genes)} shared genes, {len(genes_run1)} in run 1, {len(genes_run2)} in run 2")
    logger.info(f"{len(genes_run1 - shared_genes)} genes found only in run 1 : {genes_run1 - shared_genes}")
    logger.info(f"{len(genes_run2 - shared_genes)} genes found only in run 2 : {genes_run2 - shared_genes}")
    logger.info("---------")

    return shared_genes


def get_significant_genes(enrichment_results_df, bonferonni_threshold, run_name):

    logger = logging.getLogger("logger")

    significant_genes = set(
        enrichment_results_df.loc[enrichment_results_df["p-value"] < bonferonni_threshold, "symbol"]
    )

    logger.info(f"{len(significant_genes)} significant genes in {run_name}")

    return significant_genes


def get_shared_significant_genes(significant_genes_run1, significant_genes_run2):

    logger = logging.getLogger("logger")

    shared_significant_genes = significant_genes_run1 & significant_genes_run2
    logger.info(
        f"{len(significant_genes_run1 - shared_significant_genes)} significant genes found only in run 1 : {significant_genes_run1 - shared_significant_genes}"
    )
    logger.info(
        f"{len(significant_genes_run2 - shared_significant_genes)} significant genes found only in run 2 : {significant_genes_run2 - shared_significant_genes}"
    )
    logger.info("---------")

    return shared_significant_genes


def find_observed_and_expected_ratios_between_runs(enrichment_results_run1_df, enrichment_results_run2_df):

    logger = logging.getLogger("logger")
    ratios = dict()

    observed_ratios = enrichment_results_run1_df["observed"] / enrichment_results_run2_df["observed"]
    expected_ratios = enrichment_results_run1_df["expected"] / enrichment_results_run2_df["expected"]

    # Q1
    q1_observed_ratio = np.percentile(observed_ratios.replace([np.inf, -np.inf], np.nan).dropna(), 25)
    q1_expected_ratio = np.percentile(expected_ratios.replace([np.inf, -np.inf], np.nan).dropna(), 25)
    ratios["q1_observed"] = q1_observed_ratio
    ratios["q1_expected"] = q1_expected_ratio

    logger.info(f"Q1 observed ratio between run1 and run2 : {q1_observed_ratio}")
    logger.info(f"Q1 expected ratio between run1 and run2 : {q1_expected_ratio}")

    # Median
    median_observed_ratio = np.median(observed_ratios.replace([np.inf, -np.inf], np.nan).dropna())
    median_expected_ratio = np.median(expected_ratios.replace([np.inf, -np.inf], np.nan).dropna())
    ratios["median_observed"] = median_observed_ratio
    ratios["median_expected"] = median_expected_ratio

    logger.info(f"Median observed ratio between run1 and run2 : {median_observed_ratio}")
    logger.info(f"Median expected ratio between run1 and run2 : {median_expected_ratio}")

    # Q3
    q3_observed_ratio = np.percentile(observed_ratios.replace([np.inf, -np.inf], np.nan).dropna(), 75)
    q3_expected_ratio = np.percentile(expected_ratios.replace([np.inf, -np.inf], np.nan).dropna(), 75)
    ratios["q3_observed"] = q3_observed_ratio
    ratios["q3_expected"] = q3_expected_ratio

    logger.info(f"Q3 observed ratio between run1 and run2 : {q3_observed_ratio}")
    logger.info(f"Q3 expected ratio between run1 and run2 : {q3_expected_ratio}")
    logger.info("---------")

    return ratios


def get_gene_enrichment_infos(gene, enrichment_results_df):

    gene_enrichment_info = dict()
    gene_enrichment_info["p_value"] = enrichment_results_df.loc[enrichment_results_df.symbol == gene, "p-value"].iloc[0]
    gene_enrichment_info["observed"] = enrichment_results_df.loc[enrichment_results_df.symbol == gene, "observed"].iloc[
        0
    ]
    gene_enrichment_info["expected"] = enrichment_results_df.loc[enrichment_results_df.symbol == gene, "expected"].iloc[
        0
    ]
    gene_enrichment_info["info"] = enrichment_results_df.loc[enrichment_results_df.symbol == gene, "info"].iloc[0]

    return gene_enrichment_info


def get_gene_infos(gene_id):

    gene_biomart_infos = get_gene_infos_from_biomart(gene_id)
    gene_shet_constraint_infos = get_gene_infos_shet_constraint(gene_id)

    return {**gene_biomart_infos, **gene_shet_constraint_infos}


def get_gene_infos_from_biomart(gene_id):

    biomart_df = pd.read_csv("/nfs/ddd0/resources/biomart/biomart_protein_coding_genes_GRCh38p14.tsv", sep="\t")
    biomart_s = biomart_df.loc[biomart_df["Gene stable ID"] == gene_id.split(".")[0]].iloc[0]

    gene_biomart_infos = dict()
    gene_biomart_infos["symbol"] = biomart_s["HGNC symbol"]
    gene_biomart_infos["hgnc_id"] = biomart_s["HGNC ID"]
    gene_biomart_infos["chrom"] = biomart_s["Chromosome/scaffold name"]
    gene_biomart_infos["start"] = biomart_s["Gene start (bp)"]
    gene_biomart_infos["end"] = biomart_s["Gene end (bp)"]

    return gene_biomart_infos


def get_gene_infos_shet_constraint(gene_id):

    shet_constraint_df = pd.read_csv(
        "/lustre/scratch123/hgi/mdt1/teams/hurles/ed11/DNW/DeNovoWEST/resources/genes_shethigh_constrained_status.tsv",
        sep="\t",
    )

    gene_shet_constraint_infos = dict()

    try:
        shet_constraint_s = shet_constraint_df.loc[shet_constraint_df["ensembl_gene_id"] == gene_id.split(".")[0]].iloc[
            0
        ]

        gene_shet_constraint_infos["shethigh"] = shet_constraint_s["shethigh"]
        gene_shet_constraint_infos["constrained"] = shet_constraint_s["constrained"]

    except IndexError:

        gene_shet_constraint_infos["shethigh"] = "Not found"
        gene_shet_constraint_infos["constrained"] = "Not found"

    return gene_shet_constraint_infos


def get_gene_rates_infos(rates_file, gene_id, gene_infos):

    tabix_file = pysam.TabixFile(rates_file)

    # Get rates file columns
    res = subprocess.run(f"zcat {rates_file} | head -n 1", shell=True, capture_output=True, text=True)
    columns = res.stdout.strip().split("\t")

    records = tabix_file.fetch("chr" + gene_infos["chrom"], gene_infos["start"], gene_infos["end"])
    gene_rates = pd.DataFrame([x.split("\t") for x in list(records)], columns=columns)
    gene_rates = gene_rates.loc[gene_rates.gene_id == gene_id]
    gene_rates["ppv"] = [float(x) if x != "" else np.nan for x in gene_rates["ppv"]]
    gene_rates["prob"] = gene_rates["prob"].astype(float)

    gene_rates["weight"] = gene_rates["ppv"] * gene_rates["prob"]

    rates_infos_df = pd.DataFrame(gene_rates.consequence.value_counts())
    rates_infos_df.index.name = "consequence"
    rates_infos_df.columns = ["nb_mutations"]

    list_sum_ppv = list()
    for consequence in rates_infos_df.index:
        list_sum_ppv.append(gene_rates.loc[gene_rates.consequence == consequence, "weight"].sum())

    rates_infos_df["sum_weight"] = list_sum_ppv
    rates_infos_df["consequence"] = rates_infos_df.index
    rates_infos_df["nb_mutations"] = rates_infos_df["nb_mutations"].astype(int)

    inframe_weight = 0.506993 if str(gene_infos["shethigh"]) == "True" else 0.159436
    inframe = [
        -1,
        gene_rates.loc[gene_rates.consequence == "missense", "prob"].sum() * 0.03 * inframe_weight,
        "inframe",
    ]

    frameshift_weight = 0.907829 if str(gene_infos["shethigh"]) == "True" else 0.530357
    frameshift = [
        -1,
        gene_rates.loc[gene_rates.consequence == "nonsense", "prob"].sum() * 1.3 * frameshift_weight,
        "frameshift",
    ]

    new_rows_df = pd.DataFrame([inframe, frameshift], columns=rates_infos_df.columns)
    rates_infos_df = pd.concat([rates_infos_df, new_rows_df], ignore_index=True)

    rates_infos_df["proportion"] = rates_infos_df["sum_weight"] * 100 / rates_infos_df["sum_weight"].sum()
    rates_infos_df["proportion"] = [round(x, 2) for x in rates_infos_df["proportion"]]

    return rates_infos_df[["consequence", "nb_mutations", "sum_weight", "proportion"]]


def create_report(
    gene,
    gene_infos,
    gene_enrichment_run1,
    gene_enrichment_run2,
    ratios,
    dnm_run1_df,
    dnm_run2_df,
    rates_info_run1_df,
    rates_info_run2_df,
    outdir,
):

    env = Environment(loader=FileSystemLoader("./"))
    template = env.get_template("misc/compare_gene_template.html")

    template_vars = dict()
    template_vars["gene_identifier"] = gene
    template_vars["p_value_run1"] = gene_enrichment_run1["p_value"]
    template_vars["observed_run1"] = gene_enrichment_run1["observed"]
    template_vars["expected_run1"] = gene_enrichment_run1["expected"]
    template_vars["p_value_run2"] = gene_enrichment_run2["p_value"]
    template_vars["observed_run2"] = gene_enrichment_run2["observed"]
    template_vars["expected_run2"] = gene_enrichment_run2["expected"]

    template_vars["ratio_q1_observed"] = ratios["q1_observed"]
    template_vars["ratio_q1_expected"] = ratios["q1_expected"]
    template_vars["ratio_q2_observed"] = ratios["median_observed"]
    template_vars["ratio_q2_expected"] = ratios["median_expected"]
    template_vars["ratio_q3_observed"] = ratios["q3_observed"]
    template_vars["ratio_q3_expected"] = ratios["q3_expected"]
    template_vars["ratio_gene_observed"] = gene_enrichment_run1["observed"] / gene_enrichment_run2["observed"]
    template_vars["ratio_gene_expected"] = gene_enrichment_run1["expected"] / gene_enrichment_run2["expected"]

    template_vars["gene_shethigh"] = gene_infos["shethigh"]
    template_vars["gene_constrained"] = gene_infos["constrained"]
    template_vars["gene_symbol"] = gene_infos["symbol"]

    dnm_run1_df["varid"] = dnm_run1_df["pos"].astype(str) + "_" + dnm_run1_df["ref"] + "_" + dnm_run1_df["alt"]
    dnm_run1_df["proportion"] = dnm_run1_df["ppv"] * 100 / dnm_run1_df["ppv"].sum()
    dnm_run1_df["proportion"] = [round(x, 2) for x in dnm_run1_df["proportion"]]
    html_table_observed_dnm_run1 = dfToHtmlTable(
        dnm_run1_df[["varid", "consequence", "score", "ppv", "proportion"]],
        "run1-observed_dnm",
    )
    template_vars["table_observed_dnm_run1"] = html_table_observed_dnm_run1

    dnm_run2_df["varid"] = dnm_run2_df["pos"].astype(str) + "_" + dnm_run2_df["ref"] + "_" + dnm_run2_df["alt"]
    dnm_run2_df["proportion"] = dnm_run2_df["ppv"] * 100 / dnm_run2_df["ppv"].sum()
    dnm_run2_df["proportion"] = [round(x, 2) for x in dnm_run2_df["proportion"]]
    html_table_observed_dnm_run2 = dfToHtmlTable(
        dnm_run2_df[["varid", "consequence", "BayesDel_addAF_score", "ppv", "proportion"]],
        "run2-observed_dnm",
    )
    template_vars["table_observed_dnm_run2"] = html_table_observed_dnm_run2

    html_table_expected_dnm_run1 = dfToHtmlTable(rates_info_run1_df, "run1-expected_dnm")
    template_vars["table_expected_dnm_run1"] = html_table_expected_dnm_run1

    html_table_expected_dnm_run2 = dfToHtmlTable(rates_info_run2_df, "run2-expected_dnm")
    template_vars["table_expected_dnm_run2"] = html_table_expected_dnm_run2

    outputText = template.render(template_vars)
    with open(f"{outdir}/{gene}.html", "w") as f:
        f.write(outputText)


def dfToHtmlTable(df, id):
    """
    Create the HTML code corresponding to a dataframe

    Args:
        df (pd.DataFrame): data

    Returns:
        str: The data frame as an HTML string
    """

    header = df.columns.values.tolist()

    html = f'<table id="{id}" class="display transcript-table">'
    html += "<thead><tr>"
    for column in header:
        html += "<th>" + column + "</th>"
    html += "</tr></thead>"

    html += "<tfoot><tr>"
    for column in header:
        html += "<th>" + column + "</th>"
    html += "</tr></tfoot>"

    html += "<tbody>"
    for index, row in df.iterrows():
        rowList = row.tolist()
        html += "<tr>"
        for value in rowList:
            html += "<td>" + str(value) + "</td>"
        html += "</tr>"

    html += "</tbody>"
    html += "</table>"

    return html


def init_log():
    """Initialize logging configuration"""

    MY_LOGGING_CONFIG = {
        "version": 1,
        "disable_existing_loggers": False,
        "formatters": {
            "default_formatter": {"format": "[%(levelname)s:%(asctime)s] %(message)s"},
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


def log_configuration(conf):
    """_summary_

    Args:
        conf (_type_): _description_
    """

    logger = logging.getLogger("logger")
    logger.info("----------")
    for key, value in sorted(conf.items()):
        logger.info(f"{key} : {value}")
    logger.info("----------")


@click.command()
@click.argument("enrichment_results_run1")
@click.argument("enrichment_results_run2")
@click.option("--dnm_table_run1")
@click.option("--dnm_table_run2")
@click.option("--rates_file_run1")
@click.option("--rates_file_run2")
@click.option(
    "--bonferonni_threshold",
    type=float,
    help="Bonferonni threshold used to determine gene significance",
    default=1.33e-6,
)
@click.option(
    "--outdir",
    type=str,
    default="out_compareruns",
)
def main(
    enrichment_results_run1,
    enrichment_results_run2,
    dnm_table_run1,
    dnm_table_run2,
    rates_file_run1,
    rates_file_run2,
    bonferonni_threshold,
    outdir,
):

    init_log()
    log_configuration(click.get_current_context().params)

    os.makedirs(outdir, exist_ok=True)

    enrichment_results_run1_df = pd.read_csv(enrichment_results_run1, sep="\t")
    enrichment_results_run2_df = pd.read_csv(enrichment_results_run2, sep="\t")

    shared_genes, shared_significant_genes, significant_genes_run1, significant_genes_run2, ratios = (
        compare_enrichment_results(enrichment_results_run1_df, enrichment_results_run2_df, bonferonni_threshold)
    )

    if dnm_table_run1:
        dnm_run1_df = pd.read_csv(dnm_table_run1, sep="\t", dtype={"pos": int})

    if dnm_table_run2:
        dnm_run2_df = pd.read_csv(dnm_table_run2, sep="\t", dtype={"pos": int})

    os.makedirs(f"{outdir}/significant_run1", exist_ok=True)
    for gene in significant_genes_run1 - shared_significant_genes:
        gene_enrichment_run1 = get_gene_enrichment_infos(gene, enrichment_results_run1_df)
        gene_enrichment_run2 = get_gene_enrichment_infos(gene, enrichment_results_run2_df)

        gene_infos = get_gene_infos(gene)

        rates_info_run1_df = get_gene_rates_infos(rates_file_run1, gene, gene_infos)
        rates_info_run2_df = get_gene_rates_infos(rates_file_run2, gene, gene_infos)

        create_report(
            gene,
            gene_infos,
            gene_enrichment_run1,
            gene_enrichment_run2,
            ratios,
            dnm_run1_df.loc[dnm_run1_df.gene_id == gene].copy(),
            dnm_run2_df.loc[dnm_run2_df.gene_id == gene].copy(),
            rates_info_run1_df,
            rates_info_run2_df,
            f"{outdir}/significant_run1",
        )

    os.makedirs(f"{outdir}/significant_run2", exist_ok=True)
    for gene in significant_genes_run2 - shared_significant_genes:
        gene_enrichment_run1 = get_gene_enrichment_infos(gene, enrichment_results_run1_df)
        gene_enrichment_run2 = get_gene_enrichment_infos(gene, enrichment_results_run2_df)

        gene_infos = get_gene_infos(gene)

        rates_info_run1_df = get_gene_rates_infos(rates_file_run1, gene, gene_infos)
        rates_info_run2_df = get_gene_rates_infos(rates_file_run2, gene, gene_infos)

        create_report(
            gene,
            gene_infos,
            gene_enrichment_run1,
            gene_enrichment_run2,
            ratios,
            dnm_run1_df.loc[dnm_run1_df.gene_id == gene].copy(),
            dnm_run2_df.loc[dnm_run2_df.gene_id == gene].copy(),
            rates_info_run1_df,
            rates_info_run2_df,
            f"{outdir}/significant_run2",
        )


if __name__ == "__main__":
    main()
