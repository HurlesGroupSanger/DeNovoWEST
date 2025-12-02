import pandas as pd
import gffutils
import click
import logging
import logging.config
import tempfile


def load_dnm(dnm_file):
    """
    Load DNM file

    Args:
        dnm_file (str): Path to the de novo mutation (DNM) file (TSV format).

    """

    logger = logging.getLogger("logger")

    dnm_df = pd.read_csv(dnm_file, sep="\t")

    # Create a unique variant identifier if not present
    if "varid" not in dnm_df.columns:
        dnm_df["varid"] = [f"var_{x}" for x in range(1, dnm_df.shape[0] + 1)]
        logger.info("Created unique variant identifiers 'varid' for each DNM")

    # Rename gene_id column if already present to avoid conflicts later
    if "gene_id" in dnm_df.columns:
        logger.warning("'gene_id' column already present in DNM file. Renaming it to 'original_gene_id'.")
        dnm_df.rename(columns={"gene_id": "original_gene_id"}, inplace=True)

    return dnm_df


def load_gff(gff):
    """
    Create a gffutils database from a GFF file or
    load an existing gffutils database

    Args:
        gff (str): gff file
    """

    logger = logging.getLogger("logger")

    if gff.endswith("db"):
        gff_db = gffutils.FeatureDB(gff)
    else:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".db") as tmp:
            temp_path = tmp.name
        gff_db = gffutils.create_db(gff, temp_path, merge_strategy="create_unique")
        logger.info(f"gffutils database created here : {temp_path}")

    return gff_db


def init_log(show_time=True):
    """
    Initialise logger

    Args:
        show_time (bool, optional): print time in log entries. Defaults to True.
    """

    if show_time:
        log_format = "[%(levelname)s:%(asctime)s] %(message)s"
    else:
        log_format = "%(levelname)s | %(message)s"

    MY_LOGGING_CONFIG = {
        "version": 1,
        "disable_existing_loggers": False,
        "formatters": {
            "default_formatter": {"format": log_format},
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


def annotate_with_genes(dnm_df, gff_db, gene_offset):
    """
    Annotate DNM table with gene identifiers if the DNM falls within a gene region

    Args:
        dnm_df (pd.DataFrame): DNM table
        gff_db (gffutils.FeatureDB): gffutils database
        gene_offset (int): Offset (in bp) to consider DNMs close to gene regions.
    """

    # Build a dictionary of genes per chromosome
    logger = logging.getLogger("logger")
    genes_per_chrom = build_genes_per_chrom(gff_db)
    logger.info("Built genes per chromosome dictionnary from GFF")

    # Ensure chromosome naming consistency between DNM and GFF
    dnm_df = format_chr(dnm_df, genes_per_chrom)
    chromosomes = dnm_df["chrom"].unique()

    list_annotated_dnm_df = []
    for chrom in chromosomes:

        # Get DNMs on the current chromosome
        chrom_dnm_df = dnm_df.loc[dnm_df["chrom"] == chrom].copy()

        # If the chromosome is not found in GFF, skip it but keep the DNMs
        if chrom not in genes_per_chrom:
            list_annotated_dnm_df.append(chrom_dnm_df)
            logger.warning(f"Chromosome {chrom} not found in GFF")
            continue

        # Get genes on the current chromosome
        chrom_genes = genes_per_chrom[chrom]

        # Loop over genes and annotate DNMs falling within gene regions
        for gene_id, gene_start, gene_end in chrom_genes:
            is_in_gene = (chrom_dnm_df["pos"] >= gene_start - gene_offset) & (
                chrom_dnm_df["pos"] <= gene_end + gene_offset
            )
            gene_dnm_df = chrom_dnm_df.loc[is_in_gene].copy()

            if gene_dnm_df.empty:
                continue

            gene_dnm_df["gene_id"] = gene_id
            list_annotated_dnm_df.append(gene_dnm_df)

    annotated_dnm_df = pd.concat(list_annotated_dnm_df)
    logger.info("Annotated DNM with gene identifiers")

    # Add DNMs that are outside genes regions back to the annotated DNM dataframe
    dnm_outside_genes_df = dnm_df.loc[~dnm_df["varid"].isin(annotated_dnm_df["varid"])]
    if not dnm_outside_genes_df.empty:
        logger.info(
            f"{dnm_outside_genes_df.shape[0]} DNMs found outside gene regions. They will have empty 'gene_id' field. If you have large liminal insertions/deletions, consider increasing the GENE_OFFSET value."
        )
        annotated_dnm_df = pd.concat([annotated_dnm_df, dnm_outside_genes_df])

    # Report DNMs mapped to multiple genes
    logger.info(
        f"There are {annotated_dnm_df.loc[annotated_dnm_df.varid.duplicated(), 'varid'].nunique()} DNMs that were mapped to multiple genes."
    )
    return annotated_dnm_df


def build_genes_per_chrom(gff_db):
    """
    Build a dictionary of genes per chromosome from GFF

    Args:
        gff_db (gffutils.FeatureDB): gffutils database
    """

    genes_per_chrom = dict()
    for gene in gff_db.features_of_type("gene", order_by="start"):
        gene_id = gene.id
        gene_chrom = gene.chrom
        gene_start = gene.start
        gene_end = gene.end

        if gene_chrom not in genes_per_chrom:
            genes_per_chrom[gene_chrom] = []
        genes_per_chrom[gene_chrom].append(
            (
                gene_id,
                gene_start,
                gene_end,
            )
        )

    return genes_per_chrom


def format_chr(dnm_df, genes_per_chrom):
    """
    Ensure chromosome naming consistency between DNM and GFF

    Args:
        dnm_df (pd.DataFrame): DNM table
        genes_per_chrom (dict): genes per chromosome from GFF
    """

    # Determine if 'chr' prefix is used in GFF
    chr_prefix_in_gff = False
    for chrom in genes_per_chrom:
        if chrom.startswith("chr"):
            chr_prefix_in_gff = True
        break

    # Determine if 'chr' prefix is used in DNM
    chr_prefix_in_dnm = False
    if dnm_df.iloc[0].chrom.startswith("chr"):
        chr_prefix_in_dnm = True

    # Add 'chr' prefix to DNM chromosome column if prefix found in GFF but not in DNM
    if chr_prefix_in_gff and not chr_prefix_in_dnm:
        dnm_df["chrom"] = "chr" + dnm_df["chrom"].astype(str)
    # Remove 'chr' prefix from DNM chromosome column if prefix found in DNM but not in GFF
    elif not chr_prefix_in_gff and chr_prefix_in_dnm:
        dnm_df["chrom"] = dnm_df["chrom"].astype(str).str.replace("chr", "", regex=False)

    return dnm_df


def export_dnm_with_genes(dnm_annotated_with_gene_df, output_file):
    """
    Export DNM annotated with gene identifiers to a TSV file

    Args:
        dnm_annotated_with_gene_df (pd.DataFrame): DNM table annotated with gene identifiers
        output_file (str): Path to the output TSV file
    """

    # Bring gene_id column to the front
    dnm_annotated_with_gene_df = dnm_annotated_with_gene_df[
        ["gene_id"] + [c for c in dnm_annotated_with_gene_df.columns if c != "gene_id"]
    ]

    # Sort DNMs by chromosome and position
    dnm_annotated_with_gene_df.sort_values(by=["chrom", "pos"], inplace=True)

    # Export to TSV
    dnm_annotated_with_gene_df.to_csv(output_file, sep="\t", index=False)


@click.command()
@click.argument("dnm")
@click.argument("gff")
@click.argument("output")
@click.option("--gene-offset", default=50, help="Offset (in bp) to consider DNMs close to gene regions.")
def map_dnm_to_genes(dnm, gff, output, gene_offset):
    """
    Map de novo mutations (DNMs) to genes using a GFF annotation file.

    Args:
        dnm (str): Path to the de novo mutation (DNM) file (TSV format).
        gff (str): Path to the GFF3 annotation file or gffutils database.
        output (str): Path to the gene-annotated DNM output file (TSV format).
    """

    # Load the DNM data
    dnm_df = load_dnm(dnm)

    # Load the GFF database
    gff_db = load_gff(gff)

    # Annotate DNM with gene identifiers
    dnm_annotated_with_gene_df = annotate_with_genes(dnm_df, gff_db, gene_offset)

    # Export the results
    export_dnm_with_genes(dnm_annotated_with_gene_df, output)


if __name__ == "__main__":
    init_log()
    map_dnm_to_genes()
