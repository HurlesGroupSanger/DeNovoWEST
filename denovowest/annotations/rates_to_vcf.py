#!/usr/bin/env python
import logging
import sys

import click
import gffutils
import pandas as pd

from denovowest.utils.io_helpers import extract_ensembl_gene_id_without_version
from denovowest.utils.log import init_log
from denovowest.utils.params import CONSEQUENCES_SEVERITIES

# TODO: better handle versioned and non-versioned gene ids

#################
# Constants     #
#################

# Consequences that classify a variant as an indel
INDELS = frozenset(["inframe", "frameshift", "inframe_insertion", "inframe_deletion"])


#################
# Helpers       #
#################


def _get_info_field(info: str, key: str) -> str | None:
    """Extract a single value from a VCF INFO string by field name.

    Args:
        info: Raw VCF INFO string (e.g. ``"GENE=ENSG00000001;BCSQ=missense|..."``).
        key: INFO field name to look up (e.g. ``"GENE"`` or ``"BCSQ"``).

    Returns:
        The field value, or ``None`` if the key is absent.
    """

    prefix = key + "="
    for part in info.split(";"):
        if part.startswith(prefix):
            return part[len(prefix):]
    return None


@click.group()
def cli():
    pass


@cli.command()
@click.argument("rates")
@click.argument("fasta")
@click.argument("out_vcf")
def rates_to_vcf(rates: str, fasta: str, out_vcf: str) -> None:
    """Transform a rates file (tabular format) into VCF format for use with bcftools csq.

    Args:
        rates: Rate file.
        fasta: FASTA file or FASTA index file used to compute the rate file.
        out_vcf: Destination of the output VCF file.
    """

    # Load rates file early so we can inspect columns before writing the VCF header
    rates_df = pd.read_csv(rates, sep="\t")
    has_gene_symbol = "gene_symbol" in rates_df.columns

    # Retrieve contig length from FASTA file
    if fasta.endswith(".fai"):
        vcf_header_df = pd.read_csv(f"{fasta}", sep="\t", header=None)
    else:
        vcf_header_df = pd.read_csv(f"{fasta}.fai", sep="\t", header=None)

    # Write VCF header (SYMBOL INFO field declared only when present in the rates file)
    with open(out_vcf, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        for idx, row in vcf_header_df.iterrows():
            chrom = row[0]
            length = row[1]
            f.write(f"##contig=<ID={chrom},length={length}>\n")
        f.write('##INFO=<ID=GENE,Number=.,Type=String,Description="Gene identifier">\n')
        if has_gene_symbol:
            f.write('##INFO=<ID=SYMBOL,Number=.,Type=String,Description="Gene symbol">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    # Filter duplicate variants
    duplicates = rates_df[["gene_id", "chrom", "pos", "ref", "alt"]].duplicated()
    rates_df = rates_df[~duplicates]

    # We need to sort the rates file otherwise bcftoolscsq fails
    rates_df = rates_df.sort_values(["chrom", "pos"])

    # Write VCF (including gene id and gene symbol in INFO field to distinguish between overlapping genes)
    with open(out_vcf, "a") as f:
        for idx, row in rates_df.iterrows():
            info = f"GENE={row.gene_id}"
            if has_gene_symbol and pd.notna(row.gene_symbol):
                info += f";SYMBOL={row.gene_symbol}"
            f.write(f"{row.chrom}\t{row.pos}\t.\t{row.ref}\t{row.alt}\t.\t.\t{info}\n")


def extract_worst_consequence(
    vcf_df: pd.DataFrame,
    gff_db: gffutils.FeatureDB,
    force_indel_annotation: bool,
) -> pd.DataFrame:
    """Extract the worst bcftools csq consequence for each variant.

    BCFtools csq returns a consequence string per transcript. We extract the one
    that contains the worst consequence. Sometimes there will be multiple
    annotations for a single transcript (e.g. synonymous&splice_donor); in that
    case we retrieve the worst annotation among the multiple annotations.

    Args:
        vcf_df: Variant file annotated with bcftools csq.
        gff_db: gffutils database.
        force_indel_annotation: Force indels to be annotated as indels (i.e. inframe, frameshift).

    Returns:
        The VCF df with two new columns: consequence and full_consequence_bcftoolscsq.
    """

    both_use_version = do_both_sources_use_version(gff_db, vcf_df)
    if not both_use_version:
        ensembl_gene_id_map_version = extract_ensembl_gene_id_without_version(gff_db)

    list_worst_csq = []
    list_csqs = []
    for idx, record in vcf_df.iterrows():

        # Get gene name(s) — use key-based lookup so extra INFO fields don't shift positions
        gene_id = _get_info_field(record.INFO, "GENE")

        # TODO : Handle the case where this field is missing or has a different name, find a clever way to support all gff
        try:
            if not both_use_version:
                gene_names = gff_db[ensembl_gene_id_map_version[gene_id]].attributes["gene_name"]
            else:
                gene_names = gff_db[gene_id].attributes["gene_name"]
        except KeyError:
            if not both_use_version:
                gene_names = gff_db[ensembl_gene_id_map_version[gene_id]].attributes["Name"]
            else:
                gene_names = gff_db[gene_id].attributes["Name"]

        # Look at all the consequences returned by bcftoolscsq
        try:
            bcsq_value = _get_info_field(record.INFO, "BCSQ")
            if bcsq_value is None:
                raise IndexError
            list_csq = bcsq_value.split(",")

            # If there are multiple consequences, we take the worst one
            worst_csq_value = 1
            worst_csq = ""
            for cur_csq in list_csq:

                # Check that the annotation is on our gene of interest (overlapping genes)
                if cur_csq.split("|")[1] not in gene_names:
                    continue

                # Extract the worst consequence
                csqs = cur_csq.split("|")[0]
                for csq in csqs.split("&"):
                    if CONSEQUENCES_SEVERITIES[csq] > worst_csq_value:
                        worst_csq_value = CONSEQUENCES_SEVERITIES[csq]
                        worst_csq = csq
                        worst_csqs = csqs

            if worst_csq:

                csq = worst_csq
                csqs = worst_csqs

                # If we want all coding indels to be annotated as inframe or frameshift
                if force_indel_annotation and is_indel(record):
                    csq, csqs = assign_indel_csq(csq, csqs)

                # We keep the long bcftools string uniquely if it differs from the
                # single worst consequence
                if csqs == csq:
                    csqs = ""

            else:
                csq = ""
                csqs = ""

        except IndexError:
            # bcftoolscsq returns sometimes empty consequences
            csq = ""
            csqs = ""

        list_worst_csq.append(csq)
        list_csqs.append(csqs)

    vcf_df["consequence"] = list_worst_csq
    vcf_df["full_consequence_bcftoolscsq"] = list_csqs

    return vcf_df


def do_both_sources_use_version(gff_db: gffutils.FeatureDB, vcf_df: pd.DataFrame) -> bool:
    """Check whether both the GFF and the VCF use versioned ENSEMBL gene IDs.

    Args:
        gff_db: gffutils database.
        vcf_df: Variant file with a GENE INFO field.

    Returns:
        True if both sources include a version suffix (e.g. ENSG0000001.2).
    """

    gff_gene_id = next(gff_db.features_of_type("gene", order_by="start")).id

    vcf_gene_id = ""
    for _, record in vcf_df.iterrows():
        vcf_gene_id = _get_info_field(record.INFO, "GENE") or ""
        if vcf_gene_id:
            break

    return "." in gff_gene_id and "." in vcf_gene_id


def is_indel(record: pd.Series) -> bool:
    """Test whether the variant is an indel.

    Args:
        record: Variant as a pandas Series with REF and ALT fields.

    Returns:
        True if the variant is an indel, False if it is a SNP.
    """

    return len(record.REF) != 1 or len(record.ALT) != 1


def assign_indel_csq(csq: str, csqs: str) -> tuple[str, str]:
    """Coerce a consequence to inframe or frameshift for coding indels.

    Some indels are not called as inframe or frameshift (e.g. splice_donor).
    This means that no SNP-based score (CADD, dbNSFP...) can be assigned for
    them. Here we offer an option to force indel indels to be called as inframe
    or frameshift even if bcftoolscsq returns a different consequence.

    Args:
        csq: The worst consequence found for this variant (e.g. splice_donor).
        csqs: The full bcftoolscsq string (e.g. synonymous&splice_donor).

    Returns:
        Tuple of (csq, csqs), potentially coerced to an indel consequence.
    """

    # If the worst consequence is already inframe or frameshift, nothing to do
    if csq in INDELS:
        return csq, csqs

    # If inframe is in the list of consequences, we set the consequence as inframe
    if "inframe" in csqs:
        csq = [x for x in csqs.split("&") if x.startswith("inframe")][0]
        return csq, csqs

    # If frameshift is in the list of consequences, we set the consequence as frameshift
    if "frameshift" in csqs:
        return "frameshift", csqs

    # If the worst consequence is non coding we do nothing (non coding regions are not handled in DNW)
    if CONSEQUENCES_SEVERITIES[csq] < CONSEQUENCES_SEVERITIES["coding_sequence"]:
        return csq, csqs

    # If the worst consequence is inframe-like, we set the consequence as inframe
    if CONSEQUENCES_SEVERITIES[csq] <= CONSEQUENCES_SEVERITIES["inframe"]:
        return "inframe", csqs

    # Otherwise we consider it as frameshift
    return "frameshift", csqs


@cli.command()
@click.argument("vcf")
@click.argument("rates")
@click.argument("gff_db_path")
@click.argument("out_rates")
@click.option(
    "--force-indel-annotation",
    is_flag=True,
    default=False,
)
def vcf_to_rates(vcf: str, rates: str, gff_db_path: str, out_rates: str, force_indel_annotation: bool) -> None:
    """Merge the VCF file annotated with bcftoolscsq with a rate file (tabular).

    Args:
        vcf: Path of the VCF file.
        rates: Path to the rates file.
        gff_db_path: Path to gffutils database.
        out_rates: Path of the output rates file.
        force_indel_annotation: Force indels to be annotated as indels (i.e. inframe, frameshift).
    """

    init_log()
    logger = logging.getLogger("logger")

    # Read rates
    rates_df = pd.read_csv(rates, sep="\t", dtype={"chrom": str})

    # Read VCF
    vcf_columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    try:
        if vcf.endswith("gz"):
            vcf_df = pd.read_csv(vcf, sep="\t", comment="#", header=None, compression="gzip")
        else:
            vcf_df = pd.read_csv(vcf, sep="\t", comment="#", header=None)
        vcf_df.columns = vcf_columns
        vcf_df["#CHROM"] = vcf_df["#CHROM"].astype(str)

    except pd.errors.EmptyDataError:
        if rates_df.empty:
            logger.warning("Rates file is empty")
            rates_df["consequence"] = pd.NA
            rates_df["full_consequence_bcftoolscsq"] = pd.NA
            rates_df.to_csv(out_rates, sep="\t", index=False)
            sys.exit(0)
        else:
            logger.error("Bcftools csq file is empty")
            sys.exit(1)

    # Load gffutils database
    gff_db = gffutils.FeatureDB(gff_db_path)

    # Extract BCFtools consequence in a separate column
    vcf_df = extract_worst_consequence(vcf_df, gff_db, force_indel_annotation)

    # Extract gene id to distinguish between overlapping genes
    vcf_df["gene_id"] = [_get_info_field(x, "GENE") for x in vcf_df.INFO]

    # Extract gene symbol from INFO if encoded (written by rates_to_vcf when the rates file had gene_symbol)
    vcf_df["gene_symbol"] = vcf_df.INFO.str.extract(r"SYMBOL=([^;]+)")

    # If the rates file already has a 'consequence' column (e.g. from a previous annotation
    # run), rename it to preserve it — the new bcftools csq annotation takes the canonical name
    # because 'consequence' is a fixed column name expected downstream by the simulation.
    if "consequence" in rates_df.columns:
        logger.warning(
            "Column 'consequence' already exists in the rates file — the existing column "
            "will be renamed to 'consequence_prev'. The new bcftools csq annotation "
            "will be written as 'consequence'."
        )
        rates_df = rates_df.rename(columns={"consequence": "consequence_prev"})

    # Build right-side column list: include gene_symbol from VCF only when rates_df
    # doesn't already carry it (legacy rates files without gene_symbol column)
    right_cols = ["#CHROM", "POS", "REF", "ALT", "consequence", "full_consequence_bcftoolscsq", "gene_id"]
    if "gene_symbol" not in rates_df.columns and vcf_df["gene_symbol"].notna().any():
        right_cols.append("gene_symbol")

    # Merge files
    out_rates_df = rates_df.merge(
        vcf_df[right_cols],
        left_on=["chrom", "pos", "ref", "alt", "gene_id"],
        right_on=["#CHROM", "POS", "REF", "ALT", "gene_id"],
    )
    out_rates_df.drop(["#CHROM", "POS", "REF", "ALT"], axis=1, inplace=True)

    assert out_rates_df.shape[0] == rates_df.shape[0]

    # Export
    out_rates_df.to_csv(out_rates, sep="\t", index=False)


if __name__ == "__main__":
    cli()
