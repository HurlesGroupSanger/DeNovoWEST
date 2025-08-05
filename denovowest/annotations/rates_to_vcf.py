#!/usr/bin/env python
import click
import gffutils
import pandas as pd
import sys
import logging

from denovowest.utils.log import init_log
from denovowest.utils.params import CONSEQUENCES_SEVERITIES


@click.group()
def cli():
    pass


@cli.command()
@click.argument("rates")
@click.argument("fasta")
@click.argument("out_vcf")
def rates_to_vcf(rates, fasta, out_vcf):
    """Transforms a rate file (tabular format) in VCF format to use before calling bcftoolscsq

    Args:
        rates (file): Rate file
        fasta (file): FASTA file or fasta index file used to compute the rate file
        out_vcf (str): Destination of the output VCF file
    """

    # Retrieve contig length from FASTA file
    if fasta.endswith(".fai"):
        vcf_header_df = pd.read_csv(f"{fasta}", sep="\t", header=None)
    else:
        vcf_header_df = pd.read_csv(f"{fasta}.fai", sep="\t", header=None)

    # Write it in the VCF header format
    with open(out_vcf, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        for idx, row in vcf_header_df.iterrows():
            chrom = row[0]
            length = row[1]
            f.write(f"##contig=<ID={chrom},length={length}>\n")
        f.write('##INFO=<ID=GENE,Number=.,Type=String,Description="Gene identifier">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    # Extract information from rates file and transform it in VCF format
    rates_df = pd.read_csv(rates, sep="\t")

    # Filter duplicate variants
    duplicates = rates_df[["gene_id", "chrom", "pos", "ref", "alt"]].duplicated()
    rates_df = rates_df[~duplicates]

    # We need to sort the rates file otherwise bcftoolscsq fails
    rates_df = rates_df.sort_values(["chrom", "pos"])

    # Write VCF (including gene id in INFO field to distinguish between overlapping genes)
    with open(out_vcf, "a") as f:
        for idx, row in rates_df.iterrows():
            f.write(f"{row.chrom}\t{row.pos}\t.\t{row.ref}\t{row.alt}\t.\t.\tGENE={row.gene_id}\n")


def extract_worst_consequence(vcf_df, gff_db, force_indel_annotation):
    """
    BCFtools csq returns a consequence string per transcript.
    We extract the one that contains the worst consequence.

    Sometimes there will be multiple annotations for a single transcript (e.g. synonymous&splice_donor),
    in this case we retrieve the worst annotation among the multiple annotations

    Args:
        vcf_df (pd.DataFrame): Variant file annotated with bcftools csq
        gff_db (gffutils.FeatureDB): gffutils database
        force_indel_annotation(bool) : force indels to be annotated as indels (i.e. inframe, frameshift)

    Returns:
        pd.DataFrame: the VCF df with two new columns, one for the worst overall consequ
    """

    list_worst_csq = list()
    list_csqs = list()
    for idx, record in vcf_df.iterrows():

        # Get gene name(s)
        gene_id = record.INFO.split(";")[0].replace("GENE=", "")

        # TODO : Handle the case where this field is missing or has a different name, find a clever way to support all gff
        try:
            gene_names = gff_db[gene_id].attributes["gene_name"]
        except KeyError:
            gene_names = gff_db[gene_id].attributes["Name"]

        # Look at all the consequences returned by bcftoolscsq
        try:
            list_csq = record.INFO.split(";")[1].replace("BCSQ=", "").split(",")

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
                if force_indel_annotation & is_indel(record):
                    csq, csqs = assign_indel_csq(csq, csqs)

                # We keep the long bcftools string uniquely if it differs from the
                # single worst consequence
                if csqs == csq:
                    csqs = ""

            else:
                csq = ""
                csqs = ""

        except IndexError:
            # bcftoolscsq return sometimes empty consequences
            csq = ""
            csqs = ""

        list_worst_csq.append(csq)
        list_csqs.append(csqs)

    vcf_df["consequence"] = list_worst_csq
    vcf_df["full_consequence_bcftoolscsq"] = list_csqs

    return vcf_df


def is_indel(record):
    """
    Test whether the variant is an indel or a SNP

    Args:
        record (pd.Series): variant
    """

    if (len(record.REF) != 1) | (len(record.ALT) != 1):
        return True
    else:
        return False


def assign_indel_csq(csq, csqs):
    """
    Some indels are not called as inframe or frameshift. It means that no SNP-based score
    (CADD, dbNSFP...) can be assigned for them. In the simulation script, the user can decide
    to infer a score for those, but if he does not choose this option, these variants will be removed
    from the analysis.

    Here we offer an option to call indel indels even if bcftoolscsq returns a different consequence (e.g. splice_donor)

    Args:
        csq (str): the worst consequence found for this variant (e.g. splice_donor)
        csqs (str): the bcftoolscsq string including the worst consequence  (e.g. synonymous&splice_donor)
    """

    INDELS = ["inframe", "frameshift", "inframe_insertion", "inframe_deletion"]

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
@click.argument("gff_db")
@click.argument("out_rates")
@click.option(
    "--force-indel-annotation",
    is_flag=True,
    default=False,
)
def vcf_to_rates(vcf, rates, gff_db, out_rates, force_indel_annotation):
    """
    Merge the VCF file annotated with bcftoolscsq with a rate file (tabular)

    Args:
        vcf (str): Path of the VCF file
        rates (str): Path to the rates file
        gff_db (str) : Path to gffutils database
        out_rates (str): Path of the output rates file
        force_indel_annotation(bool) : force indels to be annotated as indels (i.e. inframe, frameshift)
    """

    init_log()
    logger = logging.getLogger("logger")

    # Read rates
    rates_df = pd.read_csv(rates, sep="\t", dtype={"chrom": str})

    # Read VCF
    try:
        if vcf.endswith("gz"):
            vcf_df = pd.read_csv(vcf, sep="\t", comment="#", header=None, compression="gzip")
        else:
            vcf_df = pd.read_csv(vcf, sep="\t", comment="#", header=None)
        vcf_df.columns = "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO".split("\t")
        vcf_df["#CHROM"] = vcf_df["#CHROM"].astype(str)

    except pd.errors.EmptyDataError:
        if rates_df.empty:
            logger.warning("Rates file is empty")
            rates_df["consequence"] = None
            rates_df.to_csv(out_rates, sep="\t", index=False)
            sys.exit(0)
        else:
            logger.error("Bcftools csq file is empty")
            sys.exit(1)

    # Load gffutils database
    gff_db = gffutils.FeatureDB(gff_db)

    # Extract BCFtools consequence in a separate column
    vcf_df = extract_worst_consequence(vcf_df, gff_db, force_indel_annotation)

    # Extract gene id to distinguish between overlapping genes
    vcf_df["gene_id"] = [x.split(";")[0].replace("GENE=", "") for x in vcf_df.INFO]

    # Merge files
    out_rates_df = rates_df.merge(
        vcf_df[["#CHROM", "POS", "REF", "ALT", "consequence", "full_consequence_bcftoolscsq", "gene_id"]],
        left_on=["chrom", "pos", "ref", "alt", "gene_id"],
        right_on=["#CHROM", "POS", "REF", "ALT", "gene_id"],
    )
    out_rates_df.drop(["#CHROM", "POS", "REF", "ALT"], axis=1, inplace=True)

    assert out_rates_df.shape[0] == rates_df.shape[0]

    # Export
    out_rates_df.to_csv(out_rates, sep="\t", index=False)


if __name__ == "__main__":
    cli()
