#!/usr/bin/env python
import click
import gffutils
import pandas as pd

# Built from bcftools +split-vep -S -
consequences_severities = {
    "intergenic": 1,
    "feature_truncation": 2,
    "feature_elongation": 2,
    "regulatory": 3,
    "TF_binding_site": 4,
    "TFBS": 4,
    "downstream": 5,
    "upstream": 5,
    "non_coding_transcript": 6,
    "non_coding": 6,
    "intron": 7,
    "NMD_transcript": 7,
    "non_coding_transcript_exon": 8,
    "5_prime_utr": 9,
    "3_prime_utr": 9,
    "coding_sequence": 10,
    "mature_miRNA": 10,
    "stop_retained": 11,
    "start_retained": 11,
    "synonymous": 11,
    "incomplete_terminal_codon": 12,
    "splice_region": 13,
    "missense": 14,
    "inframe": 14,
    "protein_altering": 14,
    "transcript_amplification": 15,
    "exon_loss": 16,
    "disruptive": 17,
    "start_lost": 18,
    "stop_lost": 18,
    "stop_gained": 18,
    "frameshift": 18,
    "splice_acceptor": 19,
    "splice_donor": 19,
    "transcript_ablation": 20,
}


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


def extract_consequence(x, gff_db):

    # Get gene name(s)
    gene_id = x.split(";")[0].replace("GENE=", "")
    gene_names = gff_db[gene_id].attributes["gene_name"]

    # gene_transcripts_id = list()
    # for transcript in gff_db.children(gff_db[gene_id], featuretype="transcript"):
    #     gene_transcripts_id.append(transcript.id)

    # Look at all the consequences returned by bcftoolscsq
    try:
        list_csq = x.split(";")[1].replace("BCSQ=", "").split(",")

        # If there are multiple consequences, we take the worst one that
        # corresponds to our gene of interest
        # (we have multiple consequences when using multiple transcripts gff or even when using single transcripts like MANE, within splice regions)
        worst_csq_value = 1
        worst_csq = ""
        for cur_csq in list_csq:
            if cur_csq.split("|")[1] in gene_names:
                csqs = cur_csq.split("|")[0]
                for csq in csqs.split("&"):
                    if consequences_severities[csq] > worst_csq_value:
                        worst_csq_value = consequences_severities[csq]
                        worst_csq = csqs

        if worst_csq:
            csq = worst_csq
        else:
            csq = ""

    except IndexError as e:
        # bcftoolscsq return sometimes empty consequences
        csq = ""

    return csq


@cli.command()
@click.argument("vcf")
@click.argument("rates")
@click.argument("gff_db")
@click.argument("out_rates")
def vcf_to_rates(vcf, rates, gff_db, out_rates):
    """Merge the VCF file annotated with bcftoolscsq with a rate file (tabular)
    Args:
        vcf (str): Path of the VCF file
        rates (str): Path to the rates file
        gff_db (str) : Path to gffutils database
        out_rates (str): Path of the output rates file
    """

    # Read VCF
    if vcf.endswith("gz"):
        vcf_df = pd.read_csv(vcf, sep="\t", comment="#", header=None, compression="gzip")
    else:
        vcf_df = pd.read_csv(vcf, sep="\t", comment="#", header=None)
    vcf_df.columns = "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO".split("\t")

    # Read rates
    rates_df = pd.read_csv(rates, sep="\t")

    # Load gffutils database
    gff_db = gffutils.FeatureDB(gff_db)

    # Extract BCFtools consequence in a separate column
    vcf_df["consequence"] = [extract_consequence(x, gff_db) for x in vcf_df.INFO]

    # Extract gene id to distinguish between overlapping genes
    vcf_df["gene_id"] = [x.split(";")[0].replace("GENE=", "") for x in vcf_df.INFO]

    # Merge files
    out_rates_df = rates_df.merge(
        vcf_df[["#CHROM", "POS", "REF", "ALT", "consequence", "gene_id"]],
        left_on=["chrom", "pos", "ref", "alt", "gene_id"],
        right_on=["#CHROM", "POS", "REF", "ALT", "gene_id"],
    )
    out_rates_df.drop(["#CHROM", "POS", "REF", "ALT"], axis=1, inplace=True)

    assert out_rates_df.shape[0] == rates_df.shape[0]

    # Export
    out_rates_df.to_csv(out_rates, sep="\t", index=False)


if __name__ == "__main__":
    cli()
