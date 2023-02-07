#!/usr/bin/env python
import pandas as pd
import click


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
    with open(out_vcf, "a") as f:
        for idx, row in rates_df.iterrows():
            f.write(f"{row.chrom}\t{row.pos}\t.\t{row.ref}\t{row.alt}\t.\t.\t.\n")


@cli.command()
@click.argument("vcf")
@click.argument("rates")
@click.argument("out_rates")
def vcf_to_rates(vcf, rates, out_rates):
    """Merge the VCF file annotated with bcftoolscsq with a rate file (tabular)
    Args:
        vcf (str): Path of the VCF file
        rates (str): Path to the rates file
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

    # Extract BCFtools consequence in a separate column
    vcf_df["consequence"] = [x.split("|")[0].replace("BCSQ=", "") for x in vcf_df.INFO]

    # Merge files
    out_rates_df = rates_df.merge(
        vcf_df[["#CHROM", "POS", "REF", "ALT", "consequence"]],
        left_on=["chrom", "pos", "ref", "alt"],
        right_on=["#CHROM", "POS", "REF", "ALT"],
    )
    out_rates_df.drop(["#CHROM", "POS", "REF", "ALT"], axis=1, inplace=True)

    assert out_rates_df.shape[0] == rates_df.shape[0]

    # Export
    out_rates_df.to_csv(out_rates, sep="\t", index=False)


if __name__ == "__main__":
    cli()
