#!/usr/bin/env python
import pandas as pd
import click


@click.command()
@click.argument("rates")
@click.argument("fasta")
@click.argument("out_vcf")
def rates_to_vcf(rates, fasta, out_vcf):
    """Transforms a rate file (tabular format) in VCF format

    Args:
        rates (file): Rate file
        fasta (file): FASTA file or fasta index file used to compute the rate file
        out_vcf (str): Destination of the output VCF file
    """

    # Retrieve contig length from FASTA file
    if fasta.endswith(".fai") :
        vcf_header_df = pd.read_csv(f"{fasta}", sep="\t", header=None)
    else :
        vcf_header_df = pd.read_csv(f"{fasta}.fai", sep="\t", header=None)
      
    
    # Write it in the VCF header format
    with open(out_vcf, "w") as f :
        f.write("##fileformat=VCFv4.2\n")
        for idx, row in vcf_header_df.iterrows() :
            chrom = row[0]
            length = row[1]
            f.write(f"##contig=<ID={chrom},length={length}>\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        
    # Extract information from rates file and transform it in VCF format
    rates_df = pd.read_csv(rates, sep="\t")
    with open(out_vcf, "a") as f :
        for idx, row in rates_df.iterrows() :
            f.write(f"{row.chrom}\t{row.pos}\t.\t{row.ref}\t{row.alt}\t.\t.\t.\n")
    


if __name__ == "__main__":
    rates_to_vcf()
