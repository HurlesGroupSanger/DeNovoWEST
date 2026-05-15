import os
import subprocess
from datetime import date
from pathlib import Path

import click


@click.command()
@click.argument("input_dir", type=click.Path(exists=True))
@click.argument("gene_input", type=str)
@click.argument("output_dir", type=click.Path())
@click.option(
    "--gene-file",
    "is_gene_file",
    is_flag=True,
    help="Interpret gene_input as a file containing one gene name per line.",
)
def create_single_gene_dataset(input_dir, gene_input, output_dir, is_gene_file):
    """
    Create a dataset containing only the specified gene from multiple input files.

    INPUT_DIR: Directory containing input files.
    GENE_NAME: Name of the gene to extract.
    OUTPUT_DIR: Path to the output directory to save the single gene dataset.

    """

    os.makedirs(output_dir, exist_ok=True)

    if is_gene_file:
        gene_path = Path(gene_input)

        if not gene_path.exists():
            raise click.ClickException(f"Gene file does not exist: {gene_path}")

        gene_list = [line.strip() for line in gene_path.read_text().splitlines() if line.strip()]
        if len(gene_list) <= 5:
            output_suffix = "_".join(gene_list)
        else:
            output_suffix = date.today().strftime("%Y-%m-%d")

        genes = "|".join(gene_list)

    else:
        genes = gene_input
        output_suffix = gene_input

    # Extract gene-specific data from DNM annotated file
    dnm_annotated_file = os.path.join(input_dir, "dnm/dnm_annotated.tsv")
    cmd_extract_header = f"head -n 1 {dnm_annotated_file} > {output_dir}/dnm_annotated_{output_suffix}.tsv"
    subprocess.run(cmd_extract_header, shell=True, check=True)

    cmd_extract_gene_data = f'grep -wE "{genes}" {dnm_annotated_file} >> {output_dir}/dnm_annotated_{output_suffix}.tsv'
    subprocess.run(cmd_extract_gene_data, shell=True, check=True)

    # Extract gene-specific data from rates file
    rates_file = os.path.join(input_dir, "rates/merged_rates.tsv.gz")
    cmd_extract_header = f"zcat {rates_file} | head -n 1 > {output_dir}/rates_annotated_{output_suffix}.tsv"
    subprocess.run(cmd_extract_header, shell=True, check=True)

    cmd_extract_gene_data = (
        f'zcat {rates_file} | grep -wE "{genes}" >> {output_dir}/rates_annotated_{output_suffix}.tsv'
    )
    subprocess.run(cmd_extract_gene_data, shell=True, check=True)

    cmd_bgzip = f"bgzip -f {output_dir}/rates_annotated_{output_suffix}.tsv"
    subprocess.run(cmd_bgzip, shell=True, check=True)


if __name__ == "__main__":
    create_single_gene_dataset()
