import click
import os
import subprocess


@click.command()
@click.argument("input_dir", type=click.Path(exists=True))
@click.argument("gene_name", type=str)
@click.argument("output_dir", type=click.Path())
def create_single_gene_dataset(input_dir, gene_name, output_dir):
    """
    Create a dataset containing only the specified gene from multiple input files.

    INPUT_DIR: Directory containing input files.
    GENE_NAME: Name of the gene to extract.
    OUTPUT_DIR: Path to the output directory to save the single gene dataset.
    """

    os.makedirs(output_dir, exist_ok=True)

    # Extract gene-specific data from DNM annotated file
    dnm_annotated_file = os.path.join(input_dir, "dnm/dnm_annotated.tsv")
    cmd_extract_header = f"head -n 1 {dnm_annotated_file} > {output_dir}/dnm_annotated_{gene_name}.tsv"
    subprocess.run(cmd_extract_header, shell=True, check=True)

    cmd_extract_gene_data = f"grep -w {gene_name} {dnm_annotated_file} >> {output_dir}/dnm_annotated_{gene_name}.tsv"
    subprocess.run(cmd_extract_gene_data, shell=True, check=True)

    # Extract gene-specific data from rates file
    rates_file = os.path.join(input_dir, "rates/merged_rates.tsv.gz")
    cmd_extract_header = f"zcat {rates_file} | head -n 1 > {output_dir}/rates_annotated_{gene_name}.tsv"
    subprocess.run(cmd_extract_header, shell=True, check=True)

    cmd_extract_gene_data = f"zcat {rates_file} | grep -w {gene_name} >> {output_dir}/rates_annotated_{gene_name}.tsv"
    subprocess.run(cmd_extract_gene_data, shell=True, check=True)

    cmd_bgzip = f"bgzip -f {output_dir}/rates_annotated_{gene_name}.tsv"
    subprocess.run(cmd_bgzip, shell=True, check=True)


if __name__ == "__main__":
    create_single_gene_dataset()
