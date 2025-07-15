import click
from denovowest.simulation import simulation
from denovowest.annotations import annotate_cadd, annotate_custom, annotate_dbnsfp, annotate_vcf


@click.group()
def main():
    """CLI for denovowest"""
    pass


# Register subcommands
main.add_command(simulation.main, name="simulation")
main.add_command(annotate_cadd.annotate_cadd, name="annotate-cadd")
main.add_command(annotate_custom.cli, name="annotate-custom")
main.add_command(annotate_dbnsfp.annotate_dbnsfp, name="annotate-dbnsfp")
main.add_command(annotate_vcf.cli, name="annotate-vcf")
