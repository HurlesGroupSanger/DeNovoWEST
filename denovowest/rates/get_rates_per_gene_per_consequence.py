#!/usr/bin/env python

import click
import pandas as pd

LIST_CSQ = [
    "synonymous",
    "missense",
    "start_lost",
    "stop_lost",
    "stop_retained",
    "stop_gained",
    "splice_acceptor",
    "splice_donor",
    "splice_region",
]


def iter_rates(rates):

    reader = pd.read_csv(rates, sep="\t", chunksize=100000, low_memory=False)

    for chunk in reader:

        rates_dfs = chunk.groupby(chunk.gene_id, sort=False)
        gene_ids = list(rates_dfs.groups.keys())

        if not gene_ids:
            continue

        # Yield complete groups
        for gene_id in gene_ids:
            yield gene_id, rates_dfs.get_group(gene_id)


@click.command()
@click.argument("rates", type=str)
@click.option("--output", type=str, default="rates_per_gene_per_consequence.tsv", show_default=True)
@click.option("--output-nb", type=str, default="nb_records_per_gene_per_consequence.tsv", show_default=True)
def main(rates, output, output_nb):
    """
    Compute expected number of mutation per gene per consequence type

    Arguments
    ---------
    RATES (TSV):
        A tab-delimited file with at least the following three columns

            gene_id : gene identifier
            consequence : consequence from bcftoolscsq
            prob : mutation probability

    OUTPUT:
        Where to write the output matrix with expected number of mutations per gene per consequence

    OUTPUT_NB:
        Where to write the output matrix with number of records per gene per consequence
    """

    prob_per_gene_per_consequence = dict()
    nb_mutations_per_gene_per_consequence = dict()

    cpt = 0
    for gene_id, rates_df in iter_rates(rates):

        gene_already_seen = False

        # Sometimes the same gene appears multiple times in the chunked reading (e.g. end of chunk, overlapping genes)
        if gene_id in prob_per_gene_per_consequence.keys():
            gene_already_seen = True
            prob_per_consequence = prob_per_gene_per_consequence[gene_id]
            nb_mutations_per_consequence = nb_mutations_per_gene_per_consequence[gene_id]
        else:
            prob_per_consequence = dict()
            nb_mutations_per_consequence = dict()

        # We keep track of the sum of probabilities per consequence
        for csq in LIST_CSQ:
            # Sometimes the same gene appears multiple times in the chunked reading (e.g. end of chunk, overlapping genes)
            if gene_already_seen:
                prob_per_consequence[csq] += rates_df.loc[rates_df.consequence == csq, "prob"].sum()
                nb_mutations_per_consequence[csq] += rates_df.loc[rates_df.consequence == csq].shape[0]

            else:
                prob_per_consequence[csq] = rates_df.loc[rates_df.consequence == csq, "prob"].sum()
                nb_mutations_per_consequence[csq] = rates_df.loc[rates_df.consequence == csq].shape[0]

        prob_per_gene_per_consequence[gene_id] = prob_per_consequence
        nb_mutations_per_gene_per_consequence[gene_id] = nb_mutations_per_consequence

        cpt += 1
        if cpt % 100 == 0:
            print(f"Processed {cpt} genes")

    # Write the number of expected mutations per gene per consequence
    res_df = pd.DataFrame(prob_per_gene_per_consequence).transpose()
    res_df.to_csv(output, sep="\t")

    # Write the number of records per gene per consequence
    count_mutations_df = pd.DataFrame(nb_mutations_per_gene_per_consequence).transpose()
    count_mutations_df.to_csv(output_nb, sep="\t")


if __name__ == "__main__":
    main()
