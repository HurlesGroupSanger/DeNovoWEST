#!/usr/bin/env python

import pandas as pd
import click
from scipy.stats import chi2


def set_constraints_loci(region, gene_rates, threshold, ratio):
    """_summary_

    Args:
        region (pd.Series): Informations about a constrained region
    """

    # Compute the probability that a random variable with chi-square distribution is greater than chisq_diff_null
    p_value = chi2.sf(region["chisq_diff_null"], df=1)

    # Leave constrained loci status to False if p_value or observed/expected ratio above thresholds
    if (p_value > threshold) or (region["obs_exp"] > ratio):
        return

    # Set constrained status to True to every loci in the region which has as missense consequence
    gene_rates.loc[
        (gene_rates.pos >= region.cds_start)
        & ((gene_rates.pos <= region.cds_end))
        & ((gene_rates.consequence == "missense")),
        "constrained",
    ] = True


@click.command()
@click.argument("rates")
@click.argument("geneconstraint")
@click.argument("regionconstraint")
@click.argument("outrates")
def annotate_constraint(rates, geneconstraint, regionconstraint, outrates, threshold=1e-3, ratio=0.4):
    """Annotate rates file with the constrained status for each missense mutation loci in constrained regions

    Args:
        rates (str): path to rates file (already annotated with consequence)
        geneconstraint (str): 18K set inc different tx
        regionconstraint (str): >3k genes single entry
        outrates (str) : output path to constrained annotated rates file
        threshold (float, optional): _description_. Defaults to 1e-3.
        ratio (float, optional): _description_. Defaults to 0.4.

    Returns:
        _type_: _description_
    """

    # Read constraints files
    geneconstraint_df = pd.read_table(geneconstraint)
    regionconstraint_df = pd.read_table(regionconstraint)

    # Read rates file and initialize all positions as non constrained
    rates_df = pd.read_csv(rates, sep="\t")
    rates_df["constrained"] = False

    # For each gene in rates file
    for id_gene, gene_rates in rates_df.groupby("gene_id"):
        # We do not consider the version of the gene here
        ensembl_id = id_gene.split(".")[0]
        ensembl_id = ensembl_id.replace("gene:", "")

        # If the gene is not found in the constraint file, print as a warning
        if ensembl_id not in set(geneconstraint_df["ensembl_gene_id"]):
            print(f"Could not find id {ensembl_id} in full constraint file")
        # If the gene is found in the constraint file
        else:
            # Restrict the constraint data frame to current gene
            cur_geneconstraint_df = geneconstraint_df[geneconstraint_df["ensembl_gene_id"] == ensembl_id]
            assert cur_geneconstraint_df.shape[0] == 1

            # If there is only one constained region in this gene
            if cur_geneconstraint_df.iloc[0]["n_regions"] == 1:
                # Compute the probability that a random variable with chi-square distribution is greater than overall_chisq
                genep_value = chi2.sf(cur_geneconstraint_df["overall_chisq"], df=1)

                # If both this probability and the ratio of observed vs expected mutations is under their respective thresholds
                # we set all missense mutation loci constraint status to true for the current gene
                obs_exp_ratio = (
                    cur_geneconstraint_df.iloc[0]["obs_missense"] / cur_geneconstraint_df.iloc[0]["exp_missense"]
                )
                if (genep_value <= threshold) and (obs_exp_ratio <= ratio):
                    gene_rates.loc[gene_rates.consequence == "missense", "constrained"] = True
            # Otherwise if there are several constained regions in the gene
            else:
                # Restrict regional dataframe to current transcript
                transcript_id = cur_geneconstraint_df.iloc[0]["ensembl_transcript_enst"]
                cur_regionconstraint_df = regionconstraint_df[
                    regionconstraint_df["ensembl_transcript_enst"] == transcript_id
                ]

                # Set constrained status for each region
                for idx, region in cur_regionconstraint_df.iterrows():
                    set_constraints_loci(region, gene_rates, threshold, ratio)

            # Modify the rates data frame with constraints status
            rates_df.loc[gene_rates.index, "constrained"] = gene_rates.constrained

        # Export rates file (updated with constrained status)
        rates_df.to_csv(outrates, sep="\t", index=False)


if __name__ == "__main__":
    annotate_constraint()
