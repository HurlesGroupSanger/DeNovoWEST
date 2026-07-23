#!/usr/bin/env python
import click
import gffutils
import pandas as pd

consequences_dict = {"missense": "missense_variant", "missense_variant": "missense_variant"}


def snp_or_indel(row):
    if (len(row.ref) - len(row.alt)) != 0:
        return "DENOVO-INDEL"
    else:
        return "DENOVO-SNP"


def handle_duplicated_gene_symbols(dnm_df):
    """
    Some gene symbols are associated with multiple ENSEMBL gene identifiers (e.g. DUSP13B).
    DNN works with gene symbols and ends up with duplicated mutations for those genes, which influences the results.
    This function retains DNMs from the gene identifier with the most DNMs and discard the others.

    #TODO : This solution is far from perfect but very few genes are affected by this issue. A better solution would be to use the gene identifier instead of the gene symbol in DNN.
    """

    duplicated_gene_symbols = list()
    for gene_symbol, min_df in dnm_df.groupby("gene_name"):
        if min_df.gene_id.nunique() != 1:
            duplicated_gene_symbols.append(gene_symbol)

    for gene_symbol in duplicated_gene_symbols:
        min_df = dnm_df[dnm_df.gene_name == gene_symbol]
        gene_id_to_keep = min_df.gene_id.value_counts().idxmax()
        gene_id_discarded = set(min_df.gene_id.unique()) - {gene_id_to_keep}

        dnm_df = dnm_df[~((dnm_df.gene_name == gene_symbol) & (dnm_df.gene_id != gene_id_to_keep))]

        print(
            f"Gene symbol {gene_symbol} is associated with multiple gene identifiers. Retaining DNMs from gene identifier {gene_id_to_keep} and discarding the others ({','.join(list(gene_id_discarded))})."
        )

    return dnm_df


@click.command()
@click.argument("dnms", type=click.Path(exists=True))
@click.argument("gff", type=click.Path(exists=True))
@click.argument("out", type=str)
def prepare_dnms(dnms, gff, out):
    """
    Prepare the DNM file in order to be used with DeNovoNear

    Args:
        dnms (str): DNM file
        gff (str): GFF or gffutils database
        out(str): DeNovoNear ready DNM file

    """

    # Load DNM file
    dnm_df = pd.read_csv(dnms, sep="\t")

    # Load GFF
    if not gff.endswith("db"):
        gff_db = gffutils.create_db(gff, f"{gff}.db", merge_strategy="create_unique")
    else:
        gff_db = gffutils.FeatureDB(gff)

    # DeNovoNear uses gene name rather than gene identifiers
    gene_dict = dict()
    for gene in gff_db.features_of_type("gene", order_by="start"):
        gene_dict[gene["ID"][0]] = gene["Name"][0]
    dnm_df["gene_name"] = dnm_df["gene_id"].map(gene_dict)

    # Handle duplicated gene symbols
    dnm_df = handle_duplicated_gene_symbols(dnm_df)

    # DeNovoNear has a column that stats whether the mutation is an indel or a SNP
    dnm_df["snp_or_indel"] = dnm_df.apply(snp_or_indel, axis=1)

    # DeNovoNear expects missense variants to be flagged as "missense_variant"
    dnm_df["consequence"] = dnm_df["consequence"].map(consequences_dict)

    # Format the data frame
    dnm_df.rename({"chrom": "chr"}, axis=1, inplace=True)
    dnm_df = dnm_df[["gene_name", "chr", "pos", "consequence", "snp_or_indel"]]

    # Export DeNovoNear ready DNM file
    dnm_df.to_csv(out, sep="\t", index=False)


if __name__ == "__main__":
    prepare_dnms()
