#!/usr/bin/env python

import logging
import os
import sys

import click
import gffutils
import pandas as pd
import pyfaidx
import utils

CDS_OFFSET = 50

def load_mutation_rate_model(mutation_rate_model_file):
    """Load mutation rate model

    Args:
        mutation_rate_model_file (str): Path to a mutation rate model file.

    Returns:
        pd.DataFrame : The mutation rate model as a pandas data frame.

    Examples:
        The mutation rate model should follow the above format.
        from	to	mu_snp
        AAA	ACA	1.6681119889870802e-09
        AAA	AGA	2.9280578256688802e-09
        AAA	ATA	1.03370814817543e-09
    """
    mutation_rate_model = pd.read_csv(mutation_rate_model_file, sep="\t")

    mutation_rate_model.index = (
        mutation_rate_model["from"] + "_" + mutation_rate_model["to"]
    )

    return mutation_rate_model


def load_gene_list(conf, gff_db,  column=0):
    """Load the optionally user provided gene of interest list.

    Args:
        conf (str): configuration
        gff_db (gffutils database) : GFF utils database
        column (str) : column to use to extract the genes identifiers 

    Returns:
         list: List of genes identifiers.

    Examples:
        The gene list can be a simple list or a TSV file like this :
        symbol	ENSG	HGNC_ID
        TARDBP	ENSG00000120948.19	HGNC:11571
    """

    logger = logging.getLogger("logger")

    if "GENE_LIST" in conf.keys() :

        logger.info(f"Getting gene list from : {conf['GENE_LIST']}")

        gene_list_file = conf["GENE_LIST"]
        # Simple list
        if column == 0 :
            gene_list = list(pd.read_csv(gene_list_file, sep="\t", header=None).iloc[:, 0].values)
        # TSV file
        else :
            gene_list = list(pd.read_csv(gene_list_file, sep="\t").loc[:, column].values)
    else :

        logger.info(f"Getting all genes in : {conf['GFF']}")
        gene_list = list()
        for gene in gff_db.all_features(featuretype= "gene", order_by = "start") :
            gene_list.append(gene.attributes["ID"][0])

    return gene_list


def load_gff(gff_file, gff_db_path):
    """Create the gff database used by gffutils.

    Args:
        gff_file (str): Path to a GFF or a gffutils database file.
        gff_db_path (str): Path to the newly created gff database.

    Returns:
        gffutils.db: GFF database.

    """

    logger = logging.getLogger("logger")

    # gffutils db input
    if gff_file.endswith(".db") :
        logger.info(f"Loading gffutils database : {gff_file}")
        gff_db = gffutils.FeatureDB(gff_file)
    # GFF input
    else :
        try:
            os.remove(gff_db_path)
            logger.info(f"Removed old gffutil database : {gff_db_path}")
        except OSError:
            pass

        logger.info(f"Creating gffutil database : {gff_db_path}")
        gff_db = gffutils.create_db(gff_file, gff_db_path, merge_strategy="create_unique")

    return gff_db


def get_alternates(seq, range_model):
    """ For a given kmer, returns all alternate codons with a different central nucleotide.

    Args:
        seq (str): kmer string
        range_model (int): length of sequence on each side of the central nucleotide

    Returns:
        list: list of alternates codons
    """

    # We generate all possible SNPs
    list_alternates = list()
    for nucleotide in ["A", "C", "T", "G"]:
        alternate_seq = seq[:range_model] + nucleotide + seq[range_model + 1 :]
        list_alternates.append(alternate_seq)

    # Remove the WT sequence from the list
    list_alternates = list(set(list_alternates) - set([seq]))
    assert len(list_alternates) == 3

    return list_alternates

def calculate_rates_cds(gene, start, mutation_rate_model, seq, range_model) :
    """Calculate the rate of each possiblle single nucleotide mutation in a given sequence according to a mutation rate model

    Args:
        gene (gffutils.feature.Feature): Gene object from gff db
        start (int): Genomic coordinate of the start of the sequence
        mutation_rate_model (pd.DataFrame): Mutation rate model
        seq (str): Sequence of interest
        range_model (int): Length of sequence on each side of the central nucleotide

    Returns:
        list: List of mutation rates
    """


    # Get reverse complement if gene on reverse strand
    if gene.strand == "-" :
        rev_seq = pyfaidx.complement(seq[::-1])
        cur_seq = rev_seq
    else :
        cur_seq = seq

    # Get the length of the current sequence
    len_seq = len(seq)
    
    list_mutations = list()
    # For each nucleotide in the sequence
    for i in range(range_model, len(cur_seq) - range_model):

        # We get the kmer centered on the current nucleotide
        ref = str(cur_seq[i - range_model : i + range_model + 1])

        # We generate all three possible alternate kmers changing the central nucleotide
        list_alternates = get_alternates(ref, range_model)

        # For the three alternate kmers we compute the mutation rate
        for alt in list_alternates:
            mutation_rate = mutation_rate_model.loc[ref + "_" + alt, "mu_snp"]

            ref_nuc = ref[range_model]
            alt_nuc = alt[range_model]

            # We update the genomic coordinates differently depending on the strand
            if gene.strand == "-" :
                ref_nuc  = pyfaidx.complement(ref_nuc)
                alt_nuc =  pyfaidx.complement(alt_nuc)
                pos = start + len_seq - 1 - i
            else :
                pos = start + i

            # We build a Serie for each possible mutation at each loci
            s = pd.Series(
                {
                    "symbol": gene.id,
                    "chrom": gene.chrom,
                    "pos": pos,
                    "ref": ref_nuc,
                    "alt": alt_nuc,
                    "prob": mutation_rate,
                }
            )

            # And append it to a list that we will use to build a data frame
            list_mutations.append(s)
        
    return list_mutations


def get_sequence(fasta, chrom, start, end ) :
    """ Returns sequence crom a fasta file according to genomic coordinates

    Args:
        fasta (pyfaidx.Fasta): Genome sequence
        chrom (string): chromosome id
        start (int): sequence start position
        end (int): sequence end position

    Returns:
        str: DNA sequence
    """

    return fasta[chrom][start:end]

def calculate_rates(mutation_rate_model, fasta, gff_db, gene_list):
    """Assign mutation rates to every considered loci using the mutation rate model.

    Args:
        mutation_rate_model (pd.DataFrame): Mutation rate for every possible kmer transition
        fasta (pyfaidx.Fasta): Genome sequence
        gff_db (gffutils.db): GFF database
        gene_list (list) : list of genes to consider

    Returns:
        pd.DataFrame: a mutation rate data frame
    """
    
    logger = logging.getLogger("logger")

    # Default mutation model is 3-mer, but other models can be used, hence we get the length here
    length_model = len(mutation_rate_model.iloc[0]["from"])

    # This variable gives the number of central nucleotide neighboring nucleotides to consider when computing rates
    range_model = length_model // 2

    # Get the number of genes in GFF
    #nb_genes = gff_db.count_features_of_type("gene")
    nb_genes = len(gene_list)
    logger.info(f"Start computing rates for {nb_genes} genes")

    list_mutation_rates = list()
    cpt = 1
    # Loop through all CDS of interest to generate mutation rates
    for gene in gff_db.all_features(featuretype="gene", order_by="start"):

        if gene.attributes["ID"][0] not in gene_list :
            continue
        
        logger.info(gene.attributes["ID"][0])

        list_mutation_rates_gene = list()
        for transcript in gff_db.children(gene, featuretype="transcript", order_by="start"):
            for cds in gff_db.children(transcript, featuretype="CDS", order_by="start"):

                # We add an offset to CDS region to consider loci such as splicing sites (default is 50),
                # the range model so that we get the neighboring nucleotides and we correct for python 0-index
                start = cds.start - CDS_OFFSET - range_model - 1
                end = cds.stop + CDS_OFFSET + range_model

                # We extract the coding sequence
                cds_seq = get_sequence(fasta, gene.chrom, start, end)
                #cds_seq = fasta[gene.chrom][start:end]
        
                # We calculate the rates for the current CDS region
                list_mutation_rates_cds = calculate_rates_cds(gene, start + 1, mutation_rate_model, cds_seq, range_model)
                list_mutation_rates_gene += list_mutation_rates_cds

            
            # TODO : We are considering only the first transcript for now
            break

        # For each gene we build a data frame and remove possible duplicated values
        mutation_rates_gene_df = pd.DataFrame(list_mutation_rates_gene)
        mutation_rates_gene_df.drop_duplicates(inplace=True, keep='first')
        list_mutation_rates.append(mutation_rates_gene_df)

        # Log progress
        if cpt % 100 == 0 :
            logger.info(f"{cpt}/{nb_genes} done")

        cpt += 1


    if (cpt - 1) != nb_genes :
        logger.warning(f"Rates computed for {cpt - 1} genes when genes in list is {nb_genes} ")
    else :
        logger.info(f"Rates computed for {nb_genes} genes")
    
    # We assemble all genes data frame together
    mutation_rates_df = pd.concat(list_mutation_rates)

    return mutation_rates_df


@click.command()
@click.option("--config")
@click.option("--gff")
@click.option("--fasta")
@click.option("--mutation_rate_model")
@click.option("--gene_list")
@click.option("--outdir")
def main(config, gff, fasta, mutation_rate_model, gene_list, outdir):
    """Generate mutation rates for genetic loci according to a mutation rate model.

    Args:
        conf_file (str): Path to a YAML configuration file.
    """

    # Initiate logger
    utils.init_log()
    logger = logging.getLogger("logger")
    logger.info("Running {}".format(__file__.split('/')[-1]))

    # Load configuration file
    if config : 
        try :
            conf = utils.load_conf(config)
        except FileNotFoundError :
            logger.error(f"No such config file {config}")
            sys.exit(1)
    else :
        conf = dict()
    
    # Superseed configuration in config file with the configuration passed through command line arguments
    ctx = click.get_current_context()
    conf = utils.superseed_conf(conf, ctx.params)

    # Log configuration
    logger.info(f"Parameters :")
    logger.info("----------")
    for key, value in conf.items() :
         logger.info(f"{key} : {value}")
    logger.info("----------")


    # Create output directory
    os.makedirs(conf["OUTDIR"], exist_ok=True)

    # Load mutation rate model
    mutation_rate_model = load_mutation_rate_model(conf["MUTATION_RATE_MODEL"])

    # Load GFF file or GFF database used by gffutils
    gff_db = load_gff(conf["GFF"], f'{conf["OUTDIR"]}/gff.db')

    # Load gene list
    gene_list = load_gene_list(conf, gff_db)    

    # Load fasta file
    fasta = pyfaidx.Fasta(conf["FASTA"])

    # Calculate mutation rates for every loci of interest
    mutation_rates_df = calculate_rates(mutation_rate_model, fasta, gff_db, gene_list)

    # Export mutation rates file
    mutation_rates_df.to_csv("{}/{}".format(conf["OUTDIR"], "mutation_rates.csv"))
    logger.info(f'Mutation rates file created here : {conf["OUTDIR"]}/mutation_rates.csv')

if __name__ == "__main__":

    main()
