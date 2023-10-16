#!/usr/bin/env python

import logging
import os
import sys

import click
import gffutils
import pandas as pd
import pyfaidx
import pysam
import utils
import glob

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

    mutation_rate_model.index = mutation_rate_model["from"] + "_" + mutation_rate_model["to"]

    return mutation_rate_model


def load_gene_list(conf, gff_db, column=0):
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

    if "GENE_LIST" in conf.keys():
        logger.info(f"Getting gene list from : {conf['GENE_LIST']}")

        gene_list_file = conf["GENE_LIST"]
        # Simple list
        if column == 0:
            gene_list = list(pd.read_csv(gene_list_file, sep="\t", header=None).iloc[:, 0].values)
        # TSV file
        else:
            gene_list = list(pd.read_csv(gene_list_file, sep="\t").loc[:, column].values)
    else:
        logger.info(f"Getting all genes in : {conf['GFF']}")
        gene_list = list()
        for gene in gff_db.all_features(featuretype="gene", order_by="start"):
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
    if gff_file.endswith(".db"):
        logger.info(f"Loading gffutils database : {gff_file}")
        gff_db = gffutils.FeatureDB(gff_file)
    # GFF input
    else:
        try:
            os.remove(gff_db_path)
            logger.info(f"Removed old gffutil database : {gff_db_path}")
        except OSError:
            pass

        logger.info(f"Creating gffutil database : {gff_db_path}")
        gff_db = gffutils.create_db(gff_file, gff_db_path, merge_strategy="create_unique")

    return gff_db


def get_alternates(seq, range_model):
    """
    For a given kmer, returns all alternative kmers with a different central nucleotide.

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

    return list_alternates


def calculate_rates_cds(gene, start, mutation_rate_model, seq, range_model):
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
    if gene.strand == "-":
        rev_seq = pyfaidx.complement(seq[::-1])
        cur_seq = rev_seq
    else:
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
            try:
                mutation_rate = mutation_rate_model.loc[ref + "_" + alt, "mu_snp"]
            except KeyError as e:
                # Ambiguous nucleotides
                continue

            ref_nuc = ref[range_model]
            alt_nuc = alt[range_model]

            # We update the genomic coordinates differently depending on the strand
            if gene.strand == "-":
                ref_nuc = pyfaidx.complement(ref_nuc)
                alt_nuc = pyfaidx.complement(alt_nuc)
                pos = start + len_seq - 1 - i
            else:
                pos = start + i

            # We build a Serie for each possible mutation at each loci
            s = pd.Series(
                {
                    "gene_id": gene.id,
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


def get_sequence(fasta, chrom, start, end):
    """Returns sequence from a fasta file according to genomic coordinates

    Args:
        fasta (pyfaidx.Fasta): Genome sequence
        chrom (string): chromosome id
        start (int): sequence start position
        end (int): sequence end position

    Returns:
        str: DNA sequence
    """

    return fasta[chrom][start:end]


def calculate_rates_kmer(mutation_rate_model, fasta, gff_db, gene_list):
    """Assign mutation rates to every considered loci using kmer mutation rate model.
    https://www.nature.com/articles/ng.3050


    Args:
        mutation_rate_model (pd.DataFrame): Mutation rate for every possible kmer transition
        fasta (pyfaidx.Fasta): Genome sequence
        gff_db (gffutils.db): GFF database
        gene_list (list) : list of genes to consider

    Returns:
        pd.DataFrame: a mutation rate data frame
    """

    logger = logging.getLogger("logger")

    # Default mutation model is 3-mer, but other models (5-mers, 7-mers...) can be used, hence we get the length here
    length_model = len(mutation_rate_model.iloc[0]["from"])

    # Number of neighboring nucleotides to consider when computing rates
    range_model = length_model // 2

    nb_genes = len(gene_list)
    logger.info(f"Start computing rates for {nb_genes} genes")

    list_mutation_rates = list()
    cpt = 1
    # Loop through all CDS of interest to generate mutation rates
    for gene in gff_db.all_features(featuretype="gene"):
        # Skip gene if not in user provided gene list
        gene_id = gene.attributes["ID"][0]
        if gene_id not in gene_list:
            continue

        # Store mutation rates for the current gene
        list_mutation_rates_gene = list()
        # Store CDS boundaries to avoid calculating rates for same CDS in several transcripts
        list_cds_boundaries = list()
        for transcript in gff_db.children(gene, level=1):
            for cds in gff_db.children(transcript, featuretype="CDS", order_by="start"):
                # We add an offset to CDS region to consider loci such as splicing sites (default is 50),
                # the range model so that we get the neighboring nucleotides and we correct for python 0-index
                start = cds.start - CDS_OFFSET - range_model - 1
                end = cds.stop + CDS_OFFSET + range_model

                # If the CDS has already been covered by another transcript we skip it
                if (start, end) in list_cds_boundaries:
                    continue
                else:
                    list_cds_boundaries.append((start, end))

                # We extract the coding sequence
                cds_seq = get_sequence(fasta, gene.chrom, start, end)
                # cds_seq = fasta[gene.chrom][start:end]

                # We calculate the rates for the current CDS region
                list_mutation_rates_cds = calculate_rates_cds(
                    gene, start + 1, mutation_rate_model, cds_seq, range_model
                )
                list_mutation_rates_gene += list_mutation_rates_cds

        # For each gene we build a data frame and remove possible duplicated values
        if list_mutation_rates_gene:
            mutation_rates_gene_df = pd.DataFrame(list_mutation_rates_gene)
            mutation_rates_gene_df.drop_duplicates(inplace=True, keep="first")
            mutation_rates_gene_df.sort_values(by=["pos", "alt"], inplace=True)
            list_mutation_rates.append(mutation_rates_gene_df)
        else:
            logger.warning(f"No mutations found for gene {gene_id}")

        # Log progress
        if cpt % 100 == 0:
            logger.info(f"{cpt}/{nb_genes} done")

        cpt += 1

    if (cpt - 1) != nb_genes:
        logger.warning(f"Rates computed for {cpt - 1} genes when number of genes in list is {nb_genes} ")
    else:
        logger.info(f"Rates computed for {nb_genes} genes")

    # We assemble all genes data frame together
    mutation_rates_df = pd.concat(list_mutation_rates)

    return mutation_rates_df


def roulette_per_chrom_files(roulette_dir):
    """
    Build a dictionnary with chromosome as key and path to roulette vcf file as value

    Args:
        roulette_dir (str): path to directory containing roulette vcf files per chromosome

    Returns:
        dict: a dictionnary with chromosome as key and path to roulette vcf file as value
    """

    roulette_per_chrom_files_list = glob.glob(f"{roulette_dir}/*.vcf.bgz")
    roulette_per_chrom_files = dict()
    for f in roulette_per_chrom_files_list:
        chrom = f.split("/")[-1].split("_")[0]
        roulette_per_chrom_files[chrom] = f

    return roulette_per_chrom_files


def calculate_rates_roulette(roulette_dir, gff_db, gene_list, model):
    """
    Use Roulette mutation rate model to generate a rates file
    https://github.com/vseplyarskiy/Roulette
    https://www.biorxiv.org/content/10.1101/2022.08.20.504670v1

    The roulette directory contains one VCF file per chromosome.
    Loci in each VCF are annotated with several mutation rates :
        - MR : roulette mutation rate
        - MR : carlson mutation rate

    Args:
        roulette_dir (string): path to directory containing roulette vcf files per chromosome
        gff_db (gffutils.db): GFF database
        gene_list (list) : list of genes to consider
        model (str) : carlson or roulette

    Returns:
        pd.DataFrame: a mutation rate data frame
    """
    logger = logging.getLogger("logger")

    # Select the mutation rate info field and the per generation mutation rate scaling factor depending on the model
    if model == "roulette":
        info_field = "MR"
        scaling_factor = utils.ROULETTE_SCALING_FACTOR
    else:
        info_field = "MC"
        scaling_factor = utils.CARLSON_SCALING_FACTOR

    nb_genes = len(gene_list)
    logger.info(f"Start computing rates for {nb_genes} genes")

    # Build a dictionnary with chromosome as key and path to roulette vcf file as value
    roulette_vcfs = roulette_per_chrom_files(roulette_dir)

    list_mutation_rates = list()
    cpt = 1
    # Loop through all CDS of interest to generate mutation rates
    for gene in gff_db.all_features(featuretype="gene"):
        # Skip gene if not in user provided gene list
        gene_id = gene.attributes["ID"][0]
        if gene_id not in gene_list:
            continue

        # Load roulette vcf file corresponding to the current gene
        chrom = gene.chrom.replace("chr", "")
        try:
            roulette_file = pysam.VariantFile(roulette_vcfs[chrom], index_filename=f"{roulette_vcfs[chrom]}.csi")
        except KeyError as e:
            # Roulette does not provide mutation rates for allosomes
            if chrom in ["X", "Y"]:
                continue
            else:
                logger.warning(f"Can't find any roulette file corresponding to chromosome {chrom} for gene {gene_id}")

        # Store mutation rates for the current gene
        list_mutation_rates_gene = list()
        # Store CDS boundaries to avoid calculating rates for same CDS in several transcripts
        list_cds_boundaries = list()
        for transcript in gff_db.children(gene, level=1):
            for cds in gff_db.children(transcript, featuretype="CDS", order_by="start"):
                # We add an offset to CDS region to consider loci such as splicing sites (default is 50),
                start = cds.start - CDS_OFFSET - 1
                end = cds.stop + CDS_OFFSET

                # If the CDS has already been covered by another transcript we skip it
                if (start, end) in list_cds_boundaries:
                    continue
                else:
                    list_cds_boundaries.append((start, end))

                # Extract mutation rates for the current CDS region and multiply by the scaling factor
                list_mutation_rates_cds = list()
                for rec in roulette_file.fetch(start=start, stop=end, region=chrom):
                    try:
                        s = pd.Series(
                            {
                                "gene_id": gene_id,
                                "chrom": chrom,
                                "pos": rec.pos,
                                "ref": rec.ref,
                                "alt": rec.alts[0],
                                "prob": float(rec.info[info_field]) * scaling_factor,
                            }
                        )
                        list_mutation_rates_cds.append(s)
                    except KeyError as e:
                        # No mutation rate info for this SNP
                        continue

                list_mutation_rates_gene += list_mutation_rates_cds

        # For each gene we build a data frame and remove possible duplicated values
        if list_mutation_rates_gene:
            mutation_rates_gene_df = pd.DataFrame(list_mutation_rates_gene)
            mutation_rates_gene_df.drop_duplicates(inplace=True, keep="first")
            mutation_rates_gene_df.sort_values(by=["pos", "alt"], inplace=True)
            list_mutation_rates.append(mutation_rates_gene_df)
        else:
            logger.warning(f"No mutations found for gene {gene_id}")

        # Log progress
        if cpt % 100 == 0:
            logger.info(f"{cpt}/{nb_genes} done")

        cpt += 1

    if (cpt - 1) != nb_genes:
        logger.warning(f"Rates computed for {cpt - 1} genes when number of genes in list is {nb_genes} ")
    else:
        logger.info(f"Rates computed for {nb_genes} genes")

    # We assemble all genes data frame together
    mutation_rates_df = pd.concat(list_mutation_rates)

    return mutation_rates_df


def validate_model(ctx, param, value):
    """
    Validatin function called at the program start to check model is correct.
    Model can be either kmer, carlson or roulette

    Args:
        ctx (click.Context):  context
        param (click.Parameter): parameter
        value (str): value provided for the model

    Raises:
        click.BadParameter: raised if user provided a value not in the allowed values

    Returns:
        str: model value
    """
    allowed_models = ["kmer", "carlson", "roulette"]
    if value not in allowed_models:
        raise click.BadParameter(f'Invalid model. Allowed values: {", ".join(allowed_models)}')
    return value


@click.command()
@click.option("--config")
@click.option("--gff")
@click.option("--fasta")
@click.option("--mutation_rate_model", help="Path to a k-mer mutation rate model or roulette directory")
@click.option("--gene_list")
@click.option("--outdir")
@click.option("--model", type=click.STRING, callback=validate_model)
def main(config, gff, fasta, mutation_rate_model, gene_list, outdir, model):
    """
    Generate per generation mutation rates according to a mutation rate model.
    Loci to be considered are provided as a GFF file or a gffutils database file, and a gene list.

    Args:
        conf_file (str): Use a YAML configuration file instead of command line arguments.
        gff (str): Annotations provided as GFF or a gffutils database file.
        fasta (str): Genome sequence in fasta format to retrieve kmer sequences (kmer only)
        mutation_rate_model (str): Path to a k-mer mutation rate model or roulette directory
        gene_list (str): File containing a list of genes of interest.
        outdir (str): Output directory.
        model (str): Mutation rate model to be used. Allowed values: kmer, carlson or roulette
    """

    # Initiate logger
    utils.init_log()
    logger = logging.getLogger("logger")
    logger.info("Running {}".format(__file__.split("/")[-1]))

    # Load configuration file
    if config:
        try:
            conf = utils.load_conf(config)
        except FileNotFoundError:
            logger.error(f"No such config file {config}")
            sys.exit(1)
    else:
        conf = dict()

    # Superseed configuration in config file with the configuration passed through command line arguments
    ctx = click.get_current_context()
    conf = utils.superseed_conf(conf, ctx.params)

    # Log configuration
    logger.info(f"Parameters :")
    logger.info("----------")
    for key, value in conf.items():
        logger.info(f"{key} : {value}")
    logger.info("----------")

    # Create output directory
    os.makedirs(conf["OUTDIR"], exist_ok=True)

    # Load GFF file or GFF database used by gffutils
    gff_db = load_gff(conf["GFF"], f'{conf["OUTDIR"]}/gff.db')

    # Load gene list
    gene_list = load_gene_list(conf, gff_db)

    if model == "kmer":
        # Load mutation rate model
        mutation_rate_model = load_mutation_rate_model(conf["MUTATION_RATE_MODEL"])

        # Load fasta file
        fasta = pyfaidx.Fasta(conf["FASTA"])

        # Calculate mutation rates for every loci of interest
        mutation_rates_df = calculate_rates_kmer(mutation_rate_model, fasta, gff_db, gene_list)
    else:
        mutation_rates_df = calculate_rates_roulette(mutation_rate_model, gff_db, gene_list, model)

    # Export mutation rates file
    mutation_rates_df.to_csv("{}/{}".format(conf["OUTDIR"], "mutation_rates.tsv"), sep="\t", index=None)
    logger.info(f'Mutation rates file created here : {conf["OUTDIR"]}/mutation_rates.tsv')


if __name__ == "__main__":
    main()
