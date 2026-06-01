#!/usr/bin/env python

import logging
from pathlib import Path

import click
import gffutils
import pandas as pd
import pyfaidx
import pysam

from denovowest.utils.io_helpers import load_gff
from denovowest.utils.log import init_log
from denovowest.utils.params import CDS_OFFSET

#############################
# Data loading              #
#############################


def load_kmer_mutation_rate_model(mutation_rate_model_file: Path) -> pd.DataFrame:
    """Load a kmer mutation rate model from a tab-separated file.

    Args:
        mutation_rate_model_file: Path to the mutation rate model file.

    Returns:
        The mutation rate model as a pandas DataFrame indexed by
        ``<from>_<to>`` kmer transitions.

    Note:
        Expected format (tab-separated)::

            from    to      mu_snp
            AAA     ACA     1.668e-09
    """

    mutation_rate_model = pd.read_csv(mutation_rate_model_file, sep="\t")
    mutation_rate_model.index = mutation_rate_model["from"] + "_" + mutation_rate_model["to"]
    return mutation_rate_model


def load_gene_list(gene_list: Path | None, gff_db: gffutils.FeatureDB) -> list[str]:
    """Load the gene-of-interest list.

    Args:
        gene_list: Path to a plain text file with one gene identifier per line.
            When ``None``, all genes present in ``gff_db`` are returned.
        gff_db: GFF utils database (used when ``gene_list`` is ``None``).

    Returns:
        List of gene identifiers.
    """

    logger = logging.getLogger("logger")

    if gene_list is not None:
        logger.info(f"Getting gene list from : {gene_list}")
        return list(pd.read_csv(gene_list, sep="\t", header=None).iloc[:, 0].values)

    logger.info("Getting all genes from GFF database")
    return [gene.attributes["ID"][0] for gene in gff_db.all_features(featuretype="gene", order_by="start")]


#############################
# Rate calculation helpers  #
#############################


def get_alternates(seq: str, range_model: int) -> list[str]:
    """Return all alternative kmers with a different central nucleotide.

    Args:
        seq: kmer string.
        range_model: Number of flanking nucleotides on each side of the central base.

    Returns:
        List of alternate kmers (three per input kmer).
    """

    # Generate all single-nucleotide substitutions at the central position
    alternates = {seq[:range_model] + nuc + seq[range_model + 1 :] for nuc in ("A", "C", "T", "G")}
    # Remove the reference sequence
    alternates.discard(seq)
    return list(alternates)


def calculate_rates_cds(
    gene: gffutils.Feature,
    start: int,
    mutation_rate_model: pd.DataFrame,
    seq: str,
    range_model: int,
) -> list[pd.Series]:
    """Calculate per-nucleotide mutation rates for a CDS region using a kmer model.

    Args:
        gene: Gene feature from the GFF database.
        start: Genomic coordinate of the first base of ``seq``.
        mutation_rate_model: Mutation rate model indexed by ``<ref>_<alt>`` kmer keys.
        seq: Nucleotide sequence covering the CDS region (plus flanking context).
        range_model: Number of flanking nucleotides on each side of the central base.

    Returns:
        List of ``pd.Series`` objects, each representing one possible SNV.
    """

    # Reverse-complement for genes on the minus strand
    if gene.strand == "-":
        cur_seq = pyfaidx.complement(seq[::-1])
    else:
        cur_seq = seq

    len_seq = len(seq)
    list_mutations: list[pd.Series] = []

    for i in range(range_model, len(cur_seq) - range_model):
        ref = str(cur_seq[i - range_model : i + range_model + 1])
        list_alternates = get_alternates(ref, range_model)

        for alt in list_alternates:
            try:
                mutation_rate = mutation_rate_model.loc[ref + "_" + alt, "mu_snp"]
            except KeyError:
                # Ambiguous nucleotides — skip
                continue

            ref_nuc = ref[range_model]
            alt_nuc = alt[range_model]

            # Coordinates differ by strand orientation
            if gene.strand == "-":
                ref_nuc = pyfaidx.complement(ref_nuc)
                alt_nuc = pyfaidx.complement(alt_nuc)
                pos = start + len_seq - 1 - i
            else:
                pos = start + i

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
            list_mutations.append(s)

    return list_mutations


def get_sequence(fasta: pyfaidx.Fasta, chrom: str, start: int, end: int) -> str:
    """Return sequence from a FASTA file for the given genomic coordinates.

    Args:
        fasta: Genome sequence opened with pyfaidx.
        chrom: Chromosome identifier.
        start: Sequence start position (0-based).
        end: Sequence end position (exclusive).

    Returns:
        DNA sequence string.
    """

    return fasta[chrom][start:end]


#############################
# kmer model                #
#############################


def calculate_rates_kmer(
    mutation_rate_model: pd.DataFrame,
    fasta: pyfaidx.Fasta,
    gff_db: gffutils.FeatureDB,
    gene_list: list[str],
) -> pd.DataFrame:
    """Assign mutation rates using the kmer mutation rate model.

    Reference: https://www.nature.com/articles/ng.3050

    Args:
        mutation_rate_model: Mutation rate for every possible kmer transition.
        fasta: Genome sequence.
        gff_db: GFF database.
        gene_list: List of gene identifiers to process.

    Returns:
        DataFrame with columns ``gene_id``, ``chrom``, ``pos``, ``ref``, ``alt``, ``prob``.
    """

    logger = logging.getLogger("logger")

    # Default mutation model is 3-mer; other sizes (5-mer, 7-mer) derive the range automatically
    length_model = len(mutation_rate_model.iloc[0]["from"])
    range_model = length_model // 2

    nb_genes = len(gene_list)
    logger.info(f"Start computing rates for {nb_genes} genes")

    list_mutation_rates: list[pd.DataFrame] = []
    cpt = 1

    # Process one gene at a time
    for gene in gff_db.all_features(featuretype="gene"):

        # Skip genes not in the user-provided list (if given)
        gene_id = gene.attributes["ID"][0]
        if gene_id not in gene_list:
            continue

        list_mutation_rates_gene: list[pd.Series] = []
        seen_cds_boundaries: set[tuple[int, int]] = set()

        # Loop over all CDS regions
        for transcript in gff_db.children(gene, level=1):
            for cds in gff_db.children(transcript, featuretype="CDS", order_by="start"):

                # We also want to include possible mutations in the splice region, so we extend the CDS region by a fixed number of bases on each side (CDS_OFFSET)
                start = cds.start - CDS_OFFSET - range_model - 1
                end = cds.stop + CDS_OFFSET + range_model

                if (start, end) in seen_cds_boundaries:
                    continue
                seen_cds_boundaries.add((start, end))

                # Get the rates for all possible SNVs in the current CDS region
                cds_seq = get_sequence(fasta, gene.chrom, start, end)
                list_mutation_rates_gene += calculate_rates_cds(
                    gene, start + 1, mutation_rate_model, cds_seq, range_model
                )

        # Remove duplicates that may arise from overlapping CDS regions and sort by position and alternate allele
        if list_mutation_rates_gene:
            gene_df = pd.DataFrame(list_mutation_rates_gene)
            gene_df.drop_duplicates(inplace=True, keep="first")
            gene_df.sort_values(by=["pos", "alt"], inplace=True)
            list_mutation_rates.append(gene_df)
        else:
            logger.warning(f"No mutations found for gene {gene_id}")

        if cpt % 100 == 0:
            logger.info(f"{cpt}/{nb_genes} done")
        cpt += 1

    if not list_mutation_rates:
        logger.error("Rates dataframe is empty")
        return pd.DataFrame(columns=["gene_id", "chrom", "pos", "ref", "alt", "prob"])

    rates_df = pd.concat(list_mutation_rates)

    if (cpt - 1) != nb_genes:
        logger.warning(f"Rates computed for {cpt - 1} genes when number of genes in list is {nb_genes}")
        logger.warning(
            f"The following genes were not found in the GFF database : {set(gene_list) - set(rates_df['gene_id'])}"
        )
    else:
        logger.info(f"Rates computed for {nb_genes} genes")

    return rates_df


#############################
# Roulette model            #
#############################


def _roulette_per_chrom_files(roulette_dir: Path) -> dict[str, Path]:
    """Build a mapping from chromosome name to the corresponding Roulette VCF path.

    Args:
        roulette_dir: Directory containing per-chromosome Roulette VCF files.

    Returns:
        Dict keyed by chromosome name (without ``chr`` prefix) mapping to the file path.
    """

    result: dict[str, Path] = {}
    for f in roulette_dir.glob("*all.vcf*gz"):
        chrom = f.name.split("_")[0]
        result[chrom] = f
    return result


def calculate_rates_roulette(
    roulette_dir: Path,
    gff_db: gffutils.FeatureDB,
    gene_list: list[str],
    model: str,
    scaling_factor: float,
) -> pd.DataFrame:
    """Assign mutation rates using the Roulette or Carlson mutation rate model.

    References:
        - https://github.com/vseplyarskiy/Roulette
        - https://www.biorxiv.org/content/10.1101/2022.08.20.504670v1

    Args:
        roulette_dir: Directory containing per-chromosome Roulette VCF files.
            Each VCF is annotated with ``MR`` (Roulette) and ``MC`` (Carlson) INFO fields.
        gff_db: GFF database.
        gene_list: List of gene identifiers to process.
        model: Either ``"roulette"`` or ``"carlson"``.
        scaling_factor: Per-generation mutation rate scaling factor applied to every
            rate value extracted from the VCF INFO field (``MR`` for Roulette, ``MC`` for Carlson).

    Returns:
        DataFrame with columns ``gene_id``, ``chrom``, ``pos``, ``ref``, ``alt``, ``prob``.
    """

    logger = logging.getLogger("logger")

    if model == "roulette":
        info_field = "MR"
    else:
        info_field = "MC"

    nb_genes = len(gene_list)
    logger.info(f"Start computing rates for {nb_genes} genes")

    roulette_vcfs = _roulette_per_chrom_files(roulette_dir)

    list_mutation_rates: list[pd.DataFrame] = []
    cpt = 1

    for gene in gff_db.all_features(featuretype="gene"):
        gene_id = gene.attributes["ID"][0]
        if gene_id not in gene_list:
            continue

        gff_uses_chr_prefix = gene.chrom.startswith("chr")
        chrom = gene.chrom.replace("chr", "")

        try:
            roulette_index = None
            for ext in (".csi", ".tbi"):
                candidate_index = Path(str(roulette_vcfs[chrom]) + ext)
                if candidate_index.exists():
                    roulette_index = str(candidate_index)
                    break
            if roulette_index is None:
                raise FileNotFoundError(f"No index file found for {roulette_vcfs[chrom]}")
            roulette_file = pysam.VariantFile(str(roulette_vcfs[chrom]), index_filename=roulette_index)
        except KeyError:
            # Roulette does not provide mutation rates for allosomes
            if chrom not in ("X", "Y"):
                logger.warning(f"Can't find any roulette file for chromosome {chrom} (gene {gene_id})")
            continue

        list_mutation_rates_gene: list[pd.Series] = []
        seen_cds_boundaries: set[tuple[int, int]] = set()

        for transcript in gff_db.children(gene, level=1):
            for cds in gff_db.children(transcript, featuretype="CDS", order_by="start"):
                start = cds.start - CDS_OFFSET - 1
                end = cds.stop + CDS_OFFSET

                if (start, end) in seen_cds_boundaries:
                    continue
                seen_cds_boundaries.add((start, end))

                for rec in roulette_file.fetch(start=start, stop=end, region=chrom):
                    try:
                        s = pd.Series(
                            {
                                "gene_id": gene_id,
                                "chrom": f"chr{chrom}" if gff_uses_chr_prefix else chrom,
                                "pos": rec.pos,
                                "ref": rec.ref,
                                "alt": rec.alts[0],
                                "prob": float(rec.info[info_field]) * scaling_factor,
                            }
                        )
                        list_mutation_rates_gene.append(s)
                    except KeyError:
                        # No mutation rate info for this SNP
                        continue

        if list_mutation_rates_gene:
            gene_df = pd.DataFrame(list_mutation_rates_gene)
            gene_df.drop_duplicates(inplace=True, keep="first")
            gene_df.sort_values(by=["pos", "alt"], inplace=True)
            list_mutation_rates.append(gene_df)
        else:
            logger.warning(f"No mutations found for gene {gene_id}")

        if cpt % 100 == 0:
            logger.info(f"{cpt}/{nb_genes} done")
        cpt += 1

    if not list_mutation_rates:
        logger.error("Rates dataframe is empty")
        return pd.DataFrame(columns=["gene_id", "chrom", "pos", "ref", "alt", "prob"])

    rates_df = pd.concat(list_mutation_rates)

    if (cpt - 1) != nb_genes:
        logger.warning(f"Rates computed for {cpt - 1} genes when number of genes in list is {nb_genes}")
        logger.warning(
            f"The following genes were not found in the GFF database : {set(gene_list) - set(rates_df['gene_id'])}"
        )
    else:
        logger.info(f"Rates computed for {nb_genes} genes")

    return rates_df


#############################
# CLI                       #
#############################


def _validate_params(model: str, fasta: Path | None, scaling_factor: float | None) -> None:
    """Validate cross-parameter dependencies that cannot be expressed as individual option constraints.

    Args:
        model: Selected mutation rate model.
        fasta: Path to the genome FASTA file (required for the kmer model).
        scaling_factor: Scaling factor (required for carlson/roulette models).

    Raises:
        click.UsageError: When a required companion parameter is missing.
    """

    if model in ("carlson", "roulette") and scaling_factor is None:
        raise click.UsageError(f'--model "{model}" requires --scaling_factor to be provided')
    if model == "kmer" and fasta is None:
        raise click.UsageError('--model "kmer" requires --fasta to be provided')


def _setup_conf(ctx: click.Context) -> None:
    """Initialise logging, log all parameters, and create the output directory.

    Args:
        ctx: Click context holding the fully-resolved command parameters.
    """

    init_log()
    logger = logging.getLogger("logger")
    logger.info("Running {}".format(Path(__file__).name))

    logger.info("Parameters :")
    logger.info("----------")
    for key, value in ctx.params.items():
        if value is not None:
            logger.info(f"{key} : {value}")
    logger.info("----------")

    Path(ctx.params["outdir"]).mkdir(parents=True, exist_ok=True)


@click.command()
@click.option(
    "--rates_model_path",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Path to a k-mer mutation rate model or roulette directory",
)
@click.option(
    "--model",
    required=True,
    type=click.Choice(["kmer", "carlson", "roulette"]),
    help="Mutation rate model to use.",
)
@click.option(
    "--gff",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="GFF file with gene annotations or an existing gffutils .db file",
)
@click.option("--outdir", required=True, type=click.Path(path_type=Path), help="Output directory")
@click.option(
    "--gene_list",
    type=click.Path(exists=True, path_type=Path),
    help="Optional file containing gene identifiers to process (one per line)",
)
@click.option(
    "--fasta",
    type=click.Path(exists=True, path_type=Path),
    help="Genome sequence in FASTA format (required for the kmer model)",
)
@click.option(
    "--scaling_factor",
    type=click.FLOAT,
    help="Scaling factor to apply to mutation rates when using the Roulette or Carlson model (ignored for kmer model)",
)
@click.option(
    "--output-rates",
    type=str,
    default="mutation_rates.tsv",
    help="Name of the output mutation rates file",
)
def main(
    rates_model_path: Path,
    model: str,
    gff: Path,
    outdir: Path,
    gene_list: Path,
    fasta: Path,
    scaling_factor: float,
    output_rates: str,
) -> None:
    """Generate a per-site per-generation haploid mutation rates file.

    The output file (``mutation_rates.tsv``) contains one row per possible SNV
    across all CDS regions (extended by ``CDS_OFFSET`` bases on each side to
    capture splicing-relevant positions) of the requested genes. Each row holds
    the gene identifier, genomic coordinates, reference and alternate alleles,
    and the estimated per-generation mutation probability (``prob``).

    This file feeds directly into the annotation and enrichment-test steps of
    the pipeline: annotations are added column-wise, and the ``prob`` column is
    used by the simulation to weight randomly drawn mutations.

    Three mutation rate models are supported:

    * **kmer** – trinucleotide (or higher-order) model; requires a FASTA genome
      and a tab-separated model file with columns ``from``, ``to``, ``mu_snp``.
    * **carlson** / **roulette** – pre-computed rates from the Roulette resource;
      require a directory of per-chromosome indexed VCF files.

    Args:
        rates_model_path: Path to a k-mer model TSV or the Roulette VCF directory.
        model: Mutation rate model.  One of ``kmer``, ``carlson``, ``roulette``.
        gff: Gene annotations as a GFF file or an existing gffutils ``.db`` file.
        outdir: Output directory; created automatically if absent.
        gene_list: File containing gene identifiers to process (one per line).
            When omitted, all genes in the GFF are processed.
        fasta: Genome sequence in FASTA format (required for the kmer model).
        scaling_factor: Per-generation scaling factor applied to Roulette/Carlson rates.
        output_rates: Name of the output mutation rates file.
    """

    _setup_conf(click.get_current_context())
    _validate_params(model=model, fasta=fasta, scaling_factor=scaling_factor)

    # Load GFF file or database
    gff_db = load_gff(gff, outdir / "gff.db")

    # Load gene list if provided, otherwise get all genes in GFF
    genes = load_gene_list(gene_list, gff_db)

    if model == "kmer":
        rate_model = load_kmer_mutation_rate_model(rates_model_path)
        genome = pyfaidx.Fasta(fasta)
        mutation_rates_df = calculate_rates_kmer(rate_model, genome, gff_db, genes)
    else:
        mutation_rates_df = calculate_rates_roulette(rates_model_path, gff_db, genes, model, scaling_factor)

    # Export mutation rates file
    out_path = outdir / output_rates
    mutation_rates_df.to_csv(out_path, sep="\t", index=False)

    logging.getLogger("logger").info(f"Mutation rates file created here : {out_path}")


if __name__ == "__main__":
    main()
