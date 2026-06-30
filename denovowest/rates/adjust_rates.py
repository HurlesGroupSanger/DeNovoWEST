#!/usr/bin/env python

import json
import logging
from pathlib import Path
from typing import Any

import click
import pandas as pd

from denovowest.utils.helpers import ensembl_ids_map_version
from denovowest.utils.log import init_log

#########
# Utils #
#########


def _load_adjustment_json(
    path: Path,
) -> tuple[dict[str, dict[str, float]], dict[str, dict[str, Any]]]:
    """Load and validate a sperm selection adjustment JSON file.

    Args:
        path: Path to the JSON file containing at least one of the 'genes' or
            'sites' top-level keys.

    Returns:
        Tuple of (adjustment_genes, adjustment_sites) dicts.

    Raises:
        click.UsageError: If the file is missing required top-level keys or any site
            entry is missing required coordinate or multiplier fields.
    """

    with open(path, "r") as f:
        d = json.load(f)

    # At least one of 'genes' or 'sites' must be present
    if "genes" not in d and "sites" not in d:
        raise click.UsageError("Adjustment JSON must contain at least one of the keys: 'genes', 'sites'")

    # Validate required fields in each site entry
    required_site_fields = {"chrom", "pos", "ref", "alt", "multiplier"}
    for site_id, site_info in d.get("sites", {}).items():
        missing_fields = required_site_fields - set(site_info.keys())
        if missing_fields:
            raise click.UsageError(
                f"Site '{site_id}' in adjustment JSON is missing required fields: "
                f"{', '.join(sorted(missing_fields))}"
            )

    return d.get("genes", {}), d.get("sites", {})


def _build_gene_version_mapping(
    rates: Path,
    adjustment_genes: dict[str, dict[str, float]],
    chunk_size: int,
) -> dict[str, str]:
    """Scan the rates file to build an ENSEMBL version mapping for the genes.

    Reads only the gene_id column to minimise memory use, then delegates to
    ensembl_ids_map_version to resolve versioned/unversioned identifier mismatches.
    Logs a warning for any gene that has no match in the rates file.

    Args:
        rates: Path to the mutation rates file.
        adjustment_genes: Gene adjustment dict whose keys are the source gene IDs.
        chunk_size: Number of rows to read per chunk during the scan.

    Returns:
        Dict mapping each source gene ID to the best-matching gene ID in the rates file.
        Gene IDs with no match are absent from the result.
    """

    log = logging.getLogger("logger")

    if not adjustment_genes:
        return {}

    # Collect all unique gene IDs from the rates file, reading only the gene_id column
    gene_ids: set[str] = set()
    for chunk in pd.read_csv(rates, sep="\t", usecols=["gene_id"], chunksize=chunk_size):
        gene_ids.update(chunk["gene_id"].unique().tolist())

    mapping = ensembl_ids_map_version(list(adjustment_genes.keys()), list(gene_ids))

    return mapping


##############################
# Sperm selection adjustment #
##############################


def apply_sperm_selection_adjustment(
    chunk: pd.DataFrame,
    adjustment_genes: dict[str, dict[str, float]],
    adjustment_sites: dict[str, dict[str, Any]],
    mapping_gene_ids: dict[str, str],
) -> tuple[pd.DataFrame, set[str]]:
    """Apply per-gene-per-consequence and per-site multipliers to a DataFrame chunk.

    Args:
        chunk: A slice of the mutation rates DataFrame with columns
            gene_id, consequence, prob, chrom, pos, ref, alt.
        adjustment_genes: Mapping of gene_id to {consequence: multiplier}.
        adjustment_sites: Mapping of site_id to {chrom, pos, ref, alt, multiplier}.
        mapping_gene_ids: Pre-built ENSEMBL version mapping from adjustment gene IDs
            to rates file gene IDs, as returned by _build_gene_version_mapping.

    Returns:
        Tuple of (adjusted chunk, set of site IDs matched in this chunk).
    """

    # Apply per-gene adjustments
    for gene_id, gene_multipliers in adjustment_genes.items():
        if gene_id not in mapping_gene_ids:
            continue
        mapped_id = mapping_gene_ids[gene_id]
        for consequence, multiplier in gene_multipliers.items():
            mask = (chunk["gene_id"] == mapped_id) & (chunk["consequence"] == consequence)
            chunk.loc[mask, "prob"] *= multiplier

    # Apply per-site adjustments, tracking which site IDs are matched in this chunk
    matched_site_ids: set[str] = set()
    for site_id, site_info in adjustment_sites.items():
        mask = (
            (chunk["chrom"] == site_info["chrom"])
            & (chunk["pos"] == site_info["pos"])
            & (chunk["ref"] == site_info["ref"])
            & (chunk["alt"] == site_info["alt"])
        )
        if mask.any():
            chunk.loc[mask, "prob"] *= site_info["multiplier"]
            matched_site_ids.add(site_id)

    return chunk, matched_site_ids


##############
# Entry point#
##############


@click.command()
@click.option(
    "--rates",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Path to a mutation rates file (plain TSV or bgzipped)",
)
@click.option(
    "--sperm-selection-adjustment",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="JSON file with gene-based or site-based multipliers for spermatogonial selection",
)
@click.option("--outdir", default=Path("."), type=click.Path(path_type=Path), help="Output directory")
@click.option(
    "--output-rates",
    type=str,
    default="mutation_rates_adjusted.tsv",
    help="Name of the output mutation rates file",
)
@click.option(
    "--chunk-size",
    type=int,
    default=100_000,
    show_default=True,
    help="Number of rows to process at a time",
)
def main(
    rates: Path,
    sperm_selection_adjustment: Path,
    outdir: Path,
    output_rates: str,
    chunk_size: int,
) -> None:
    """Adjust mutation rates by applying per-gene and per-site multipliers.

    Reads a mutation rates TSV, applies the multipliers from a JSON adjustment file,
    and writes the adjusted rates to the output directory. The file is processed in
    chunks to keep memory use bounded regardless of file size.
    """

    init_log()
    log = logging.getLogger("logger")

    outdir.mkdir(parents=True, exist_ok=True)

    adjustment_genes, adjustment_sites = _load_adjustment_json(sperm_selection_adjustment)

    # Preliminary scan: read only the gene_id column to build the ENSEMBL version mapping
    mapping_gene_ids = _build_gene_version_mapping(rates, adjustment_genes, chunk_size)

    # Stream through the rates file, applying adjustments and writing output incrementally
    out_path = outdir / output_rates
    all_matched_site_ids: set[str] = set()

    for i, chunk in enumerate(pd.read_csv(rates, sep="\t", chunksize=chunk_size)):
        chunk, matched = apply_sperm_selection_adjustment(chunk, adjustment_genes, adjustment_sites, mapping_gene_ids)
        all_matched_site_ids.update(matched)
        chunk.to_csv(out_path, sep="\t", index=False, mode="w" if i == 0 else "a", header=(i == 0))

    log.info(f"Adjusted mutation rates file created here: {out_path}")


if __name__ == "__main__":
    main()
