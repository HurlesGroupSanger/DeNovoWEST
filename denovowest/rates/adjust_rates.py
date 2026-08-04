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


def _update_gene_mutation_rates_before_after(
    gene_mutation_rates_before_after: dict[str, dict[str, float]],
    chunk: pd.DataFrame,
    chunk_adjusted: pd.DataFrame,
    adjusted_genes: list[str],
) -> dict[str, dict[str, float]]:
    """Update the gene mutation rates before/after adjustment.

    Args:
        gene_mutation_rates_before_after: Dictionary to update with mutation rates.
        chunk: DataFrame chunk containing the mutation rates.
        stage: Either "before" or "after" to indicate which stage of adjustment.

    Returns:
        Updated gene_mutation_rates_before_after dictionary.
    """
    for gene_id in adjusted_genes:
        total_prob_before = chunk.loc[chunk.gene_id == gene_id, "prob"].sum()
        total_prob_after = chunk_adjusted.loc[chunk_adjusted.gene_id == gene_id, "prob"].sum()

        entry = gene_mutation_rates_before_after.setdefault(gene_id, {})
        entry["before"] = entry.get("before", 0.0) + total_prob_before
        entry["after"] = entry.get("after", 0.0) + total_prob_after

    return gene_mutation_rates_before_after


def _export_adjustment_stats(gene_mutation_rates_before_after: dict, path: Path) -> None:
    """Export per-gene mutation rates before and after adjustment to a JSON file.

    Only genes whose total mutation rate changed are included in the output.

    Args:
        gene_mutation_rates_before_after: Mapping of gene_id to {"before": float, "after": float}.
        path: Destination path for the JSON file.
    """

    adjusted_genes = []
    for gene_id, rates in gene_mutation_rates_before_after.items():
        before = rates.get("before", 0.0)
        after = rates.get("after", 0.0)

        if before != after:
            adjusted_genes.append(gene_id)

    gene_mutation_rates_before_after = {
        gene_id: rates for gene_id, rates in gene_mutation_rates_before_after.items() if gene_id in adjusted_genes
    }

    with open(path, "w") as f:
        json.dump(gene_mutation_rates_before_after, f, indent=2)


##############################
# Sperm selection adjustment #
##############################


def apply_sperm_selection_adjustment(
    chunk: pd.DataFrame,
    adjustment_genes: dict[str, dict[str, float]],
    adjustment_sites: dict[str, dict[str, Any]],
    mapping_gene_ids: dict[str, str],
) -> pd.DataFrame:
    """Apply per-gene-per-consequence and per-site multipliers to a DataFrame chunk.

    Args:
        chunk: A slice of the mutation rates DataFrame with columns
            gene_id, consequence, prob, chrom, pos, ref, alt.
        adjustment_genes: Mapping of gene_id to {consequence: multiplier}.
        adjustment_sites: Mapping of site_id to {chrom, pos, ref, alt, multiplier}.
        mapping_gene_ids: ENSEMBL version mapping from adjustment gene IDs to rates
            file gene IDs, computed per chunk via ensembl_ids_map_version.

    Returns:
        Adjusted chunk.
    """

    adjusted_genes = list()  # Track which genes were adjusted in this chunk

    # Apply per-gene adjustments
    for gene_id, gene_multipliers in adjustment_genes.items():
        if gene_id not in mapping_gene_ids:
            continue
        mapped_id = mapping_gene_ids[gene_id]
        adjusted_genes.append(mapped_id)
        for consequence, multiplier in gene_multipliers.items():
            mask = (chunk["gene_id"] == mapped_id) & (chunk["consequence"] == consequence)
            chunk.loc[mask, "prob"] *= multiplier

    # Apply per-site adjustments, tracking which site IDs are matched in this chunk
    for site_id, site_info in adjustment_sites.items():
        mask = (
            (chunk["chrom"] == site_info["chrom"])
            & (chunk["pos"] == site_info["pos"])
            & (chunk["ref"] == site_info["ref"])
            & (chunk["alt"] == site_info["alt"])
        )
        if mask.any():
            chunk.loc[mask, "prob"] *= site_info["multiplier"]
            adjusted_genes += list(chunk.loc[mask, "gene_id"].unique())

    return chunk, list(set(adjusted_genes))


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
    "--output",
    type=str,
    default="mutation_rates_adjusted.tsv",
    help="Name of the output mutation rates file",
)
@click.option(
    "--chunk-size",
    type=int,
    default=1_000_000,
    show_default=True,
    help="Number of rows to process at a time",
)
def main(
    rates: Path,
    sperm_selection_adjustment: Path,
    outdir: Path,
    output: str,
    chunk_size: int,
) -> None:
    """Adjust mutation rates by applying per-gene and per-site multipliers.

    Reads a mutation rates TSV, applies the multipliers from a JSON adjustment file,
    and writes the adjusted rates to the output directory. The file is processed in
    chunks to keep memory use bounded regardless of file size.
    """

    init_log()
    log = logging.getLogger("logger")
    log.info(f"Adjusting mutation rates using sperm selection adjustment file: {sperm_selection_adjustment}")

    outdir.mkdir(parents=True, exist_ok=True)

    # Load the adjustment JSON file and validate its structure
    adjustment_genes, adjustment_sites = _load_adjustment_json(sperm_selection_adjustment)

    nb_genes_seen = 0
    out_path = outdir / output
    gene_mutation_rates_before_after = dict()
    # Read the mutation rates file in chunks, apply adjustments, and write to output
    for i, chunk in enumerate(
        pd.read_csv(
            rates,
            sep="\t",
            chunksize=chunk_size,
            low_memory=False,
            dtype={"chrom": str, "pos": int, "ref": str, "alt": str, "gene_id": str},
        )
    ):

        log.info("Processing chunk %d (rows %d to %d)", i + 1, i * chunk_size + 1, (i + 1) * chunk_size)

        # In case the gene IDs in the adjustment file are from a different ENSEMBL version than the rates file (or do not use the version)
        mapping_gene_ids = ensembl_ids_map_version(list(adjustment_genes.keys()), chunk["gene_id"].unique().tolist())

        # Adjustment
        chunk_adjusted, adjusted_genes = apply_sperm_selection_adjustment(
            chunk.copy(), adjustment_genes, adjustment_sites, mapping_gene_ids
        )

        # Write
        chunk_adjusted.to_csv(out_path, sep="\t", index=False, mode="w" if i == 0 else "a", header=(i == 0))

        # Log adjusted genes
        if adjusted_genes:
            gene_mutation_rates_before_after = _update_gene_mutation_rates_before_after(
                gene_mutation_rates_before_after, chunk, chunk_adjusted, adjusted_genes
            )

    log.info(f"Adjusted mutation rates file created here: {out_path}")

    _export_adjustment_stats(gene_mutation_rates_before_after, outdir / "adjusted_rates.json")


if __name__ == "__main__":
    main()
