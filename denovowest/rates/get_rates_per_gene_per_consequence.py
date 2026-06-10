#!/usr/bin/env python

from collections import defaultdict
from pathlib import Path
from typing import Generator

import click
import pandas as pd

from denovowest.utils.params import ConsequenceGroups

##############
# CONSTANTS  #
##############

CONSEQUENCES = ConsequenceGroups.OE_RATIO

###########
# HELPERS #
###########


def _iter_rates(rates: Path) -> Generator[tuple[str, pd.DataFrame], None, None]:
    """Yield (gene_id, sub-dataframe) pairs from a rates file, chunk by chunk.

    Note:
        The same gene_id may be yielded more than once when a gene straddles a
        chunk boundary; callers must accumulate across repeated occurrences.

    Args:
        rates: Path to a tab-delimited rates file with at least ``gene_id``,
            ``consequence``, and ``prob`` columns.

    Yields:
        Tuples of (gene_id, rows for that gene in this chunk).
    """
    reader = pd.read_csv(rates, sep="\t", chunksize=100_000, low_memory=False)
    for chunk in reader:
        for gene_id, group in chunk.groupby("gene_id", sort=False):
            yield str(gene_id), group


###########
# COMPUTE #
###########


def compute_rates_per_gene_per_consequence(
    rates: Path,
) -> tuple[dict[str, dict[str, float]], dict[str, dict[str, int]]]:
    """Accumulate mutation probabilities and record counts per gene per consequence.

    Args:
        rates: Path to the rates TSV file.

    Returns:
        A pair ``(prob, count)`` where:

        - ``prob``: gene_id → consequence → sum of ``prob`` values.
        - ``count``: gene_id → consequence → number of records.
    """
    prob: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))
    count: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))

    n_chunks = 0
    for gene_id, df in _iter_rates(rates):
        for csq in CONSEQUENCES:
            mask = df["consequence"] == csq
            prob[gene_id][csq] += float(df.loc[mask, "prob"].sum())
            count[gene_id][csq] += int(mask.sum())

        n_chunks += 1
        if n_chunks % 100 == 0:
            print(f"Processed {n_chunks} gene chunks")

    return dict(prob), dict(count)


##########
# OUTPUT #
##########


def _write_matrix(data: dict[str, dict[str, float]], output: Path) -> None:
    pd.DataFrame(data).transpose().to_csv(output, sep="\t")


#######
# CLI #
#######


@click.command()
@click.argument("rates", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--output",
    type=click.Path(path_type=Path),
    default="rates_per_gene_per_consequence.tsv",
    show_default=True,
)
@click.option(
    "--output-nb",
    type=click.Path(path_type=Path),
    default="nb_records_per_gene_per_consequence.tsv",
    show_default=True,
)
def main(rates: Path, output: Path, output_nb: Path) -> None:
    """Compute expected number of mutations per gene per consequence type.

    \b
    Arguments:
        RATES: Tab-delimited file with at least gene_id, consequence, and prob columns.
    """
    prob, count = compute_rates_per_gene_per_consequence(rates)
    _write_matrix(prob, output)
    _write_matrix(count, output_nb)


if __name__ == "__main__":
    main()
