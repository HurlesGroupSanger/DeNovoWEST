#!/usr/bin/env python
"""
Generate a self-contained HTML report from DeNovoWEST pipeline results.

Reads combined_tests_results.tsv and any available per-gene stats files from
the results directory, then writes a single HTML file with a sortable gene
table and expandable per-gene detail panels.
"""

####################
# Imports          #
####################

import json
import logging
from datetime import date
from pathlib import Path
from typing import Optional

import click
import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape

from denovowest.utils.log import init_log

####################
# Constants        #
####################

COMBINED_TESTS_RELPATH = "combined_tests/combined_tests_results.tsv"
NB_RECORDS_RELPATH = "rates/stats/nb_records_per_gene_per_csq.tsv"
RATES_RELPATH = "rates/stats/rates_per_gene_per_csq.tsv"

PVAL_COLUMNS = [
    "min_pval",
    "combined_mis_pval",
    "enrichment_all_pval",
    "enrichment_mis_pval",
    "clustering_pval",
]

FUNCTIONAL_CONSEQUENCES = [
    "missense",
    "stop_gained",
    "start_lost",
    "stop_lost",
    "splice_region",
    "splice_acceptor",
    "splice_donor",
]

ALL_CONSEQUENCES = ["synonymous"] + FUNCTIONAL_CONSEQUENCES

# Fixed nominal threshold; the corrected threshold is supplied by the user at runtime.
NOMINAL_THRESHOLD = 5e-2

####################
# Data loading     #
####################


def _load_wide_stats(path: Path) -> Optional[pd.DataFrame]:
    """Load a wide-format stats TSV (genes as rows, consequences as columns)."""
    if not path.exists():
        return None
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index.name = "gene_id"
    return df


def load_combined_tests(results_dir: Path) -> pd.DataFrame:
    """Load the combined tests results TSV.

    Args:
        results_dir: pipeline output directory.

    Returns:
        DataFrame with one row per gene.

    Raises:
        FileNotFoundError: if the expected TSV is absent.
    """

    path = results_dir / COMBINED_TESTS_RELPATH
    if not path.exists():
        raise FileNotFoundError(f"Combined tests results not found: {path}")
    return pd.read_csv(path, sep="\t")


def load_nb_records(results_dir: Path) -> Optional[pd.DataFrame]:
    """Load the per-gene per-consequence record count matrix, if present."""
    return _load_wide_stats(results_dir / NB_RECORDS_RELPATH)


def load_rates(results_dir: Path) -> Optional[pd.DataFrame]:
    """Load the per-gene per-consequence mutation rate matrix, if present."""
    return _load_wide_stats(results_dir / RATES_RELPATH)


####################
# Caution flags    #
####################


def _check_missing_pvals(row: pd.Series) -> list[str]:
    """Flag genes where one or more p-value tests produced no result."""
    missing = [col for col in PVAL_COLUMNS if pd.isna(row.get(col))]
    if missing:
        return [f"Missing p-value(s): {', '.join(missing)}"]
    return []


def _check_zero_rates(gene_id: str, rates_df: Optional[pd.DataFrame]) -> list[str]:
    """Flag genes with a zero total mutation rate for any functional consequence."""
    if rates_df is None or gene_id not in rates_df.index:
        return []
    row = rates_df.loc[gene_id]
    zero_csqs = [
        csq for csq in FUNCTIONAL_CONSEQUENCES
        if csq in row.index and row[csq] == 0.0
    ]
    if zero_csqs:
        return [f"Zero mutation rate for: {', '.join(zero_csqs)}"]
    return []


def _check_zero_records(gene_id: str, nb_records_df: Optional[pd.DataFrame]) -> list[str]:
    """Flag genes with no annotated records for any functional consequence."""
    if nb_records_df is None or gene_id not in nb_records_df.index:
        return []
    row = nb_records_df.loc[gene_id]
    zero_csqs = [
        csq for csq in FUNCTIONAL_CONSEQUENCES
        if csq in row.index and row[csq] == 0
    ]
    if zero_csqs:
        return [f"No annotated records for: {', '.join(zero_csqs)}"]
    return []


def compute_caution_flags(
    df: pd.DataFrame,
    nb_records_df: Optional[pd.DataFrame],
    rates_df: Optional[pd.DataFrame],
) -> pd.DataFrame:
    """Add caution flag and reasons columns to the results dataframe.

    Args:
        df: combined tests results.
        nb_records_df: per-gene per-consequence record counts (optional).
        rates_df: per-gene per-consequence mutation rates (optional).

    Returns:
        Copy of df with added 'caution' (bool) and 'caution_reasons' (list[str]) columns.
    """

    reasons_list = []
    for _, row in df.iterrows():
        gene_id = row["gene_id"]
        reasons = (
            _check_missing_pvals(row)
            + _check_zero_rates(gene_id, rates_df)
            + _check_zero_records(gene_id, nb_records_df)
        )
        reasons_list.append(reasons)

    df = df.copy()
    df["caution_reasons"] = reasons_list
    df["caution"] = df["caution_reasons"].map(lambda r: len(r) > 0)
    return df


####################
# Data serialisation
####################


def _safe_float(val) -> Optional[float]:
    """Return float or None for NaN/missing values."""
    try:
        if pd.isna(val):
            return None
    except (TypeError, ValueError):
        pass
    return float(val)


def _safe_int(val) -> Optional[int]:
    """Return int or None for NaN/missing values."""
    try:
        if pd.isna(val):
            return None
    except (TypeError, ValueError):
        pass
    return int(val)


def _build_gene_records(
    df: pd.DataFrame,
    nb_records_df: Optional[pd.DataFrame],
    rates_df: Optional[pd.DataFrame],
) -> list[dict]:
    """Serialise the results dataframe into a JSON-friendly list of dicts.

    Args:
        df: combined tests results with caution columns.
        nb_records_df: optional per-consequence record counts.
        rates_df: optional per-consequence mutation rates.

    Returns:
        List of dicts, one per gene, ready for json.dumps().
    """

    records = []
    for _, row in df.iterrows():
        gene_id = row["gene_id"]
        base_id = gene_id.split(".")[0]

        # Build per-consequence breakdown for detail panel
        csq_data = {}
        for csq in ALL_CONSEQUENCES:
            entry = {}
            if (
                nb_records_df is not None
                and gene_id in nb_records_df.index
                and csq in nb_records_df.columns
            ):
                v = nb_records_df.loc[gene_id, csq]
                entry["nb_records"] = _safe_int(v)
            if (
                rates_df is not None
                and gene_id in rates_df.index
                and csq in rates_df.columns
            ):
                v = rates_df.loc[gene_id, csq]
                entry["rate"] = _safe_float(v)
            if entry:
                csq_data[csq] = entry

        records.append(
            {
                "gene_id": gene_id,
                "base_id": base_id,
                "min_pval": _safe_float(row.get("min_pval")),
                "combined_mis_pval": _safe_float(row.get("combined_mis_pval")),
                "enrichment_all_pval": _safe_float(row.get("enrichment_all_pval")),
                "enrichment_mis_pval": _safe_float(row.get("enrichment_mis_pval")),
                "clustering_pval": _safe_float(row.get("clustering_pval")),
                "observed_all": _safe_float(row.get("observed_all")),
                "expected_all": _safe_float(row.get("expected_all")),
                "observed_mis": _safe_float(row.get("observed_mis")),
                "expected_mis": _safe_float(row.get("expected_mis")),
                "nb_mis_variants": _safe_int(row.get("nb_mis_variants")),
                "clustering_dist": _safe_float(row.get("clustering_dist")),
                "caution": bool(row["caution"]),
                "caution_reasons": row["caution_reasons"],
                "csq_data": csq_data,
            }
        )
    return records


####################
# HTML generation  #
####################

_TEMPLATE_DIR = Path(__file__).parent / "templates"


def build_html(
    df: pd.DataFrame,
    nb_records_df: Optional[pd.DataFrame],
    rates_df: Optional[pd.DataFrame],
    results_dir: Path,
    run_date: str,
    corrected_threshold: float,
) -> str:
    """Render the complete HTML report via the Jinja2 template.

    Args:
        df: combined tests results with caution columns.
        nb_records_df: optional per-consequence record counts.
        rates_df: optional per-consequence mutation rates.
        results_dir: pipeline output directory (used in the report title).
        run_date: ISO-format date string for the report header.
        corrected_threshold: multiple-testing corrected significance threshold.

    Returns:
        Self-contained HTML document as a string.

    Note:
        gene_json is passed as a plain string and rendered with the ``| safe``
        filter in the template, so autoescape does not double-encode it.
    """

    n_total = len(df)
    n_significant = int((df["min_pval"].fillna(1.0) < corrected_threshold).sum())
    n_caution = int(df["caution"].sum())

    gene_records = _build_gene_records(df, nb_records_df, rates_df)
    # Escape </ to prevent a premature </script> tag inside the JSON block
    gene_json = json.dumps(gene_records, indent=2).replace("</", "<\\/")

    env = Environment(
        loader=FileSystemLoader(str(_TEMPLATE_DIR)),
        autoescape=select_autoescape(["html"]),
    )
    template = env.get_template("report.html.j2")
    return template.render(
        title=f"DeNovoWEST Report — {results_dir.name}",
        results_dir=str(results_dir),
        run_date=run_date,
        n_total=n_total,
        n_significant=n_significant,
        n_caution=n_caution,
        corrected_threshold=corrected_threshold,
        nominal_threshold=NOMINAL_THRESHOLD,
        gene_json=gene_json,
    )


####################
# CLI entry point  #
####################


@click.command()
@click.argument(
    "results_dir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option(
    "--threshold",
    "-t",
    type=float,
    required=True,
    help="Multiple-testing corrected significance threshold (e.g. 2.5e-6 for Bonferroni).",
)
@click.option(
    "--output",
    "-o",
    type=click.Path(path_type=Path),
    default=None,
    help="Output HTML file. Defaults to <results_dir>/report.html.",
)
def report(results_dir: Path, threshold: float, output: Optional[Path]) -> None:
    """Generate an HTML report summarising DeNovoWEST pipeline results.

    RESULTS_DIR is the root output directory produced by the pipeline
    (e.g. /path/to/benchmarking_output/2026-06-10).

    Args:
        results_dir: pipeline output root.
        threshold: multiple-testing corrected significance threshold.
        output: destination HTML file (default: <results_dir>/report.html).
    """

    init_log()
    logger = logging.getLogger("logger")

    if output is None:
        output = results_dir / "report.html"

    run_date = date.today().isoformat()
    logger.info(f"Loading results from {results_dir}")
    logger.info(f"Corrected significance threshold: {threshold:.2e}")

    # Load all available data sources
    df = load_combined_tests(results_dir)
    nb_records_df = load_nb_records(results_dir)
    rates_df = load_rates(results_dir)

    if nb_records_df is not None:
        logger.info(f"Loaded record counts for {len(nb_records_df)} genes")
    else:
        logger.warning("Record count stats not found — caution flags will be limited")

    if rates_df is not None:
        logger.info(f"Loaded mutation rates for {len(rates_df)} genes")
    else:
        logger.warning("Mutation rate stats not found — caution flags will be limited")

    # Compute caution flags and generate HTML
    df = compute_caution_flags(df, nb_records_df, rates_df)
    n_caution = int(df["caution"].sum())
    logger.info(f"{n_caution} gene(s) flagged with caution")

    html = build_html(df, nb_records_df, rates_df, results_dir, run_date, threshold)
    output.write_text(html, encoding="utf-8")
    logger.info(f"Report written to {output}")


if __name__ == "__main__":
    report()
