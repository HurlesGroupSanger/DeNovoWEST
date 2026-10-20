# HTML report

`denovowest report` summarises the results of a workflow run in a single HTML page: a
sortable table of all tested genes, with details and quality warnings for each one.

```bash
denovowest report /analyses/out/full -t 1.3e-6 \
    --g2p G2P_DD_2026-06-11.csv \
    --gff MANE.GRCh38.v1.0.ensembl_genomic.gff3.gz \
    -o report.html
```

The report is not produced by the workflow: run it once the workflow has finished.
Options are described in the [CLI reference](../reference/cli.md#report).

## Requirements

The report reads the workflow's output directory. It needs
`combined_tests/combined_tests_results.tsv`, so the workflow must have run the enrichment
test with `runtype = "both"` **and** the clustering test.

It also uses, when present:

| File | Used for |
|------|----------|
| `simulation/{all-coding,mis}/simulation_logs.json` | Expected and observed counts and scores per consequence |
| `rates/stats/rates_stats.json` | Fraction of sites with a missing score |
| `rates/stats/filtered_rates_stats.json` | Fraction of sites removed by `excluded_regions` |
| `rates/stats/adjusted_rates.json` | Genes whose rates were adjusted |

!!! note
    Per-consequence observed scores, missing-annotation warnings and the list of DNMs
    in significant genes rely on information only written to `simulation_logs.json`
    in debug mode (`enrichment.debug = true`).

## Contents

- **Summary**: number of genes tested, significant (`min_pval` below the threshold) and
  flagged.
- **Gene table**: one row per gene with all p-values, observed / expected ratios, and
  - its [Gene2Phenotype](https://www.ebi.ac.uk/gene2phenotype/) confidence (with `--g2p`):
    filled badges for monoallelic, dashed for biallelic-only genes;
  - the gene panels it belongs to (with `--source-lists`).
- **Gene details** (click a row): DNM counts and scores per consequence, clustering
  details, quality metrics, and, for significant genes not already known as monoallelic
  G2P genes, the list of observed DNMs.

## Caution flags

Genes are flagged when their result should be checked before interpretation:

| Flag | Raised when |
|------|-------------|
| Missing p-values | A p-value is missing where it was expected. |
| Missing annotation | More than 10% of the gene's sites have no score. |
| Filtered regions | More than 10% of the gene's sites were removed by `excluded_regions`. |
| Adjusted rates | The gene's mutation rate was adjusted for positive selection in sperm. |
| Symbol ambiguity | Several gene identifiers share the gene's symbol, which may affect the clustering test (needs `--gff`). |
| G2P overlap | A significant gene that is not in G2P overlaps a G2P gene (needs `--g2p` and `--gff`). |

The page loads its styles and scripts from public CDNs, so it needs internet access to
display properly.
