# Command-line interface

The `denovowest` command runs individual steps of the analysis without Nextflow.
The workflow calls the same commands internally.

```bash
denovowest --help
denovowest <command> --help
```

| Command | Description |
|---------|-------------|
| [`simulation`](#simulation) | Gene-level enrichment test. |
| [`annotate-cadd`](#annotate-cadd) | Add CADD scores to a variants file. |
| [`annotate-dbnsfp`](#annotate-dbnsfp) | Add dbNSFP scores to a variants file. |
| [`annotate-custom`](#annotate-custom) | Add columns from a tabix-indexed TSV file. |
| [`annotate-vcf`](#annotate-vcf) | Add INFO fields from a VCF file. |
| [`report`](#report) | Build an HTML report from workflow results. |

"Variants file" means a rates file or a DNM file. Both are tab-separated files with at
least the columns `gene_id`, `chrom`, `pos`, `ref` and `alt` (see [Input files](../user-guide/inputs.md)).

## `simulation`

Tests each gene for an enrichment of DNMs, weighted by a variant score.
See [Enrichment test](../method/enrichment.md).

```bash
denovowest simulation [OPTIONS] DNM RATES SCORE_COLUMN
```

| Argument | Description |
|----------|-------------|
| `DNM` | Annotated DNM file. |
| `RATES` | Annotated rates file. |
| `SCORE_COLUMN` | Column holding the variant score, present in both files. |

| Option | Default | Description |
|--------|---------|-------------|
| `--nmales INTEGER` | **required** | Number of male individuals in the cohort. |
| `--nfemales INTEGER` | **required** | Number of female individuals in the cohort. |
| `--gene-list TEXT` | | Restrict the analysis to the genes in this file. |
| `--runtype` | `all-coding` | Variants tested: `all-coding`, `mis` (missense), `syn` (synonymous), `ptv` (protein-truncating), `protein-altering`. See [Run types](../method/enrichment.md#run-types). |
| `--impute-missing-scores` | off | Replace missing scores by the median score of the possible variants with the same consequence in the same gene, instead of discarding them. |
| `--inframe-missense-ratio FLOAT` | `0.03` | Observed ratio of in-frame indels to missense variants, used to derive the expected rate of in-frame indels. |
| `--inframe-score-quantile FLOAT` | `0.5` | In-frame indels get the score at this quantile of the gene's missense scores. |
| `--frameshift-nonsense-ratio FLOAT` | `1.3` | Observed ratio of frameshift to nonsense variants, used to derive the expected rate of frameshifts. |
| `--frameshift-score-quantile FLOAT` | `0.6` | Frameshifts get the score at this quantile of the gene's nonsense scores. |
| `--nsim INTEGER` | `10000000` | Number of simulations for each number of DNMs considered in a gene. |
| `--pvalcap FLOAT` | `0.01` | Stop simulating a gene once its p-value is known to be above this value. |
| `--jobs INTEGER` | `1` | Number of CPU cores. |
| `--debug` | off | Log detailed per-gene information and fix the random seed. Results are reproducible with `--jobs 1`. |
| `--outdir TEXT` | `./` | Output directory. |
| `--outfile TEXT` | `enrichment_results.tsv` | Name of the results file. |

Writes the results table and `simulation_logs.json`, described in
[Output files](../user-guide/outputs.md#simulation).

```bash
denovowest simulation \
    --nmales 17422 --nfemales 13636 \
    --runtype mis --jobs 4 \
    dnms_annotated.tsv merged_rates.tsv.gz AlphaMissense_score
```

## `annotate-cadd`

```bash
denovowest annotate-cadd RATES CADD OUTPUT
```

| Argument | Description |
|----------|-------------|
| `RATES` | Variants file to annotate. |
| `CADD` | CADD whole-genome scores file, bgzipped and tabix-indexed. |
| `OUTPUT` | Annotated output file. |

## `annotate-dbnsfp`

```bash
denovowest annotate-dbnsfp [OPTIONS] RATES_DNM DBNSFP OUTPUT
```

When a variant has several scores (several transcripts or overlapping genes), the
maximum score among the transcripts present in the GFF is kept. Without `--gff`, the
maximum score for the matching gene is kept.

| Option | Description |
|--------|-------------|
| `-c, --columns TEXT` | dbNSFP columns to extract (comma-separated). |
| `-C, --columns-file TEXT` | File listing the dbNSFP columns to extract, one per line. |
| `--gff TEXT` | GFF file or gffutils database, for transcript-level score selection. |

## `annotate-custom`

```bash
denovowest annotate-custom [OPTIONS] VARIANTS ANNOTATION OUTPUT
```

`ANNOTATION` must be tab-separated, start with the columns `chrom`, `pos`, `ref`, `alt`,
and be indexed with [tabix](https://www.htslib.org/doc/tabix.html).

| Option | Description |
|--------|-------------|
| `-c, --columns TEXT` | Columns to extract (comma-separated). |
| `-C, --columns-file TEXT` | File listing the columns to extract, one per line. |

## `annotate-vcf`

```bash
denovowest annotate-vcf [OPTIONS] VARIANTS ANNOTATION OUTPUT
```

Adds fields from the INFO column of a bgzipped, tabix-indexed VCF.

!!! warning
    This command has been tested with variants files using Ensembl gene identifiers
    and with VCFs such as popEVE and COSMIC. Check the output carefully with other files.

| Option | Description |
|--------|-------------|
| `-c, --columns TEXT` | INFO fields to extract (comma-separated). |
| `-C, --columns-file TEXT` | File listing the INFO fields to extract, one per line. |
| `--match-gene` | Match on gene as well as on position and alleles. |
| `--gff TEXT` | GFF used to build the rates file, to map gene symbols to Ensembl identifiers when the VCF only has symbols. |

## `report`

```bash
denovowest report [OPTIONS] RESULTS_DIR
```

Builds an HTML summary of a workflow run. `RESULTS_DIR` is the workflow's `outdir`.
See [HTML report](../user-guide/report.md).

| Option | Description |
|--------|-------------|
| `-t, --threshold FLOAT` | **Required.** Significance threshold after multiple-testing correction (e.g. `2.5e-6` for a Bonferroni correction over ~20,000 genes). |
| `--g2p FILE` | [Gene2Phenotype](https://www.ebi.ac.uk/gene2phenotype/) CSV file, to show which significant genes are already known. |
| `--gff FILE` | GFF3 file. Avoids false splice-site warnings for single-exon genes and, with `--g2p`, flags new genes that overlap a known one. |
| `--source-lists FILE` | TSV with columns `gene_id` and `source_lists` (semicolon-separated names of gene panels), shown next to each gene. |
| `-o, --output PATH` | Output HTML file. Default: `<RESULTS_DIR>/report.html`. |
