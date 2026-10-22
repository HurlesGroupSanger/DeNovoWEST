# Output files

All results are written under `params.outdir`. The layout below is for a full run;
each directory only appears if the corresponding step ran.

```
<outdir>/
├── configs/
├── rates/
│   ├── stats/
│   └── split/
├── dnm/
├── simulation/
│   ├── all-coding/
│   └── mis/
├── clustering/
│   ├── linear/
│   ├── 3D/
│   └── combined/
├── combined_tests/
└── oe_ratios/
```

## `configs/`

A copy of each configuration file used for the run, named
`<config file>.<YYYY-MM-DD_HH-MM>`. The git commit of DeNovoWEST used for the run is
appended as the last line. Keep this file with your results: it is the record of how
they were produced.

## `rates/`

| File | When | Description |
|------|------|-------------|
| `gene_list.tsv` | always | Genes analysed. |
| `merged_rates.tsv.gz` (+ `.tbi`) | `create_rates` | The rates file: every possible substitution in the analysed genes, with its mutation rate and, if `annotate_rates` is set, its annotations. This is the file to pass as `rates` to later runs. |
| `rates_adjusted.tsv.gz` (+ `.tbi`) | `sperm_positive_selection_adjustment` | Rates file after adjustment for positive selection in sperm. |
| `rates_adjusted.vcf.gz` (+ `.tbi`) | `annotate_rates` | The final rates file in VCF format, usable as `clustering.rates`. Named `rates_adjusted` even when no adjustment was applied. |
| `split/` | `create_rates` | Unannotated rates for each batch of genes (symbolic links to the work directory). |
| `stats/rates_stats.json` | rates file annotated or provided | Summary statistics of the rates file, including the fraction of variants with missing annotations. |
| `stats/filtered_rates_stats.json` | `excluded_regions` | Number of loci removed by region filtering, overall and per gene. |
| `stats/adjusted_rates.json` | `sperm_positive_selection_adjustment` | Per-gene mutation rates before and after adjustment. |
| `stats/rates_per_gene_per_csq.tsv`, `stats/nb_records_per_gene_per_csq.tsv` | `plot_oe_ratios` | Total mutation rate and number of variants per gene and consequence. |

## `dnm/`

| File | When | Description |
|------|------|-------------|
| `dnms_kept.tsv` | `dnm` set | DNMs kept after filtering on the gene list, coding regions, sex chromosomes and excluded regions. |
| `dnms_discarded.tsv` | `dnm` set | DNMs removed by filtering, with the reason. |
| `dnms_annotated.tsv` | `annotate_dnm` | Kept DNMs with their annotations. This is the file to pass as `dnm` to later runs. |

## `simulation/`

One subdirectory per enrichment run type: `all-coding/` and/or `mis/`.

| File | Description |
|------|-------------|
| `enrichment_results.tsv` | One row per gene with at least one observed DNM. See below. |
| `simulation_logs.json` | Per-gene details of the simulation. See below. |
| `batches/` | Results for each batch of genes (symbolic links to the work directory). |

### `enrichment_results.tsv`

| Column | Description |
|--------|-------------|
| `gene_id` | Gene identifier, as in the rates file. |
| `gene_symbol` | Gene symbol (if the rates file has a `gene_symbol` column). |
| `expected_score` | Expected gene score under the null model. |
| `observed_score` | Sum of the scores of the observed DNMs in the gene. |
| `p-value` | Probability of a gene score at least as high as the observed one under the null model. |

### `simulation_logs.json`

A JSON object keyed by `gene_id`. For each gene:

| Key | Description |
|-----|-------------|
| `nb_observed_dnms`, `nb_expected_dnms` | Number of observed DNMs and expected number of DNMs. |
| `observed_score`, `expected_score` | Observed and expected gene scores. |
| `expected_prob_per_consequence` | Expected number of DNMs per consequence. |
| `expected_score_per_consequence` | Contribution of each consequence to the expected score. |
| `nb_simulations`, `sequential_simulation`, `last_iteration_n` | How many simulations were run for the gene. |
| `info` | A note when the p-value was not simulated, e.g. `observed < expected, pvalue set at 1`. |
| `nb_missing_observed_scores`, `nb_missing_expected_scores`, `nb_missing_mutation_rate` | Variants with missing data. |
| `nb_rates_variants_imputed`, `nb_observed_variants_imputed`, … | Variants whose score was imputed or that were discarded (see `impute_missing_scores`). |
| `wall_time` | Time spent on the gene. |

With `debug = true`, the log also contains the observed DNMs and the result of every
simulation round.

## `clustering/`

One subdirectory per clustering run type: `linear/`, `3D/`, and `combined/` when
`runtype = "both"`.

| Column | Description |
|--------|-------------|
| `gene_id`, `gene_symbol` | Gene identifiers. |
| `mutation_category` | Category of DNMs tested (missense, nonsense…). |
| `events_n` | Number of DNMs tested. |
| `dist` | Observed clustering distance. |
| `probability` | Clustering p-value. |
| `mode` | `combined/` only: whether the result comes from the `3D` or the `linear` test (3D is used when available). |

## `combined_tests/`

Produced when the enrichment test runs with `runtype = "both"` and the clustering test
also runs. `combined_tests_results.tsv` has one row per gene:

| Column | Description |
|--------|-------------|
| `gene_id`, `gene_symbol` | Gene identifiers. |
| `min_pval` | Minimum of `enrichment_all_pval` and `combined_mis_pval`. **This is the DeNovoWEST p-value.** |
| `enrichment_all_pval` | Enrichment test p-value, all coding variants. |
| `combined_mis_pval` | Missense enrichment and clustering p-values combined with Fisher's method. |
| `enrichment_mis_pval` | Enrichment test p-value, missense variants. |
| `clustering_pval` | Clustering test p-value, missense variants. |
| `observed_all`, `expected_all` | Observed and expected gene scores, all coding variants. |
| `observed_mis`, `expected_mis` | Observed and expected gene scores, missense variants. |
| `clustering_nb_mis_variants`, `clustering_dist` | Number of missense DNMs and distance used by the clustering test. |
| `clustering_mode` | `3D` or `linear` (clustering `runtype = "both"` only). |

`min_pval` is the minimum of two tests, so correct for multiple testing accordingly:
the 2020 publication used a Bonferroni threshold of 0.05 / (2 × number of genes tested).
See [Concepts](../method/concepts.md#combining-the-tests).

## `oe_ratios/`

Produced when `plot_oe_ratios = true`.

| File | Description |
|------|-------------|
| `oe_ratios.json` | Observed / expected DNM ratio per consequence. |
| `oe_ratios_per_category.json` | Same, per gene category (if `oe_ratios.gene_categories` is set). |
| `oe_ratios_per_gnomad_loeuf_bins.json`, `oe_ratios_per_gnomad_missense_bins.json` | Same, per gnomAD constraint bin (if `oe_ratios.gnomad_constraints` is set). |
| `plots/*.png` | One plot per consequence and grouping. |

## Outside `outdir`

The example configuration files write Nextflow's execution timeline and resource report
to `${basedir}/trace/${run_name}/`. Intermediate files are in `workDir` and can be deleted
once you no longer need to `-resume` the run.

<!-- Hidden until the HTML report is integrated in the workflow

## HTML report

Once the workflow has finished, summarise the results in an HTML report with
`denovowest report`. See [HTML report](report.md).
-->
