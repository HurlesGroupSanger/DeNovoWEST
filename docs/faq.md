# FAQ and troubleshooting

## Results

### Why do I get slightly different p-values when I rerun an analysis?

The enrichment test estimates p-values by simulation, so they vary a little between
runs. The variation is small compared with the significance thresholds used in practice.
To get identical results, set `enrichment.debug = true` (with `enrichment.cpus = 1`),
which fixes the random seed.

### Why do some genes have a p-value of exactly 1?

When the observed score of a gene is lower than its expected score, there
is no evidence of enrichment and the simulation is skipped. The `info` field of
`simulation_logs.json` says so for each gene.

### Why is a gene missing from the results?

The enrichment test only reports genes with at least one observed DNM that are present
in both the rates and the DNM files. Check:

- `dnm/dnms_discarded.tsv`, which lists the DNMs removed by filtering and why;
- that the gene identifiers in your DNM file match those of the GFF / rates file
  (e.g. Ensembl gene IDs, with or without version);
- `simulation_logs.json` for missing scores or mutation rates.

### What significance threshold should I use?

`min_pval` in `combined_tests_results.tsv` is the minimum of two tests. The 2020
publication used a Bonferroni correction for twice the number of genes tested:
0.05 / (2 × 18,762) ≈ 1.3 × 10<sup>-6</sup>.

### How many simulations do I need?

`nsim` is the number of simulations run for each possible number of DNMs in a gene.
The resolution of each simulated term is about 1 / `nsim`, and each term is weighted by
the Poisson probability of observing that many DNMs, so the default of 10<sup>7</sup>
gives p-values well below exome-wide significance. Lower `nsim` only for testing.

## Installation and running

### `bcftools csq` fails or produces no consequence

DeNovoWEST needs bcftools ≥ 1.19. Check the version available to the jobs
(`bcftools --version`), or the module set in `bcftoolsModule`.

The GFF must be in the format expected by `bcftools csq`. Ensembl GFF3 files work
directly; for other GFF files, see `misc/gff2gff` and [Input files](user-guide/inputs.md#gff).

### Nextflow does not use conda

Since Nextflow 22, conda must be enabled explicitly: add `conda.enabled = true` to
your configuration file.

### Conda cannot reach the default channels

Some secure environments (e.g. Genomics England) mirror conda channels under other
names. Edit the `channels` section of `misc/conda/denovowest.yml`. See
[Running on a cluster](user-guide/hpc.md#genomics-england-research-environment).

### Jobs fail with an out-of-memory or time-limit error

Jobs are retried automatically with more resources with `conf/lsf.conf`. If they still
fail, reduce `split_step` so each job handles fewer genes, or raise the limits in your
profile. See [Running on a cluster](user-guide/hpc.md).

### How do I restart a failed run?

Run the same command with `-resume`. Nextflow reuses the results of the jobs that
completed. Do not delete `workDir` before the run has finished.

### The workflow says "You have to select a pipeline execution mode"

None of the [mode flags](user-guide/configuration.md#modes) is set to `true` in your
configuration file.
