# Clustering test

The clustering test asks whether the missense DNMs of a gene are closer to each other
than expected by chance. It complements the [enrichment test](enrichment.md): a gene with
a few missense DNMs that all hit the same domain may not be enriched in DNMs, but is
clearly clustered.

DeNovoWEST runs it with [denovonear](https://github.com/jeremymcrae/denovonear).

## Linear and 3D clustering

| `clustering.runtype` | denovonear command | Distance between DNMs |
|----------------------|--------------------|-----------------------|
| `linear` | `denovonear cluster` | along the coding sequence |
| `3D` (default) | `denovonear cluster-structure` | in the predicted protein structure (AlphaFold) |
| `both` | both | 3D result used when available, linear otherwise |

For each gene, denovonear computes a dispersion statistic of the observed DNMs (the
geometric mean of the pairwise distances), and compares it with the distribution
obtained by placing the same number of DNMs at random, weighted by their mutation rates.

3D clustering needs the AlphaFold predicted structures of the human proteome
(`clustering.protein_structures`), and currently assumes GRCh38.

With `runtype = "both"`, `clustering/combined/clustering_results.tsv` reports, for each
gene, the 3D result if there is one (e.g. a structure was available) and the linear
result otherwise. The `mode` column says which.

## Mutation rates

By default denovonear uses its own trinucleotide mutation rate model. Set
`clustering.rates` to a genome-wide mutation rates VCF to use the same rates as the
enrichment test, for example:

- the Roulette rates VCFs;
- `rates/rates_adjusted.vcf.gz`, the VCF version of the rates file produced by the
  workflow (after region filtering and rate adjustment).

## Preparing the DNMs

Before running denovonear, the workflow:

- converts gene identifiers to gene symbols, which denovonear uses (when several gene
  identifiers share a symbol, the one with the most DNMs is kept);
- keeps only missense DNMs for testing;
- maps the results back to gene identifiers.

!!! tip
    DNMs on sex chromosomes can cause errors in denovonear. Set
    `exclude_sex_chromosomes = true` if you run into them.

## Combination with the enrichment test

The missense clustering p-value is combined with the missense enrichment p-value using
Fisher's method. See [Concepts](concepts.md#combining-the-tests) and
[Output files](../user-guide/outputs.md#combined_tests).
