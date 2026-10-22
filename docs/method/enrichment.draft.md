# Enrichment test

This page describes how `denovowest simulation` tests a gene. For the general idea, read
[Concepts](concepts.md) first.

## Inputs

- the **rates file**, with for each possible substitution of a gene its mutation
  probability (`prob`), its consequence and its score;
- the **DNM file**, with for each observed DNM its consequence and its score;
- the **score column** to use;
- the numbers of **males and females** in the cohort.

Only genes present in both files are tested. Genes are matched on `gene_id`.

## Preparing the variants

1. **Consequences.** When bcftools csq reports several consequences for a variant
   (joined by `&`), the most severe one is kept. Consequences are then grouped into
   classes: `stop_gained` → nonsense, `splice_donor` and `splice_acceptor` → splice LoF,
   in-frame insertions and deletions → in-frame, and so on.
2. **Run type.** Only the consequences of the chosen [run type](#run-types) are kept.
3. **Expected number of DNMs.** Each probability is multiplied by the number of
   haploid genomes transmitted, 2 × (males + females). Positions on chromosome X get an
   additional correction for the sex composition of the cohort and the higher mutation
   rate of the male germline (male-to-female ratio α = 3.4).
4. **Missing scores** are removed or imputed (see [Missing scores](#missing-scores)).
5. **Indels** are added (see [Indels](#indels)).

## Run types

| Run type | Consequences tested |
|----------|---------------------|
| `all-coding` | frameshift, in-frame insertion and deletion, missense, stop gained, splice acceptor and donor, splice region, start lost, stop lost, stop retained, synonymous |
| `mis` | missense |
| `syn` | synonymous |
| `ptv` | stop gained, frameshift, splice acceptor and donor |
| `protein-altering` | missense, in-frame insertion and deletion, stop gained, frameshift, splice acceptor and donor, start lost, stop lost |

The workflow exposes `all-coding`, `mis` and `both` (which runs the first two). The
other run types are available from the [CLI](../reference/cli.md#simulation).

## Indels

Rates files only contain single-nucleotide substitutions, so DeNovoWEST adds, for each
gene, synthetic indel entries whose rates are derived from those of SNVs, with the ratios
observed in large DNM cohorts:

| Indel class | Expected rate | Score |
|-------------|---------------|-------|
| in-frame | 0.03 × missense rate (`--inframe-missense-ratio`) | median of the gene's missense scores (`--inframe-score-quantile 0.5`) |
| frameshift | 1.3 × nonsense rate (`--frameshift-nonsense-ratio`) | 60th percentile of the gene's nonsense scores (`--frameshift-score-quantile 0.6`) |
| splice site | 0.16 × splice LoF SNV rate | median of the gene's splice LoF scores |
| splice region | 0.1 × splice region SNV rate (`all-coding` only) | median of the gene's splice region scores |

Indel entries are only added for the `all-coding`, `ptv` and `protein-altering` run types.
Observed indel DNMs get the score of the synthetic entry of their class in the same gene.

## Missing scores

Variant effect predictors do not score every variant. By default, rates entries and DNMs
without a score are **removed**. Removed DNMs are counted in `simulation_logs.json`
(`nb_observed_dnms_discarded`); a gene whose DNMs are all removed is not tested.

With `impute_missing_scores` (`--impute-missing-scores`), a missing score is instead
**replaced by the median score of the possible substitutions with the same consequence
class in the same gene**, taken from the rates file. If no variant of a gene has a score,
all its variants get a score of 1, which turns the test for that gene into a plain count
of DNMs (`all_scores_missing` in the logs).

## How the p-value is computed

Let *S* be the gene score and *N* the number of DNMs in the gene. Under the null model,
*N* follows a Poisson distribution with mean μ, the expected number of DNMs, and each DNM
is drawn among the possible variants in proportion to their expected rates. The p-value is

```
p = Σn  P(N = n) × P(S ≥ s_obs | N = n)
```

- **Shortcuts.** If the observed score is below the expected score, Σ prob × score, or
  if no DNM is observed, p = 1 and nothing is simulated.
- **n = 0 and n = 1** are computed exactly.
- **n ≥ 2**: *P(S ≥ s<sub>obs</sub> | N = n)* is estimated from `nsim` simulations of
  *n* variants drawn with replacement, weighted by their rates. Terms that are 0 or 1 by
  construction (the *n* highest scores sum to less than s<sub>obs</sub>, or the *n* lowest
  to more), or whose Poisson weight is below 10<sup>-12</sup>, are not simulated.
- **Pseudocount.** Each simulated term is (*k* + 1) / (`nsim` + 1), where *k* is the number
  of simulations at least as extreme as observed, so no term is exactly 0.
- **Range of n.** n goes up to max(2 × observed DNMs, 250). For genes with μ ≤ 50, n
  increases from 2 and the loop stops once the Poisson tail is negligible; for genes with
  μ > 50, n is explored outwards from μ.
- **Early stopping.** Simulations stop once the running p-value exceeds `--pvalcap`
  (default 0.01): the gene is then clearly not significant, and the reported p-value is a
  lower bound just above the cap.
- **Extremely enriched genes.** If no simulation reaches the observed score at all, the
  p-value is set to an upper bound derived from the Poisson probability of the smallest
  number of DNMs that could reach it (floored at 10<sup>-300</sup>).

`nsim` simulations are run **for each value of n** that needs simulating, so the total
number of simulations for a gene (`nb_simulations` in the logs) is usually much larger
than `nsim`.

## Randomness and reproducibility

Simulations are not seeded by default, so p-values vary slightly between runs. With
`--debug` (`enrichment.debug`), the random seed is fixed, which makes results reproducible
when running on a single core (`--jobs 1`).

## Outputs

See [Output files](../user-guide/outputs.md#simulation).
