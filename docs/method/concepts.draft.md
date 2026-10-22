# Concepts

DeNovoWEST asks, for each gene: *are there more damaging de novo mutations (DNMs) in this
gene than we would expect by chance, given how mutable it is?*

## The null model: expected mutation rates

Every position of a gene can mutate, at a rate that depends mostly on its local sequence
context. A **mutation rate model** gives, for every possible single-nucleotide substitution,
the probability that it arises *de novo* in one generation. Summed over all the possible
substitutions of a gene, and scaled by the number of individuals in the cohort, it gives the
number of DNMs expected in that gene under the null hypothesis that the gene plays no
role in the disease.

DeNovoWEST stores these probabilities in a **rates file**: one row per possible substitution
in the coding sequence and splice regions of each gene. See [Mutation rates](rates.md).

The expected number of DNMs depends on the sex composition of the cohort, because
chromosome X is present in one copy in males and mutates at a different rate in the male
and female germlines.

## Weighting DNMs by severity

Counting DNMs treats a synonymous variant and a stop-gained variant the same way.
DeNovoWEST instead gives each variant a **severity score**, so that one highly damaging DNM
counts more than one likely benign DNM. The **gene score** is the sum of the scores of the
DNMs observed in the gene.

In the 2020 publication, scores were positive predictive values (PPV) estimated from the
enrichment of DNMs per consequence class and CADD score. DeNovoWEST 2 accepts **any
numerical score** present in the rates and DNM files, typically the output of a
computational variant effect predictor such as CADD, AlphaMissense, REVEL or popEVE.

!!! note "Choosing a score"
    The score must be defined for the variants included in the test: a missense-only
    predictor (e.g. AlphaMissense) suits the missense test, while the all-coding test needs
    a score for all consequence classes (e.g. CADD). Use a list of two scores in
    `enrichment.score` to use a different score for each test. Higher scores must mean
    more damaging.

## The enrichment test

The enrichment test computes the probability, under the null model, of a gene score at
least as high as the one observed. This distribution has no closed form: it depends on how
many DNMs occur (a Poisson variable) and on which variants they are (drawn according to their
mutation rates). DeNovoWEST estimates it by simulation.
See [Enrichment test](enrichment.md).

It is run on two sets of variants:

- **all coding** variants, which has most power for genes where loss-of-function
  variants cause disease;
- **missense** variants only, which has more power for genes where disease is caused by
  missense variants altering protein function.

## The clustering test

Disease-causing missense variants acting through gain of function or dominant negative
effects often cluster in specific protein domains. The clustering test, run with
[denovonear](https://github.com/jeremymcrae/denovonear), asks whether the missense DNMs
of a gene are closer to each other than expected by chance, along the protein sequence or
in its predicted 3D structure. See [Clustering test](clustering.md).

## Combining the tests

For each gene, the workflow produces:

1. `enrichment_all_pval`: enrichment of all coding DNMs;
2. `combined_mis_pval`: the missense enrichment and the clustering p-values combined with
   Fisher's method;
3. `min_pval`: the minimum of (1) and (2), the **DeNovoWEST p-value**.

Because `min_pval` is the best of two tests, the multiple-testing correction must account
for both: the 2020 publication used a Bonferroni threshold of 0.05 / (2 × number of genes
tested).

## Assumptions and caveats

- **Calibration of the mutation rate model.** P-values are only as good as the expected
  rates. Use the [O/E ratio QC](../user-guide/configuration.md#observed-expected-ratios-oe_ratios)
  to check that synonymous DNMs, which are mostly neutral, occur at the expected rate.
- **Regions with unreliable calls.** DNMs in low-complexity regions or segmental
  duplications are enriched in false positives. Remove them from both files with
  `excluded_regions`.
- **Positive selection in the male germline.** A few genes (e.g. in the RAS–MAPK
  pathway) carry mutations that give spermatogonial stem cells a growth advantage, and
  so have more DNMs than their mutation rate predicts. Adjust their rates with
  `sperm_positive_selection_adjustment`. See [Mutation rates](rates.md#adjusting-for-positive-selection-in-sperm).
- **Gene definitions.** The rates file and the DNM file must use the same gene
  identifiers and the same transcripts; the workflow builds both from the same GFF.
