# What's new in DeNovoWEST 2

DeNovoWEST 2 is a rewrite of the code used in the
[2020 publication](https://www.nature.com/articles/s41586-020-2832-5). The statistical
approach is the same; what changed is how the inputs are produced, how variants are
scored, and how the analysis is run.

The 2020 code remains available on the
[`v1.0.0`](https://github.com/HurlesGroupSanger/DeNovoWEST/tree/v1.0.0) branch.

## Highlights

- **A complete workflow.** v1 needed a pre-computed rates file and scripts submitted by hand.
  v2 builds the rates file, annotates it, filters and annotates the DNMs, and runs both tests
  in one Nextflow workflow.
- **Any variant score.** v1 used PPV weights estimated from CADD, s<sub>het</sub> and missense
  constraint in the DDD cohort. v2 uses any numerical score from the rates and DNM files
  (CADD, AlphaMissense, REVEL, popEVE…), with possibly a different score for the missense test.
- **Several mutation rate models**: trinucleotide context, Roulette and Carlson.
- **Integrated clustering test**, in linear or 3D protein space.
- **New options**: excluded regions, adjustment for positive selection in sperm, imputation
  of missing scores, new run types (`syn`, `ptv`, `protein-altering`).
- **Quality control and reporting**: mutation rate statistics, O/E ratios, HTML report.

## Detailed differences

| | v1 (2020) | v2 |
|---|---|---|
| Running | `DNE_test.py` jobs submitted by `submit_DNE_test_pergene.py`; clustering and combination (R) run separately | Nextflow workflow, or `denovowest` CLI |
| Rates file | Pre-computed genome-wide file (GRCh37 and GRCh38 on the Sanger FTP) | Built by the workflow from a GFF, a genome and a mutation rate model |
| Gene identifiers | HGNC symbol | Ensembl `gene_id` (with `gene_symbol`) |
| Consequences | VEP | bcftools csq |
| Variant weights | PPV look-up table (`weights_ppv_*.tab`) | Raw score from a column of the rates/DNM files |
| Missense test | Missense-only weights file | `runtype = "mis"` |
| Indel scores | PPV of the indel class | Quantile of the gene's missense (in-frame) or nonsense (frameshift) scores |
| Splice indels | Not modelled | Modelled from splice SNV rates |
| Missing scores | Variants dropped | Dropped, or imputed with `impute_missing_scores` |
| Number of simulations | 10<sup>9</sup> by default, scaled per number of DNMs, at least 10<sup>7</sup> | 10<sup>7</sup> per number of DNMs |
| Early stopping (`pvalcap`) | Off by default | Stops at p > 0.01 |
| Smallest p-values | Simulated zero replaced by 10<sup>-14</sup> | Pseudocount, and an analytical upper bound for extremely enriched genes |
| Random seed | Always fixed | Fixed only in debug mode |
| Clustering | Linear, run separately | Linear or 3D (AlphaFold), integrated; can use the same rates as the enrichment test |
| Results | `symbol, hgnc_id, expected, observed, p-value, info` | `gene_id, gene_symbol, expected_score, observed_score, p-value` + JSON logs |

## Moving from v1

- **Rates file.** The v1 rates files are keyed by gene symbol and annotated with VEP;
  build a new rates file with the workflow (`create_rates`, `annotate_rates`) instead.
- **Weights.** There is no weights file any more. To reproduce a PPV-like weighting,
  add your weights as a column of the rates and DNM files (e.g. with a
  [custom annotation](../user-guide/configuration.md#dbnsfp-custom-and-vcf-sources)) and use
  it as `enrichment.score`.
- **Results.** Expect differences from v1 results, due to the new scores, consequences,
  gene definitions and simulation settings. Genes that were significant with a large
  margin should remain so.
