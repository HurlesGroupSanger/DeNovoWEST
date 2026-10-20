# Recipes

Each recipe shows the `params` block for a common analysis. Combine it with the header
described in [Configuration reference](configuration.md#structure-of-a-configuration-file).

## Build an annotated rates file once

Building and annotating the rates file for all genes is the most expensive step. Do it
once, then reuse the result for every cohort and enrichment setting.

```groovy
params {
  create_rates = true
  annotate_rates = true

  gff = "/data/MANE.GRCh38.v1.0.ensembl_genomic.gff3.gz"
  genome_fasta = "/data/GRCh38.fa"
  excluded_regions = "/data/excluded_regions.bed"

  rates_creation {
    model = "kmer"
    path = "../resources/rates_model_dnn.tsv"
  }

  annotation {
    cadd_file = "/data/cadd/whole_genome_SNVs.tsv.gz"
    dbnsfp = [
      database: "/data/dbNSFP5.3.1a_grch38.gz",
      columns: "AlphaMissense_score,REVEL_score"
    ]
  }
}
```

The result is `rates/merged_rates.tsv.gz` (and `rates/rates_adjusted.vcf.gz` for clustering).

## Analyse a cohort with an existing rates file

```groovy
params {
  annotate_dnm = true
  run_enrichment = true
  run_clustering = true

  rates = "/analyses/out/rates_kmer/rates/merged_rates.tsv.gz"
  dnm = "/data/cohort_dnms.tsv"
  gff = "/data/MANE.GRCh38.v1.0.ensembl_genomic.gff3.gz"
  genome_fasta = "/data/GRCh38.fa"
  excluded_regions = "/data/excluded_regions.bed"

  // Same sources as for the rates file, so that the DNMs get the same columns
  annotation {
    cadd_file = "/data/cadd/whole_genome_SNVs.tsv.gz"
    dbnsfp = [
      database: "/data/dbNSFP5.3.1a_grch38.gz",
      columns: "AlphaMissense_score,REVEL_score"
    ]
  }

  enrichment {
    nmales = 17422
    nfemales = 13636
    score = ["score", "AlphaMissense_score"]   // CADD for all coding, AlphaMissense for missense
    runtype = "both"
    cpus = 4
  }

  clustering {
    runtype = "both"
    protein_structures = "/data/alphafold/UP000005640_9606_HUMAN_v4.tar"
    rates = "/analyses/out/rates_kmer/rates/rates_adjusted.vcf.gz"
  }
}
```

!!! important
    Annotate the DNMs with the **same sources** as the rates file. The enrichment test
    reads the score from the same column in both files.

## Use the Roulette mutation rate model

```groovy
params {
  create_rates = true
  annotate_rates = true

  rates_creation {
    model = "roulette"
    path = "/data/roulette/"     // directory of per-chromosome Roulette VCFs
    scaling_factor = 5.075e-8    // required: default factor, or calibrate it on your cohort (see Mutation rates)
  }
  // ...
}
```

Use the Roulette rates in the clustering test too, with `clustering.rates`.

## Compare several settings

Keep `basedir` and change `run_name` for each setting, then compare the
`enrichment_results.tsv` files:

```groovy
params.basedir = "/analyses/cohort_2026"
params.run_name = "alphamissense_imputed"

params {
  run_enrichment = true
  rates = "/analyses/cohort_2026/out/annotation/rates/merged_rates.tsv.gz"
  dnm = "/analyses/cohort_2026/out/annotation/dnm/dnms_annotated.tsv"

  enrichment {
    nmales = 17422
    nfemales = 13636
    score = "AlphaMissense_score"
    runtype = "mis"
    impute_missing_scores = true
  }
}
```

## Test a configuration quickly

Before a full run, check your configuration on a few genes with few simulations:

```groovy
params {
  gene_list = "my_10_genes.tsv"
  enrichment {
    nsim = 10000
    // ...
  }
}
```

## Reproduce a run exactly

The enrichment test draws random numbers, so p-values vary slightly between runs.
Set `debug = true` and `cpus = 1` in `enrichment` to use a fixed seed. This also writes much more
detailed logs, so only use it on a subset of genes.

## Run a single step without Nextflow

Every step of the workflow is also a [command-line tool](../reference/cli.md), for example:

```bash
denovowest simulation --nmales 17422 --nfemales 13636 --runtype mis \
    dnms_annotated.tsv merged_rates.tsv.gz AlphaMissense_score
```
