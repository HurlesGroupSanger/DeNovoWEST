# Configuration reference

The workflow is configured with a [Nextflow configuration file](https://www.nextflow.io/docs/latest/config.html)
passed with `-c`:

```bash
nextflow run denovowest.nf -c my_analysis.config
```

Start from one of the example files in `nextflow/` (see [Quick start](../getting-started/quickstart.md#other-example-configurations)).
The workflow prints all the parameters it received when it starts, and saves a copy
of the configuration file in the output directory at the end of the run.

## Structure of a configuration file

```groovy
includeConfig 'conf/local.conf'     // or conf/lsf.conf, see "Running on a cluster"
includeConfig 'conf/default.conf'

conda.enabled = true                // see Installation for containers

params.basedir  = "/path/to/analyses"
params.run_name = "kmer_cadd"

workDir       = "${params.basedir}/work/${params.run_name}"
params.outdir = "${params.basedir}/out/${params.run_name}"

params {
  // 1. modes: which steps to run
  // 2. shared parameters: input files used by several steps
  // 3. one block per step: rates_creation { }, annotation { }, enrichment { }, clustering { }
}
```

!!! tip
    Keep the same `basedir` for a project and change only `run_name` between
    analyses (e.g. `kmer` vs `roulette` mutation rate model). Each run then gets its own
    output and work directory.

## Modes

Steps are switched on with boolean flags in `params`. They are all `false` by default,
and at least one must be set.

| Parameter | Runs |
|-----------|------|
| `create_rates` | [Rates file creation](#rates-file-creation-rates_creation) from the GFF and a mutation rate model |
| `annotate_rates` | [Annotation](#annotation-annotation) of the rates file |
| `annotate_dnm` | [Annotation](#annotation-annotation) of the DNM file |
| `run_enrichment` | [Enrichment test](#enrichment-test-enrichment) |
| `run_clustering` | [Clustering test](#clustering-test-clustering) |
| `plot_oe_ratios` | [Observed / expected ratio QC](#observed-expected-ratios-oe_ratios) |

Steps can be combined in a single run (e.g. `nextflow.config.full`) or run separately,
feeding the outputs of one run to the next (e.g. annotate once, then try several
enrichment settings).

When both enrichment (with `runtype = "both"`) and clustering are switched on, the
workflow also combines their p-values (see [Output files](outputs.md#combined_tests)).

## Shared parameters

Set at the top level of `params`.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `dnm` | path | — | Observed DNM file. See [Input files](inputs.md#dnm-file). |
| `rates` | path | — | An existing rates file (bgzipped, with a `.tbi` index). When set, `create_rates` is skipped. |
| `gff` | path | — | Gene annotation in GFF3 format. Required to create or annotate rates, annotate DNMs, and for clustering. |
| `genome_fasta` | path | — | Reference genome FASTA, with a `.fai` index next to it. |
| `gene_list` | path | — | Restrict the analysis to these genes (one `gene_id` per line, as in the GFF, version included). If not set, the list is taken from the rates file, or from the GFF. |
| `excluded_regions` | path | — | BED file of regions to remove from both the rates and the DNM files (e.g. low-complexity regions, segmental duplications). |
| `exclude_sex_chromosomes` | bool | `false` | Discard DNMs on chromosomes X and Y. |
| `sperm_positive_selection_adjustment` | path | — | JSON file of factors used to scale mutation rates for genes or sites under positive selection in the male germline. See [Mutation rates](../method/rates.md#adjusting-for-positive-selection-in-sperm). |
| `split_step` | int | `100` | Number of genes per batch. Each batch is processed as a separate job. |
| `gff_db` | path | — | A gffutils database already built from `gff`, to skip rebuilding it. Mostly useful for development and tests. |

### Tools from environment modules

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `useModules` | bool | `false` | Load tools with `module load` before the processes that need them. |
| `bcftoolsModule` | str | — | Module providing bcftools ≥ 1.19. |
| `tabixModule` | str | — | Module providing tabix and bgzip. |
| `bedtoolsModule` | str | — | Module providing bedtools (only needed with `excluded_regions`). |

## Rates file creation (`rates_creation`)

Used when `create_rates = true` and no `rates` file is given. A rates file lists every
possible single-nucleotide substitution in the coding sequence and splice regions of each
gene, with its expected mutation rate. See [Mutation rates](../method/rates.md).

| Parameter | Type | Allowed values | Default | Description |
|-----------|------|----------------|---------|-------------|
| `model` | str | `kmer`, `roulette`, `carlson` | `kmer` | Mutation rate model. |
| `path` | path | — | — | Model data: a trinucleotide rates table for `kmer` (e.g. `resources/rates_model_dnn.tsv`), or the directory of per-chromosome VCFs for `roulette` and `carlson`. |
| `scaling_factor` | float | — | — | Factor by which `roulette` and `carlson` rates are multiplied to obtain per-generation, per-haploid-genome probabilities. **Required** for these two models. See [Mutation rates](../method/rates.md#scaling-factor). |

```groovy
rates_creation {
  model = "roulette"
  path = "/data/roulette/"
  scaling_factor = 5.075e-8   // Roulette default, see Mutation rates
}
```

## Annotation (`annotation`)

Used when `annotate_rates` and/or `annotate_dnm` is set. The rates file and the DNM
file are annotated with the same sources, so that any annotation can later be used as
the score in the enrichment test. Each source is optional.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `annotate_bcftoolscsq` | bool | `true` | Annotate the consequence of each variant (e.g. `missense_variant`, `stop_gained`) with `bcftools csq` and the `gff` file. The enrichment test needs a consequence column. |
| `cadd_file` | path | — | CADD scores (tabix-indexed TSV). Adds the columns `raw` and `score` (PHRED-scaled). |
| `dbnsfp` | map | — | Scores from [dbNSFP](https://www.dbnsfp.org/), see below. |
| `custom` | map | — | Columns from up to three tab-separated, tabix-indexed files, see below. |
| `vcf` | map | — | Fields from the INFO column of up to three VCF files, see below. |

The formats of these files are described in [Input files](inputs.md#annotation-sources).

!!! warning "Contig names"
    CADD, dbNSFP, custom and VCF annotation files must name chromosomes **without** the
    `chr` prefix (`1`, not `chr1`), as distributed by CADD and dbNSFP. The prefix of the
    rates and DNM files does not matter.

### dbNSFP, custom and VCF sources

These sources accept a list of columns to extract. Give it either inline (`columns`,
comma-separated) or in a file with one column name per line (`columns_file`). If neither
is given, all columns are extracted.

| Parameter | Description |
|-----------|-------------|
| `dbnsfp.database` | Path to the genome-wide dbNSFP file (bgzipped, with `.tbi`). |
| `dbnsfp.columns` / `dbnsfp.columns_file` | dbNSFP columns to extract. |
| `custom.<name>.path` | Path to a tabix-indexed TSV whose first columns are `chrom`, `pos`, `ref`, `alt`. |
| `custom.<name>.columns` / `custom.<name>.columns_file` | Columns to extract. |
| `vcf.<name>.path` | Path to a bgzipped, tabix-indexed VCF. |
| `vcf.<name>.columns` / `vcf.<name>.columns_file` | INFO fields to extract. |

`<name>` is a label of your choice. At most **three** `custom` and three `vcf` sources can
be used in one run.

```groovy
annotation {
  cadd_file = "/data/cadd/whole_genome_SNVs.tsv.gz"

  dbnsfp = [
    database: "/data/dbNSFP5.3.1a_grch38.gz",
    columns: "AlphaMissense_score,REVEL_score"
  ]

  custom = [
    "itan": [
      "path": "Itan_missense_full_annotations.tsv.gz",
      "columns_file": "Itan_annotations_to_keep.csv"
    ]
  ]

  vcf = [
    "cosmic": [
      "path": "Cosmic_GenomeScreensMutant_Normal_v101_GRCh38.vcf.gz",
      "columns": "GENOME_SCREEN_SAMPLE_COUNT"
    ],
    "popeve": [
      "path": "popevemap.vcf.gz",
      "columns": "popEVE,EVE,ESM1v"
    ]
  ]
}
```

## Enrichment test (`enrichment`)

Used when `run_enrichment = true`. Requires an annotated rates file and DNM file, either
produced in the same run or given with `rates` and `dnm`. See
[Enrichment test](../method/enrichment.md).

| Parameter | Type | Allowed values | Default | Description |
|-----------|------|----------------|---------|-------------|
| `nmales` | int | — | **required** | Number of male individuals in the cohort. |
| `nfemales` | int | — | **required** | Number of female individuals in the cohort. |
| `score` | str or list | — | **required** | Column of the rates and DNM files holding the variant severity score. With `runtype = "both"`, a list of two columns `[all_coding_score, missense_score]` uses a different score for each test. |
| `runtype` | str | `all-coding`, `mis`, `both` | `all-coding` | Test all coding and splice variants, missense variants only, or run both tests. See [Run types](../method/enrichment.md#run-types). |
| `nsim` | int | — | `10000000` | Number of simulations for each number of DNMs considered in a gene (see [Enrichment test](../method/enrichment.md#how-the-p-value-is-computed)). |
| `impute_missing_scores` | bool | — | `false` | Replace missing scores by the median score of the possible variants with the same consequence in the same gene, instead of discarding them. See [Missing scores](../method/enrichment.md#missing-scores). |
| `cpus` | int | — | `1` | Number of CPUs used by each enrichment job. |
| `debug` | bool | — | `false` | Write detailed per-gene logs and use a fixed random seed, so that results are reproducible (with `cpus = 1`). |

```groovy
enrichment {
  nmales = 17422
  nfemales = 13636
  score = ["CADD_phred", "AlphaMissense_score"]
  runtype = "both"
}
```

!!! note
    The CLI also offers `syn`, `ptv` and `protein-altering` run types and lets you tune
    how indels are scored. These options are not yet exposed in the workflow; use
    [`denovowest simulation`](../reference/cli.md#simulation) directly.

## Clustering test (`clustering`)

Used when `run_clustering = true`. Requires a DNM file and `gff`. See
[Clustering test](../method/clustering.md).

| Parameter | Type | Allowed values | Default | Description |
|-----------|------|----------------|---------|-------------|
| `runtype` | str | `linear`, `3D`, `both` | `3D` | Test clustering along the protein sequence, in the 3D protein structure, or both. |
| `protein_structures` | path | — | — | Archive of AlphaFold predicted structures. **Required** for `3D` and `both`. |
| `rates` | path | — | — | Genome-wide mutation rates VCF (with `.tbi`) used by denovonear instead of its built-in model, e.g. Roulette rates, or the `rates/rates_adjusted.vcf.gz` file produced by the workflow. |

## Observed / expected ratios (`oe_ratios`)

Used when `plot_oe_ratios = true`. Compares the number of observed DNMs with the number
expected from the rates file, per consequence, as a quality check of the mutation rate
model. Uses `enrichment.nmales` and `enrichment.nfemales`.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `gene_categories` | path | — | Optional table assigning genes to categories; ratios are also computed per category. |
| `gnomad_constraints` | path | — | Optional gnomAD constraint table; ratios are also computed per LOEUF and missense constraint bin. |

## Nextflow settings

The example configuration files also set:

- `workDir`: where Nextflow keeps intermediate files. Use a scratch file system and
  delete it once a run is finished.
- `timeline` and `report`: Nextflow execution reports, written to
  `${basedir}/trace/${run_name}/`.

Use `-resume` to restart a run without recomputing completed steps:

```bash
nextflow run denovowest.nf -c my_analysis.config -resume
```
