# Quick start

This page runs the whole workflow on a toy dataset: a single gene (*OR4F5*) and a
handful of DNMs bundled with the repository. It takes a few minutes and is the best
way to check your [installation](installation.md).

## Run the toy example

```bash
cd DeNovoWEST/nextflow
nextflow run denovowest.nf -c nextflow.config.full
```

`nextflow.config.full` switches on every step:

1. build a mutation rates file for *OR4F5* with the trinucleotide (`kmer`) model;
2. annotate the rates and the DNMs with their consequence (bcftools csq) and a CADD score;
3. run the enrichment test on all coding variants and on missense variants;
4. run the linear clustering test;
5. combine enrichment and clustering p-values.

!!! note "Environment modules"
    The example configuration files load bcftools and tabix through environment modules
    with Sanger-specific module names (`useModules = true`). If these tools are on your
    `PATH`, set `useModules = false`; otherwise change the module names.
    See [Using environment modules](installation.md#using-environment-modules).

## Look at the results

Results are written to `nextflow/test/out/full/`:

```
test/out/full/
├── configs/                 copy of the configuration used, with the git commit
├── rates/                   mutation rates file, annotated and indexed
├── dnm/                     DNMs kept / discarded by filtering, and annotated DNMs
├── simulation/
│   ├── all-coding/          enrichment test on all coding variants
│   └── mis/                 enrichment test on missense variants
├── clustering/linear/       clustering test results
└── combined_tests/          final per-gene p-values
```

The file most people want is `combined_tests/combined_tests_results.tsv`, with one row per
gene and its minimum p-value across tests. All files are described in
[Output files](../user-guide/outputs.md).

Nextflow's intermediate files go to `test/work/full/`, and execution reports
(timeline, resource usage) to `test/trace/full/`.

## Other example configurations

`nextflow/` holds one example configuration per typical use. Copy one and adapt it
for your data:

| File | What it runs |
|------|--------------|
| `nextflow.config.annotate` | create and annotate a rates file, annotate DNMs |
| `nextflow.config.enrichment` | enrichment test from already annotated rates and DNM files |
| `nextflow.config.clustering` | clustering test from already annotated rates and DNM files |
| `nextflow.config.full` | everything above in one run |

## Next steps

- Understand the steps of the workflow: [Workflow overview](../user-guide/overview.md)
- Prepare your own data: [Input files](../user-guide/inputs.md)
- Write your own configuration: [Configuration reference](../user-guide/configuration.md)
