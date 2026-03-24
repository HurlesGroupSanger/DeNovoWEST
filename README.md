# DeNovoWEST

DeNovoWEST is a simulation-based framework for testing whether genes show an enrichment of damaging *de novo* mutations (DNMs). The method places different classes of coding and splice variants on a shared severity scale and compares the observed mutational burden in each gene against a gene-specific null model built from expected mutation rates.

The original method is described in [Kaplanis, Samocha, Wiel, Zhang et al. Nature 2020](https://www.nature.com/articles/s41586-020-2832-5). This repository contains the refactored codebase used to prepare inputs, annotate variants, and run the gene-level enrichment simulation.

The historical publication-era implementation is available in the [`v1.0.0` branch](https://github.com/HurlesGroupSanger/DeNovoWEST/tree/v1.0.0).

## Repository structure

The repository has two complementary entrypoints:

- `denovowest` Python package and CLI: reusable commands for annotation and simulation.
- `nextflow/` workflow: orchestration layer for end-to-end or multi-step analyses.

Useful directories:

- `denovowest/`: Python package containing the simulation logic, annotation commands, rates utilities, and shared helpers.
- `nextflow/`: Nextflow pipeline, module definitions, and config files.
- `examples/`: small example inputs and outputs for selected steps.
- `tests/`: Python tests and a few command-line smoke scripts.
- `misc/`: standalone helper scripts that are not part of the main CLI surface.

## Installation

### Workflow-oriented setup

If you mainly want to run the Nextflow workflow, use the provided Conda environment:

```bash
git clone https://github.com/HurlesGroupSanger/DeNovoWEST.git
cd DeNovoWEST
conda env create -f misc/conda/denovowest.yml
conda activate denovowest
```

This is the expected setup for the `nextflow/` workflow.

### Python package setup

If you want to use the Python CLI directly, install the package in your environment:

```bash
pip install -e .
```

This exposes the `denovowest` command defined in `pyproject.toml`.

## How to run DeNovoWEST

### Option 1: Nextflow workflow

Use the Nextflow workflow when you want to run a coordinated multi-step analysis.

```bash
cd nextflow
nextflow run denovowest.nf -c nextflow.config.annotate
```

The exact config file depends on the part of the workflow you want to run.

### Option 2: Python CLI

Use the CLI when you want to run an individual annotation or simulation step.

```bash
denovowest --help
denovowest simulation --help
denovowest annotate-cadd --help
denovowest annotate-custom --help
denovowest annotate-dbnsfp --help
denovowest annotate-vcf --help
```

Currently exposed commands:

- `simulation`: gene-level enrichment testing from observed DNMs and annotated mutation rates.
- `annotate-cadd`: add CADD scores to a rates table.
- `annotate-custom`: add columns from a generic tabix-indexed annotation table.
- `annotate-dbnsfp`: add dbNSFP-derived scores or annotations.
- `annotate-vcf`: add annotations stored in a VCF INFO field.

## Data expectations

Most commands operate on tabular variant files that include at least these columns:

- `gene_id`
- `chrom`
- `pos`
- `ref`
- `alt`

Some steps also expect columns such as `consequence`, `symbol`, or a user-specified score column. Annotation commands generally match records by genomic coordinates and alleles, with optional gene-aware filtering for specific resources.

The simulation command expects:

- an observed DNM table
- a rates table containing all possible mutations for the analyzed genes
- a score column present in both tables
- cohort composition via `--nmales` and `--nfemales`

## Outputs

The simulation command writes:

- a tab-separated results table with gene-level expected score, observed score, and p-value
- a JSON log file containing per-gene simulation details

Annotation commands write tab-separated files that preserve the input rows and append the requested annotation columns.

## Stochastic behaviour

The enrichment test uses stochastic simulation. Small differences between runs are therefore expected unless you deliberately configure the runtime for debugging or replay purposes.

## Development notes

Run the Python tests with:

```bash
python -m pytest -q
```

The most useful fast checks during development are:

```bash
python -m pytest -q tests/test_cli_imports.py
python -m pytest -q tests/test_simulation_core.py
python -m pytest -q tests/misc
```

## Documentation goals for the codebase

A useful mental model for the package is:

- `denovowest.simulation`: prepares DNM and rates inputs, then runs the enrichment model.
- `denovowest.annotations`: attaches external scores or metadata to variant tables.
- `denovowest.rates`: builds and summarizes the expected mutation-rate tables used by the simulation.
- `denovowest.utils`: shared I/O, logging, and parameter definitions.

If you are extending the workflow, prefer treating the CLI surface as the stable public interface and the individual Python modules as implementation details unless a module is explicitly designed for reuse.
