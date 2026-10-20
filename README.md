# DeNovoWEST

**DeNovoWEST** (*De Novo* Weighted Enrichment Simulation Test) identifies genes that carry
more damaging *de novo* mutations (DNMs) than expected by chance in a cohort of
parent–offspring trios.

It combines a gene-level **enrichment test**, which weights each DNM by a variant severity
score and compares the observed burden with simulations from expected mutation rates,
with a **clustering test** of missense DNMs along the protein sequence or in its 3D structure.

📖 **Documentation: <https://hurlesgroupsanger.github.io/DeNovoWEST/>**

> ⚠️ DeNovoWEST **2.0.0** is in beta and developed on the
> [`develop`](https://github.com/HurlesGroupSanger/DeNovoWEST/tree/develop) branch.
> The implementation used in the 2020 publication is on the
> [`v1.0.0`](https://github.com/HurlesGroupSanger/DeNovoWEST/tree/v1.0.0) branch.

## Features

- End-to-end **Nextflow workflow**: from a GFF file and a mutation rate model to gene-level p-values.
- Choice of **mutation rate models**: trinucleotide context, Roulette or Carlson.
- **Any variant score** as severity weight: CADD, AlphaMissense, REVEL, popEVE… from CADD, dbNSFP, VCF or TSV files.
- **Enrichment tests** on all coding variants and on missense variants.
- **Clustering tests** in linear or 3D protein space, through [denovonear](https://github.com/jeremymcrae/denovonear).
- Optional filtering of excluded regions and adjustment of rates for positive selection in sperm.
- An **HTML report** summarising significant genes.
- A **Python command-line tool** to run each step on its own.

## Installation

The workflow needs [Nextflow](https://www.nextflow.io/) ≥ 23.10, conda (or Singularity / Docker),
and bcftools ≥ 1.19, tabix and jq.

```bash
git clone https://github.com/HurlesGroupSanger/DeNovoWEST.git
cd DeNovoWEST
git checkout develop
```

To use the command-line tool on its own:

```bash
pip install denovowest
```

See the [installation guide](https://hurlesgroupsanger.github.io/DeNovoWEST/dev/getting-started/installation/)
for containers and HPC setups.

## Quick start

Run the full workflow on the bundled toy dataset:

```bash
cd nextflow
nextflow run denovowest.nf -c nextflow.config.full
```

Results are written to `nextflow/test/out/full/`. The final per-gene p-values are in
`combined_tests/combined_tests_results.tsv`.

Next, read the [quick start](https://hurlesgroupsanger.github.io/DeNovoWEST/dev/getting-started/quickstart/)
and the [configuration reference](https://hurlesgroupsanger.github.io/DeNovoWEST/dev/user-guide/configuration/)
to run DeNovoWEST on your own data.

## Documentation

| Section | Contents |
|---------|----------|
| [Getting started](https://hurlesgroupsanger.github.io/DeNovoWEST/dev/getting-started/installation/) | Installation, quick start, what's new in v2 |
| [User guide](https://hurlesgroupsanger.github.io/DeNovoWEST/dev/user-guide/overview/) | Workflow steps, configuration, input and output files, cluster usage |
| [Method](https://hurlesgroupsanger.github.io/DeNovoWEST/dev/method/concepts/) | How the enrichment and clustering tests work |
| [CLI reference](https://hurlesgroupsanger.github.io/DeNovoWEST/dev/reference/cli/) | The `denovowest` command-line tool |

## Citation

If you use DeNovoWEST, please cite:

> Kaplanis J, Samocha KE, Wiel L, Zhang Z, *et al.* Evidence for 28 genetic disorders
> discovered by combining healthcare and research data. *Nature* **586**, 757–762 (2020).
> [doi:10.1038/s41586-020-2832-5](https://doi.org/10.1038/s41586-020-2832-5)

## Contributing and support

Bug reports and questions are welcome on the
[issue tracker](https://github.com/HurlesGroupSanger/DeNovoWEST/issues). See the
[development guide](https://hurlesgroupsanger.github.io/DeNovoWEST/dev/development/)
to set up a development environment and run the tests.

## Licence

DeNovoWEST is released under the [MIT licence](https://github.com/HurlesGroupSanger/DeNovoWEST/blob/develop/LICENSE).
