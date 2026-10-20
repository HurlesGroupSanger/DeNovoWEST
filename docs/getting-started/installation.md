# Installation

There are two ways to use DeNovoWEST:

- the **Nextflow workflow**, recommended for full analyses: it handles splitting the
  work by gene, running steps in parallel and keeping track of outputs;
- the **Python command-line tool**, to run a single step (e.g. annotate a file or run the
  enrichment test on files you prepared yourself).

## Get the code

```bash
git clone https://github.com/HurlesGroupSanger/DeNovoWEST.git
cd DeNovoWEST
git checkout develop   # DeNovoWEST 2 is currently developed on this branch
```

## Nextflow workflow

### Requirements

| Tool | Version | Notes |
|------|---------|-------|
| [Nextflow](https://www.nextflow.io/docs/latest/install.html) | ≥ 23.10 | |
| conda, Singularity/Apptainer **or** Docker | | provides the Python environment for each process |
| [bcftools](https://samtools.github.io/bcftools/) | ≥ 1.19 | used for consequence annotation (`bcftools csq`) |
| tabix / bgzip ([htslib](https://www.htslib.org/)) | | sorting and indexing rates files |
| [bedtools](https://bedtools.readthedocs.io/) | | only when filtering [excluded regions](../user-guide/configuration.md#shared-parameters) |
| [jq](https://jqlang.github.io/jq/) | | merging simulation logs |

bcftools, tabix, bedtools and jq must be available on your `PATH`, or loaded through
environment modules (see [Using environment modules](#using-environment-modules)).

### Software environment

Choose one of the following and set it in your Nextflow configuration file.

=== "Conda"

    Nextflow builds the environment from `misc/conda/denovowest.yml` the first time the
    workflow runs; `conf/local.conf` already points to it. Enable conda in your
    configuration file:

    ```groovy
    conda.enabled = true
    ```

    !!! tip
        On a cluster, set `conda.cacheDir` (or the `NXF_CONDA_CACHEDIR` environment
        variable) to a shared location so the environment is only built once.

=== "Singularity / Apptainer"

    ```groovy
    conda.enabled = false
    singularity.enabled = true
    process.container = 'docker://delagee/denovowest:4'
    ```

=== "Docker"

    ```groovy
    conda.enabled = false
    docker.enabled = true
    process.container = 'delagee/denovowest:4'
    ```

### Using environment modules

On HPC systems where bcftools, tabix and bedtools are provided as
[environment modules](https://modules.readthedocs.io/), add the module names to the
`params` block of your configuration file instead of installing them:

```groovy
params {
  useModules     = true
  bcftoolsModule = 'bcftools/1.20'
  tabixModule    = 'samtools/1.20'
  bedtoolsModule = 'bedtools/2.31'
}
```

Module names are site-specific; check `module avail` on your system.

## Python package

The `denovowest` package requires Python ≥ 3.10 and is published on
[PyPI](https://pypi.org/project/denovowest/).

!!! info "Version numbers"
    The Python package has its own version number, separate from the version of
    DeNovoWEST (the method and workflow, currently 2.0.0-beta). The conda environment
    used by the workflow installs the package from PyPI.

=== "From PyPI"

    ```bash
    pip install denovowest
    ```

=== "From a clone (development)"

    ```bash
    pip install -e .
    ```

    Changes to the source code then take effect without reinstalling.
    See [Development](../development.md).

Check the installation:

```bash
denovowest --help
```

The clustering test relies on [denovonear](https://github.com/jeremymcrae/denovonear),
which is installed by the conda environment. Install it separately (`pip install denovonear`)
if you use the Python package on its own and want to run clustering.

## Next step

Run the [quick start](quickstart.md) to check that everything works.
