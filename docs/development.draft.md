# Development

## Setup

```bash
git clone https://github.com/HurlesGroupSanger/DeNovoWEST.git
cd DeNovoWEST
git checkout develop
conda env create -f misc/conda/denovowest.yml
conda activate denovowest
pip install -e .
```

Installing in editable mode (`-e`) makes changes to the Python code take effect
without reinstalling.

## Repository layout

| Path | Contents |
|------|----------|
| `denovowest/` | The Python package. |
| `denovowest/simulation/` | Enrichment test: input preparation, score handling, simulations. |
| `denovowest/annotations/` | Annotation commands (CADD, dbNSFP, custom TSV, VCF). |
| `denovowest/rates/` | Rates file creation, adjustment and statistics. |
| `denovowest/dnm/` | DNM filtering. |
| `denovowest/report/` | HTML report and its Jinja2 template. |
| `denovowest/utils/` | Shared I/O, logging and parameter definitions. |
| `nextflow/` | Workflow (`denovowest.nf`), modules, profiles (`conf/`), example configurations and nf-tests. |
| `nextflow/bin/` | Helper scripts called by workflow processes. |
| `examples/` | Small inputs and expected outputs for each step, used by the tests. |
| `tests/` | Python tests. |
| `resources/` | Reference data used by default (trinucleotide model, constraint tables…). |
| `misc/` | Standalone scripts that are not part of the CLI. |
| `docs/` | This documentation. |
| `input/` | Data of the 2020 publication. Not used by DeNovoWEST 2. |

The `denovowest` CLI is the public interface of the package; modules are
implementation details unless documented otherwise.

## Tests

### Python

```bash
python -m pytest -q
```

Tests are tagged with markers, so you can run a subset:

| Marker | Tests |
|--------|-------|
| `unit` | Unit tests of package internals. |
| `smoke` | End-to-end CLI runs on small fixtures. |
| `regression` | Comparison of outputs with known-good reference files. |
| `misc` | Scripts under `misc/`. |
| `functional` | Broader workflow and integration tests. |

```bash
python -m pytest -q -m unit
python -m pytest -q tests/test_simulation_core.py
```

### Nextflow

Workflow processes are tested with [nf-test](https://www.nf-test.com/):

```bash
cd nextflow
./run_tests.sh
```

### Pre-push hook

`misc/hooks/pre-push` runs the nf-tests before every push to `develop`, and refuses
commits whose message starts with `WIP`. Install it with:

```bash
cp misc/hooks/pre-push .git/hooks/pre-push
chmod +x .git/hooks/pre-push
```

## Documentation

This site is built with [MkDocs](https://www.mkdocs.org/) and the
[Material](https://squidfunk.github.io/mkdocs-material/) theme from the `docs/` folder.

```bash
pip install -r docs/requirements.txt
mkdocs serve           # live preview on http://127.0.0.1:8000
mkdocs build --strict  # fails on broken links, as in CI
```

Documentation is versioned with [mike](https://github.com/jimporter/mike) and deployed
to GitHub Pages by `.github/workflows/docs.yml`:

- every push to `develop` updates the `dev` version;
- every version tag (`v2.0.0`, `v2.1.0b1`…) publishes a `<major>.<minor>` version;
  final releases also become `latest`, the default version.

When you change a parameter, an input or an output format, update the corresponding
page in the same pull request.

## Releasing

DeNovoWEST (the method and workflow) and the `denovowest` Python package have separate
version numbers.

**Python package.** The workflow's conda environment (`misc/conda/denovowest.yml`)
installs the package from PyPI, so publish a new package version whenever the workflow
depends on changes to the Python code:

1. update `version` in `pyproject.toml`;
2. build and upload the package:

    ```bash
    python -m build
    python -m twine upload dist/*
    ```

**DeNovoWEST release.** Merge `develop` into `master` and tag the release with the
DeNovoWEST version (`git tag v2.0.0`). Pushing the tag publishes the documentation for
that version.
