# Running on a cluster

A full analysis of ~20,000 genes runs a few hundred jobs, some of them long, so it is
usually run on a compute cluster. Nextflow submits the jobs for you: you only need to
tell it which scheduler to use and how many resources each job needs.

## Execution profiles

The workflow ships with two profiles in `nextflow/conf/`. Include one at the top of
your configuration file:

| File | Use |
|------|-----|
| `conf/local.conf` | Run all jobs on the current machine. |
| `conf/lsf.conf` | Submit jobs to an [LSF](https://www.ibm.com/docs/en/spectrum-lsf) cluster. Written for the Wellcome Sanger Institute farm. |

```groovy
includeConfig 'conf/lsf.conf'
includeConfig 'conf/default.conf'
```

### Adapting `lsf.conf`

`lsf.conf` defines resources per type of job:

| Jobs | Queue | Memory | Time |
|------|-------|--------|------|
| default | `normal` | 2 GB | 10 h |
| `process_long` (enrichment and clustering) | `long` | 2 GB | 48 h |
| `intensive` (rates statistics) | `long` | 40 GB | 48 h |

Failed jobs are retried up to twice with more memory and time when they were killed
for exceeding their limits.

Before using it outside the Sanger farm, edit:

- the group passed to LSF in `clusterOptions` (`-G team29-grp`) — replace it with your
  own group (`-G`) or project (`-P`);
- the queue names and the limits in the `check_max` function at the end of the file.

For other schedulers (Slurm, SGE, PBS…), copy `lsf.conf`, change `executor` and the
`clusterOptions`, and see the
[Nextflow executor documentation](https://www.nextflow.io/docs/latest/executor.html).

## Running the workflow itself

The Nextflow process must keep running for the whole analysis, which can take days.
Submit it as a job of its own, or run it inside `screen` / `tmux` on a login node if
your site allows it:

```bash
bsub -q long -R "select[mem>4000] rusage[mem=4000]" -M 4000 \
     -o nextflow.%J.out -e nextflow.%J.err \
     "nextflow run denovowest.nf -c my_analysis.config -resume"
```

## Performance

- **`split_step`** sets how many genes go in each job (default 100). Smaller batches mean
  more, shorter jobs.
- **`enrichment.cpus`** sets the number of cores per enrichment job.
- **`enrichment.nsim`**: the default of 10<sup>7</sup> simulations is needed to reach
  p-values small enough to pass exome-wide significance. Use far fewer (e.g. 10<sup>4</sup>)
  when testing a configuration.
- Put `workDir` on a fast scratch file system: it receives large intermediate files.

## Site-specific notes

### Genomics England Research Environment

A few changes are needed to run on the Genomics England (GEL) HPC:

- **Conda channels** must go through the GEL mirror: in `misc/conda/denovowest.yml`,
  prefix the channels with `conda-`:

    ```yaml
    channels:
      - conda-conda-forge
      - conda-bioconda
    ```

- **bcftools** ≥ 1.19 must be available.
- **LSF**: in `lsf.conf`, use your project code in `clusterOptions` (`-P <project>`)
  and set `queue = 'short'` in the `process` block.
- **Nextflow** ≥ 22 requires `conda.enabled = true` in the configuration file.
