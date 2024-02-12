# DeNovoWEST  

`DeNovoWEST` is a simulation-based method to test for a statistically significant enrichment of damaging *de novo* mutations (DNMs) in individual genes. This method scores all classes of variants (e.g. nonsense, missense, splice site) on a unified severity scale based on the empirically-estimated positive predictive value of being pathogenic, and incorporates a gene-based weighting derived from the deficit of protein truncating variants in the general population.  

More information on the method can be found in [Kaplanis, Samocha, Wiel, Zhang et al Nature 2020](https://www.nature.com/articles/s41586-020-2832-5).  

`DeNovoWEST` version 2.0.0 is a refactoring of the tool that, among other things, includes all the preleminary steps necessary before running the simulation analysis. It is built around nextflow and uses a modular approach to accomodate for each user's specific goal. The version used for publication can be found [here](https://github.com/HurlesGroupSanger/DeNovoWEST/tree/v1.0.0)


## Requirements  

To run `DeNovoWEST` 2.0.0, you need : 
- conda 
- nextflow
- bcftools (>= 1.19)


## Install

```
git clone https://github.com/HurlesGroupSanger/DeNovoWEST.git
git checkout develop
```

## Run

All the parameters are provided via a configuration file that has to be tuned to your own analysis. 

```
cd nextflow
nextflow run denovowest.nf -c nextflow.config.annotate
```

## Conda

By default, nextflow will create the conda environments for you in the nextflow working directory. However this approach has some pitfalls :
- If you remove or change the nextflow working directory, the conda environments will need to be recreated.
- Can be slower

Alternatively, you can install the conda environments outside of nextflow :

```
conda create -f misc/conda/denovowest.yml
conda create -f misc/conda/rpy2.yml
```

and provide the path to those environments in the configuration file (e.g [lsf.conf](nextflow/conf/lsf.conf))


```
// If the conda environments were installed outside of nextflow, you can also provide their path directly
def denovowestEnvLocation = "/software/team29/ed11/miniconda3/envs/denovowest"
def rpy2EnvLocation = "./software/team29/ed11/miniconda3/envs/rpy2"
```