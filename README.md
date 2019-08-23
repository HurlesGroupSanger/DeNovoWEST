# DeNovoWEST  

`DeNovoWEST` is a simulation-based method to test for a statistically significant enrichment of damaging *de novo* mutations (DNMs) in individual genes. This method scores all classes of variants (e.g. nonsense, missense, splice site) on a unified severity scale based on the empirically-estimated positive predictive value of being pathogenic, and incorporates a gene-based weighting derived from the deficit of protein truncating variants in the general population.  

More information on the method can be found in Kaplanis, Samocha, van de Wiel, Zhang et al [on bioRxiv](link to come).  

There are three main components of this repository:  
1. `DeNovoWEST` code  
2. Input files needed to recreate the figures in the manuscript  
3. Code to recreate figures in the manuscript  

## DeNovoWEST  

To run `DeNovoWEST`, you will need Python x.x and R. More details are provided in the `README` in the `DeNovoWEST` directory.  


## Input files  

This directory contains the files needed to recreate some of the main text figures from the Kaplanis, Samocha, van de Wiel, Zhang et al manuscript. Specifically, the following files are provided:  
1. *De novo* mutations from 31,058 individuals with a developmental disorder
2. Extended `DeNovoWEST` results file with information about per gene mutation rates  
3. Sex information for probands in the study  
4. [gnomAD](https://gnomad.broadinstitute.org/) 2.1.1 constraint file  
5. [HGNC](https://www.genenames.org/) file  
6. Positive predictive value weight files that were generated as part of `DeNovoWEST`   


## Paper figure code  

For each of the main text figures, there is a separate Rmarkdown file or directory with the code to generate the figures from the provided input files. 

