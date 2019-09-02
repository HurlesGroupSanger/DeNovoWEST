# DeNovoWEST  

`DeNovoWEST` is a simulation-based method to test for a statistically significant enrichment of damaging *de novo* mutations (DNMs) in individual genes. This method scores all classes of variants (e.g. nonsense, missense, splice site) on a unified severity scale based on the empirically-estimated positive predictive value of being pathogenic, and incorporates a gene-based weighting derived from the deficit of protein truncating variants in the general population.  

More information on the method can be found in Kaplanis, Samocha, Wiel, Zhang et al [on bioRxiv](link to come).  

There are three main components of this repository:  
1. `DeNovoWEST` code  
2. Input files needed to recreate the figures in the manuscript  
3. Code to recreate figures in the manuscript  

## DeNovoWEST  

To run `DeNovoWEST`, you will need Python 3.x and R. More details are provided in the `README` in the `DeNovoWEST` directory.  


## Input files  

This directory contains the files needed to recreate some of the main text figures from the Kaplanis, Samocha, Wiel, Zhang et al manuscript. Specifically, the following files are provided:  
1. *De novo* mutations from 31,058 individuals with a developmental disorder
2. Extended `DeNovoWEST` results file with information about per gene mutation rates  
3. Sex information for probands in the study  
4. Positive predictive value weight files that were generated as part of `DeNovoWEST`     
5. Results of downsampling analysis  
6. Results from simulations of modelling the number of remaining HI DD-associated genes  
7. Classifications for likelihood of structural malformation of ultrasound  
8. Gene features used to generate Figure 2(b) in manuscript  
9. Results file from [DeNovoNear](https://github.com/jeremymcrae/denovonear)  
10. Results file from enrichment part of `DeNovoWEST`  
11. Information file needed to apply IHW for our dataset  


## Paper figure code  

For each of the main text figures, there is a separate Rmarkdown file or directory with the code to generate the figures from the provided input files. 

