# Creation of weights used in DeNovoWEST  

This R markdown shows the process of creating the positive predictive value (PPV) weights that are used in DeNovoWEST, specifically `weights_ppv_2020_01_17.tab`.  

It is currently written based on the provided _de novo_ mutation file (`DDD_RUMC_GDX_denovos_cadd_shet_wweights_2020_01_17.txt.gz`) and rates files that were released as part of the Kaplanis, Samocha, Wiel, Zhang et al manuscript.   

As a warning: loading the rates file can take a while since they are large.

## Required user inputs  

The user will need to specify:
1. The directory where the rates files are stored   
2. The full path to the annotated _de novo_ mutation file  

## Required libraries  

* data.table  
* dplyr  
* ggplot2  
* wesanderson  

And, if you want to recreate a figure included as part of the supplement:
* cowplot  
* tidyr  

## Warnings about applications to other datasets  

**Size of the dataset**  
When we ran the predecessor of this pipeline on a dataset of ~10k trios (from the DDD cohort), we found that the method worked but the estimates of fold-enrichment were much noisier than when we used the full 31k trios in our cohort. Therefore, we warn that weights generated on smaller datasets may not be reliable.	  

**User judgement calls**  
This script cannot be run "out of the box" on a new dataset. It will require some user inputs and judgement calls, specifically around the number of CADD bins used for missense and nonsense variants. For example, if using fewer trios, you will likely need to use fewer CADD bins (e.g. 3 bins instead of 7). The exact number is up to the user.  