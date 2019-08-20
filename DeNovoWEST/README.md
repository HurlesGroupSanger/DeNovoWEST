
# DeNovoWEST

There are three steps to running DeNovoWEST:
1. Run enrichment test on DNMs
2. Run clustering test on DNMs
3. Combine results from both tests and apply IHW

## 1. Enrichment Test
Conduct gene enrichment test for de novo mutations.

Usage:
```
usage: DNE_test.py [-h] [--weightdic WEIGHTDIC] [--rates RATES]
                   [--denovos DENOVOS] [--output OUTPUT] [--nsim NSIM] 
                   [--nmales NMALES] [--nfemales NFEMALES]

optional arguments:
  -h, --help            show this help message and exit
  --weightdic WEIGHTDIC
                        path to dictionary of weights 
  --rates RATES         path to site specific rates per gene with CADD, MAF
                        and constraint info
  --denovos DENOVOS     path to denovos annotated with constraint and CADD and
                        MAF
  --output OUTPUT       output file to save to
  --nsim NSIM           number of simulations for each mutation score
  --nmales NMALES       number of males in study
  --nfemales NFEMALES   number of females in study
```
### More detailed description of inputs with current paths used:

*Weight Dictionary:* 

This is a dictionary of weights made by looking at enrichment of that specific class of variants across the cohort. We have included this in the input directory. We have included both the weights file that has been generated on all our data and the weights file that has been generated after removing variants in the most significant DD-associated genes (see manuscript for details). 

*Rates:*

This is a file with site specific mutation rates per gene which is also annotated with CADD score, MAF (gnomad) and a boolean for whether it lies in a missense contrained regions. For non-missense mutations that lie in missense constrained regions these have been annotated as 'True' if the gene has ANY missense constrained region and 'False' if the gene does not have a single missense constrained region. This has been done to calibrate the weights in a sensible way.
You can download these rates files HERE. These are split into 6 different files. 

*De novos:*

This is a list of de novo mutations from our cohort that have been annotated with CADD score, MAF (gnomad) and a boolean for whether it lies in a constrained region (see explanation in rates file description). We have provided our list of de novos in the input folder.

*Nsim:*

This is the number of simulations you would like to perform for each mutation count for each gene. The default is 10^7.

*Nmales, Nfemales:*

The number of males and females in our cohort. For the total cohort this was 17422 males and 13636 females

### Submit command
To run enrichment test on all de novo mutations in our joint cohort, the following prints the necessary commands to run:

```python submit_DNE_test.py --ratespath <directory of rates> --weightspath input/weights_ppv_2019_01_09.tab --denovospath  DDD_RUMC_GDX_denovos_2019_05_15__wweights.txt.gz --outpath <output directory> --nmale 17422 --nfemale 13636```


## 2. Clustering test

This is done by running DeNovoNear, this can be found here.

## 3. Combine tests and apply IHW

This is done by running the R script as follows:

``` Rscript combine_IHWcorrect.R [dne_file] [dnn_file] [out_file] ```

Where the dne_file is the results from the enrichment test, the dnn_file is the results from the clustering test and the out_file is the filename to save it into. To replicate this step for the data from the paper without running the previous step, this can be run as follows:

PUT IN ACTUAL PATHS HERE
``` Rscript combine_IHWcorrect.R [dne_file] [dnn_file] [out_file] ```


