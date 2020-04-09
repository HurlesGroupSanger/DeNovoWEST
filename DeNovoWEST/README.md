
# DeNovoWEST

There are three steps to running DeNovoWEST:
1. Run enrichment test on all DNMs and the separately on just missense DNMs
2. Run clustering test on DNMs
3. Combine results from tests

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
  --nsim NSIM           number of simulations for each gene
  --nmales NMALES       number of males in study
  --nfemales NFEMALES   number of females in study
  --pvalcap PVALCAP     stop simulations for gene if cumulative p-value > pvalcap (default = 1)
```
### More detailed description of inputs with current paths used:

*Weight Dictionary:* 

This is a dictionary of weights made by looking at enrichment of that specific class of variants across the cohort. We have included this in the input directory. We have included both the weights file that has been generated on all our data (```/input/weights_ppv_2020_01_17.tab``` for all variants and ```/input/weights_ppv_missense_2020_01_17.tab``` which is subset to only missense weights for the missense enrichment run) and the weights file that has been generated after removing variants in the most significant DD-associated genes (```/input/weights_ppv_notop10_2020_01_17.tab```). See manuscript for more details. 

*Rates:*

This is a file with site specific mutation rates per gene which is also annotated with CADD score, MAF (gnomad) and a boolean for whether it lies in a missense contrained regions. For non-missense mutations that lie in missense constrained regions these have been annotated as 'True' if the gene has ANY missense constrained region and 'False' if the gene does not have a single missense constrained region. This has been done to calibrate the weights in a sensible way.
You can download this file from ``` ftp://ftp.sanger.ac.uk/pub/project/ddd/rates ```.

*De novos:*

This is a list of de novo mutations from our cohort that have been annotated with CADD score, MAF (gnomad), a boolean for whether it falls in a high s_het (s_het>0.15) or low s_het (s_het<0.15) gene, and a boolean for whether it lies in a constrained region (see explanation in rates file description). We have provided our list of de novos in the input folder: ```/input/DDD_RUMC_GDX_denovos_cadd_shet_wweights_2020_01_17.txt.gz```.

*Nsim:*

This is the number of simulations you would like to perform for each mutation count for each gene. The default is 10^7.

*Nmales, Nfemales:*

The number of males and females in our cohort. For the total cohort this was 17,422 males and 13,636 females. For the undiagnosed cohort this was 14,087 males and 10,201 females.

### Submit command
To run the enrichment tests (all mutation types and missense only) on all de novo mutations in our joint cohort, we provide a script that prints the necessary commands to run.

Usage:
```
usage: submit_DNE_test_pergene.py [-h] --ratespath RATESPATH --weightspath
                                  WEIGHTSPATH --misweightspath MISWEIGHTSPATH
                                  --denovospath DENOVOSPATH --outpath OUTPATH
                                  --nmale NMALE --nfemale NFEMALE
                                  [--onlyprint]
                                  [--ngenesperfile NGENESPERFILE]
                                  [--pvalcap PVALCAP]

Generate commands per gene for DNE_test.py

optional arguments:
  -h, --help            show this help message and exit
  --ratespath RATESPATH
                        File path to rates file
  --weightspath WEIGHTSPATH
                        File path to weights file
  --misweightspath MISWEIGHTSPATH
                        File path to missense weights file
  --denovospath DENOVOSPATH
                        File path for the de novo variants
  --outpath OUTPATH     Output file path
  --nmale NMALE         Number of males
  --nfemale NFEMALE     Number of females
  --onlyprint           Only print the commands, do not submit
  --ngenesperfile NGENESPERFILE
                        Number of genes per output file
  --pvalcap PVALCAP     Cap for the p-value in DNE test (default is 1.0)
```

For this work, we ran the following command:
```python submit_DNE_test_pergene.py --ratespath <directory of rates> --weightspath /input/weights_ppv_2020_01_17.tab --misweightspath /input/weights_ppv_missense_2020_01_17.tab --denovospath  /input/DDD_RUMC_GDX_denovos_cadd_shet_wweights_2020_01_17.txt.gz --outpath <output directory> --nmale 17422 --nfemale 13636 --ngenesperfile 10```

Results were concatenated and can be found at ```/input/dne_test_2020_01_21_all.tab``` and ```/input/dne_test_2020_01_21_mis.tab```.

## 2. Clustering test

This is done by running [DeNovoNear](https://github.com/jeremymcrae/denovonear). Our results from running this on the full cohort can be found at ```/input/denovonear_out_missense_31058_ntrios_2019_05_15.txt```

## 3. Combine tests

This is done by running the R script as follows:

``` Rscript combine_dne_dnn.R [dne_file_all] [dne_file_missense] [dnn_file] [out_file] ```

Where the dne_file_all is the results from the enrichment test on all variants, the dne_file_missense is the results from the missense only enrichment, the dnn_file is the results from the clustering test, and the out_file is the filename to save it into. To replicate this step for the data from the paper without running the previous step, this can be run as follows:

``` Rscript combine_dne_dnn.R /input/dne_test_2020_01_21_all.tab /input/dne_test_2020_01_21_mis.tab /input/denovonear_out_missense_31058_ntrios_2019_05_15.txt [out_file] ```


## Required libraries  

### Python  

* argparse  
* datetime  
* itertools  
* logging  
* numpy  
* os  
* pandas  
* random  
* scipy (stats)  
* sys  

### R  

* data.table  
* ggplot2  
* gridExtra  
* lattice  
* metap  
* wesanderson  
