
# DNEtest
Conduct gene enrichment test for de novo mutations.

Usage:
```
usage: DNE_test.py [-h] [--weightdic WEIGHTDIC] [--rates RATES]
                   [--denovos DENOVOS] [--output OUTPUT] [--caddweights]
                   [--nsim NSIM] [--nmales NMALES] [--nfemales NFEMALES]

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

This is a dictionary of weights made by looking at enrichment of that specific class of variants across the cohort. 
The following weights files have been included:

*Rates:*

This is a file with site specific mutation rates per gene which is also annotated with CADD score, MAF (gnomad) and a boolean for whether it lies in a missense contrained regions. For non-missense mutations that lie in missense constrained regions these have been annotated as 'True' if the gene has ANY missense constrained region and 'False' if the gene does not have a single missense constrained region. This has been done to calibrate the weights in a sensible way.
Rates files are currently in directory ```/lustre/scratch115/projects/ddd/users/jk18/denovo_enrich/DNE_test/mutation_rates/rates_proteincoding```. They have been divided into 6 to allow for parallelisation. The test is interested in those with CADD and MAF annotation (end in cadd_maf.txt.gz)

*De novos:*

This is a list of de novo mutations from our cohort that have been annotated with CADD score, MAF (gnomad) and a boolean for whether it lies in a constrained region (see explanation in rates file description). 

*Nsim:*

This is the number of simulations you would like to perform for each gene. Generally use 10^7 as this is the decimal point we care for our pvalue

*Nmales, Nfemales:*

The number of males and females in our cohort. This is

### Submit command
To run the test as we have done here run bla: