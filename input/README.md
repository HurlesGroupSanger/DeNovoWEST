# README for provided input files  

## *De novo* mutations

`DDD_RUMC_GDX_denovos_2019_05_15__wweights.txt.gz`

45,221 *de novo* mutations (DNMs) from 31,058 developmental disorder patients.

| Column | Description |
| --- | --- |
| alt  | Alternative allele of the DNM |
| altprop_child | Fraction of reads that are from the alternative allele |
| chrom | Chromosome of the DNM |
| consequence | Mutation consequence of the DNM () ???? |
| cq | Mutation consequence of the DNM (e.g. missense, synonymous) |
| hgnc_id | HGNC ID |
| id | Proband ID |
| maf | Minor allele frequency in ??? |
| pos | Position of the DNM |
| prob | ?? |
| raw | CADD version 1.0 raw score |
| ref | Reference allele of the DNM |
| score | CADD version 1.0 Phred score |
| study | Cohort of the proband with the DNM |
| symbol | HGNC symbol |
| constrained | Whether the DNM is in a missense constrained region (if stop_gained or missense) or in a gene with evidence of regional missense constraint (other mutation types) |


## Results of DeNovoWEST  

`extended_denovoWEST_results.tab`  

Results of applying `DeNovoWEST` to the *de novo* mutations (DNMs) from 31,058 developmental disorder patients. For each gene, the count of DNMs is provided, split by mutation type and by center. This extended results file also includes the results from a previous method ((mupit)[https://github.com/jeremymcrae/mupit]).  

| Column | Description |
| --- | --- |
| symbol | HGNC symbol |
| num_hgnc_id | HGNC ID |
| chr | Chromosome of the gene |
| bp |  Number of coding base pairs in the gene |
| n_exons | Number of exons in the gene |
| cds_length | ??? |
| diag_group | |
| p_all | |
| p_syn | |
| p_mis | |
| p_lof | |
| consensus_gene | |
| discordant_gene | |
| DDD_count | |
| GDX_count | 
| RUMC_count | 
| coding_sequence_variant | 
| frameshift_variant | 
| inframe_deletion | 
| inframe_insertion | 
| initiator_codon_variant | 
| missense_variant | 
| splice_acceptor_variant | 
| splice_donor_variant | 
| stop_gained | 
| stop_lost | 
| stop_retained_variant | 
| synonymous_variant | 
| lofcount | 
| lofexpected | 
| lofratio | 
| misexpected | |
| misratio | |
| me | 
| lofe | 
| powmed | 
| s_het | 
| pLI | 
| observed | 
| expected | 
| wratio | 
| dne_p | 
| dnn_p | 
| com_p | 
| min_p | 
| denovoWEST_p | 
| undiagnosed_denovoWEST_p | P-value from `DeNovoWEST` when applied to only undiagosed cases |
| mup_pval | P-value from mupit |
| sig | Boolean indicator of whether the gene was significant in either the full of undiagnosed analysis |  
 

## Sex information for individuals in the study  

`fordist_joint_dnm_ID_sex_2019_08_22.txt`

Sex for each proband.  

| Column | Description |
| --- | --- |
| id | Proband ID |
| sex | Sex of the proband |
 

## gnomAD constraint information  

`gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz`

Downloaded from [the gnomAD downloads page](https://gnomad.broadinstitute.org/downloads#gene-constraint).  

See their website for more information.  


## HGNC file   

`protein-coding_gene.txt`  

Downloaded from [the HGNC Statistics & download files page](https://www.genenames.org/download/statistics-and-files/).  

See their website for more information.  


## Positive predictive value weights  

`weights_ppv_2019_01_09.tab`  

Positive predictive values (PPV) for the *de novo* mutations, split by mutation consequence and constraint strata. These are used as weights in `DeNovoWEST`.    

| Column | Description |
| --- | --- |
| cq | Mutation consequence class |
| score | CADD version 1.0 Phred score |
| con | PPV for variants in missense constrained regions/genes |
| uncon | PPV for variants not in missense constrained regions/genes |


`weights_ppv_notop10_2019_01_09.tab`  

The above table, but generated after removing the top 10 most mutated genes (*DDX3X*, *ARID1B*, *ANKRD11*, *KMT2A*, *MECP2*, *DYRK1A*, *SCN2A*, *STXBP1*, *MED13L*, *CREBBP*).  

