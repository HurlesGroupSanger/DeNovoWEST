# README for provided input files  

This directory contains the files needed to recreate some of the main text figures from the Kaplanis, Samocha, van de Wiel, Zhang et al manuscript. Specifically, the following files are provided:
1. *De novo* mutations from 31,058 individuals with a developmental disorder
2. Extended `DeNovoWEST` results file
3. Sex information for probands
4. Positive predictive value weight files  

## *De novo* mutations

`DDD_RUMC_GDX_denovos_2019_05_15__wweights.txt.gz`

45,221 *de novo* mutations (DNMs) from 31,058 developmental disorder patients.

| Column | Description |
| --- | --- |
| alt  | Alternative allele of the DNM |
| altprop_child | Fraction of reads that are from the alternative allele |
| chrom | Chromosome of the DNM |
| consequence | Mutation consequence of the DNM (e.g. missense_variant, synonymous_variant) |
| cq | Mutation consequence of the DNM (e.g. missense_variant, synonymous_variant) |
| hgnc_id | HGNC ID |
| id | Proband ID |
| maf | Minor allele frequency in gnomAD (release 2.0.2) |
| pos | Position of the DNM |
| prob | Empty column for now |
| raw | CADD version 1.0 raw score |
| ref | Reference allele of the DNM |
| score | CADD version 1.0 Phred score |
| study | Cohort of the proband with the DNM |
| symbol | HGNC symbol |
| constrained | Whether the DNM is in a missense constrained region (if stop_gained or missense) or in a gene with evidence of regional missense constraint (other mutation types) |


## Results of DeNovoWEST  

`extended_denovoWEST_results.tab`  

Results of applying `DeNovoWEST` to the *de novo* mutations (DNMs) from 31,058 developmental disorder patients. For each gene, the count of DNMs is provided, split by mutation type and by center. This extended results file also includes the results from a (previous method)[https://github.com/jeremymcrae/mupit].  

| Column | Description |
| --- | --- |
| symbol | HGNC symbol |
| num_hgnc_id | HGNC ID |
| chr | Chromosome of the gene |
| bp |  Number of coding base pairs in the gene |
| n_exons | Number of exons in the gene |
| cds_length | Number of coding base pairs |
| diag_group | Diagnostic gene list memberships (e.g. consensus, novel) |
| p_all | |
| p_syn | |
| p_mis | |
| p_lof | |
| consensus_gene | Boolean indicator of whether a gene is on the consensus known diagnostic gene list |
| discordant_gene | Boolean indicator of whether a gene is on the discordant known diagnostic gene list |
| DDD_count | Number of DNMs in DDD probands |
| GDX_count | Number of DNMs in GeneDx probands |
| RUMC_count | Number of DNMs in Radboud University Medical Center probands |
| coding_sequence_variant | Number of coding_sequence_variant DNMs | 
| frameshift_variant | Number of frameshift_variant DNMs |
| inframe_deletion | Number of inframe_deletion DNMs |
| inframe_insertion | Number of inframe_insertion DNMs |
| initiator_codon_variant | Number of initiator_codon_variant DNMs |
| missense_variant | Number of missense_variant DNMs |
| splice_acceptor_variant | Number of splice_acceptor_variant DNMs |
| splice_donor_variant | Number of splice_donor_variant DNMs |
| stop_gained | Number of stop_gained DNMs |
| stop_lost | Number of stop_lost DNMs |
| stop_retained_variant | Number of stop_retained_variant DNMs |
| synonymous_variant | Number of synonymous_variant DNMs |
| lofcount | Number of loss-of-function DNMs |
| lofexpected | Expected number of loss-of-function DNMs |
| lofratio | Observed/expected ratio of loss-of-function DNMs |
| misexpected | Expected number of missense DNMs |
| misratio | Observed/expected ration of missense DNMs |
| me | P-value for the enrichment of missense DNMs |
| lofe | P-value for the enrichment of loss-of-function DNMs |
| powmed | |
| s_het | s_het score |
| pLI | gnomAD pLI score |
| observed | |
| expected | |
| wratio | | |
| dne_p | Enrichment p-value generated as part of `DeNovoWEST` |
| dnn_p | Clustering p-value from `denovonear` |
| com_p | Combined p-value of `dne_p` and `dnn_p` |
| min_p | Minimum p-value of `dne_p` and `dnn_p` |
| denovoWEST_p | P-value from `DeNovoWEST` when applied to the full dataset |
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

