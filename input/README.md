# README for provided input files  

This directory contains the files needed to recreate some of the main text figures from the Kaplanis, Samocha, Wiel, Zhang et al Nature 2020 paper. Specifically, the following files are provided:
1. *De novo* mutations from 31,058 individuals with a developmental disorder
2. Extended `DeNovoWEST` results file
3. Sex information for probands
4. Positive predictive value weight files  
5. Results of downsampling analysis
6. Results from simulations of modelling the number of remaining HI DD-associated genes 
7. Classifications for likelihood of structural malformation of ultrasound 
8. Gene features used to generate Figure 2 (b) in manuscript
9. Results file from DeNovoNear
10. Results file from enrichment part of DeNovoWEST
11. List of diagnostic developmental disorders according to DDG2P

## *De novo* mutations

`DDD_RUMC_GDX_denovos_cadd_shet_wweights_2020_01_17.txt.gz`

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
| constrained | Whether the DNM is in a missense constrained region for missense DNMs (should be False for all other consequences) |
| shethigh | Whether the DNM is in a gene with a high s_het, this is defined as s_het > 0.15 |


## Results of DeNovoWEST  

`extended_denovoWEST_results.tab`  

Results of applying `DeNovoWEST` to the *de novo* mutations (DNMs) from 31,058 developmental disorder patients. For each gene, the count of DNMs is provided, split by mutation type and by center. This extended results file also includes the results from a [previous method](https://github.com/jeremymcrae/mupit).  

| Column | Description |
| --- | --- |
| symbol | HGNC symbol |
| num_hgnc_id | HGNC ID |
| chr | Chromosome of the gene |
| bp |  Number of coding base pairs in the gene |
| n_exons | Number of exons in the gene |
| diag_group | Diagnostic gene list memberships (e.g. consensus, novel) |
| p_all | Probability of a single nucleotide mutation based on a [sequence-context based mutational model](https://www.nature.com/articles/ng.3050) of the canonical Ensembl transcript |
| p_syn | Probability of a synonymous mutation based on a [sequence-context based mutational model](https://www.nature.com/articles/ng.3050) of the canonical Ensembl transcript|
| p_mis | Probability of a missense mutation based on a [sequence-context based mutational model](https://www.nature.com/articles/ng.3050) of the canonical Ensembl transcript |
| p_lof | Probability of a loss-of-function mutation based on a [sequence-context based mutational model](https://www.nature.com/articles/ng.3050) of the canonical Ensembl transcript |
| obs_lof | number of observed PTVs in gnomad taken from gnomad constraint file |
| obs_syn | number of observed synonymous in gnomad taken from gnomad constraint file  |
| constraint_flag | flag taken from gnomad constraint file |
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
| misratio | Observed/expected ratio of missense DNMs |
| me | P-value for the enrichment of missense DNMs |
| lofe | P-value for the enrichment of loss-of-function DNMs |
| totaldnm| count of DNMs in this gene |
| powmed | Power to detect median observed PTV enrichment in this gene|
| s_het | s_het score |
| pLI | gnomAD pLI score |
| observed_full | observed gene score for all variants in full cohort|
| expected_full |expected gene score for all variants in full cohort |
| enrich_all_p_full | pEnrich for all variants in full cohort |
| enrich_mis_p_full | pMisEnrich for missense variants in full cohort |
| cluster_p_full | pClustering for full cohort |
| denovoWEST_p_full | DeNovoWEST p-value for full cohort |
| observed_ud | observed gene score for all variants in undiagnosed subset|
| expected_ud |expected gene score for all variants in undiagnosed subset|
| enrich_all_p_ud | pEnrich for all variants in undiagnosed subset|
| enrich_mis_p_ud | pMisEnrich for missense variants in undiagnosed subset|
| cluster_p_ud | pClustering for undiagnosed subset |
| denovoWEST_p_ud | DeNovoWEST p-value for undiagnosed subset |
| denovoWEST_p_ud_synmod | DeNovoWEST p-value for the undiagnosed set using mutation rates modified by the observed to expected ratio of synonymous variants in gnomAD (only for 11 novel genes that are synonymous outliers in gnomAD) |
| gnomad_oe_syn | The observed to expected ratio of synonymous variants from gnomAD v2.1.1, used to modify mutation rates for the 11 novel genes that are synonymous outliers in gnomAD |
| mup_pval | P-value from mupit |
| sig | Boolean indicator of whether gene is significantly associated with DD according to our analysis: this is defined as denovoWEST_p_full< (0.05/(2 x 18,762)) if it is a consensus gene or denovoWEST_p_ud <(0.05/(2 x 18,762))  for non-consensus genes|  
 

## Sex information for individuals in the study  

`fordist_joint_dnm_ID_sex_2019_08_30.txt`

Sex for each proband.  

| Column | Description |
| --- | --- |
| id | Proband ID |
| sex | Sex of the proband |


## Positive predictive value weights  

`weights_ppv_2020_01_17.tab`  

Positive predictive values (PPV) for the *de novo* mutations, split by mutation consequence and constraint strata. These are used as weights in `DeNovoWEST`.    

| Column | Description |
| --- | --- |
| cq | Mutation consequence class |
| score | CADD version 1.0 Phred score |
| constraimed | Boolean for if variant is missense and falls in missense constrained region |
| shethigh | PPV for variants in gene with a high s_het (s_het>=0.15) |
| shetlow | PPV for variants in gene with a low s_het (s_het<0.15) |


`weights_ppv_sanstop10_2020_02_14.tab`  

The above table, but generated after removing the top 10 most mutated genes (*DDX3X*, *ARID1B*, *ANKRD11*, *KMT2A*, *MECP2*, *DYRK1A*, *SCN2A*, *STXBP1*, *MED13L*, *CREBBP*).  

`weights_ppv_missense_2020_01_17.tab`  

The above table, but only containing weights for missense variants. This is used for the generation of missense-specific enrichments.  

## Results of downsampling analysis
 
`downsampling_numbers.tab`

Results from downsampling full cohort to 5k,10k,15k,20k and 25k, rerunning DeNovoWEST and counting number of signifianct genes.

| Column | Description |
| --- | --- |
| samplesize | Sample size the cohort has been downsampled to|
| sig_bon | The number of significant genes after bonferonni correction (DeNovoWEST p-value< (0.05/18762 x 2) ) |
 
 
## Results from simulations of modelling the number of remaining HI DD-associated genes

`PTV_modelresults_2020-03-19.tab`

Results from model simulations using code found at ```/paperfigures/Figure4_code_for_model/PTVmodel.R```. Data used to generate Figure 4(b) in the manuscript.

| Column | Description |
| --- | --- |
| elof_vals | PTV enrichment |
| prophi_vals | Proportion of genes that are haplinsufficient DD genes|
| logprob | log likelihood of observing our observed number of 'de novo' PTVs under this scenario|
| prob |likelihood of observing our observed number of 'de novo' PTVs under this scenario |

## Classifications for likelihood of structural malformation of ultrasound

`USabnormality_likelihood.csv`

Column names are self explanatory. Clinician classified each consensus gene with the likelihood of presenting with a structural malformation of ultrasound. These are classified low (1), medium (2) and high (3).

## Gene features

`gene_features.tab`

Gene features compiled to generate figure 2(b). Column descriptors below:

| Column | Description |
| --- | --- |
| symbol | Gene symbol |
| dNdS_macaque | macaque dN/dS downloaded from Ensembl |
| CODING_GERP | GERP summed across coding region |
| PROMOTER_GERP |GERP summed across promoter region|
| glen | gene length |
| cancer_gene | a binary variable of whether the gene is a known somatic driver gene |
| PPI_DGR | Network degree |
| PPI_BTWN | Network betweeness |
| PPI_LLS2HI | Network distance to consensus known gene |
| median_fetalbrain_rpkm |  the median RPKM in the fetal brain taken from BrainSpan |
| relevant_go |  a binary variable of whether the gene was annotated with one of twenty GO terms that were enriched in consensus DD genes |
| pLI | gnomAD pLI score |

## DeNovoNear results

`denovonear_out_missense_31058_ntrios_2019_05_15.txt`

Output from running [DeNovoNear](https://github.com/jeremymcrae/denovonear) on full cohort. 

## De Novo Enrichment results

`merged_all_dne_test_ppv_2020_03_09.tab`

Output from running enrichment test from DeNovoWEST on all variants in full cohort. See README for DeNovoWEST section for more details.

`merged_mis_dne_test_ppv_2020_03_09.tab`

Output from running enrichment test from DeNovoWEST on only missense variants in full cohort. See README for DeNovoWEST section for more details.

## List of diagnostic developmental disorder genes from DDG2P  

`DDG2P_22_9_2019.csv`  

A curated list of genes associated to developmental disorders downloaded from [DDG2P](https://www.ebi.ac.uk/gene2phenotype/downloads) on 22 September 2019. For each gene-disease pair, the list provides information on the level of certainty of the association, the types of mutations associated, and the allelic requirement. More information can be found [here](https://www.ebi.ac.uk/gene2phenotype/README).  
