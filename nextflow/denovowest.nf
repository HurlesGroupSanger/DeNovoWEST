#!/usr/bin/env nextflow     
// bsub -R 'select[mem>5000] rusage[mem=5000]' -M 5000 -J denovowest -o /lustre/scratch123/hgi/mdt1/teams/hurles/ed11/DNW/out_rate_creation/after_migration/denovowest.%J.o -e /lustre/scratch123/hgi/mdt1/teams/hurles/ed11/DNW/out_rate_creation/after_migration/denovowest.%J.e "nextflow run denovowest.nf -resume -w /lustre/scratch123/hgi/mdt1/teams/hurles/ed11/DNW/out_rate_creation/after_migration"

nextflow.enable.dsl=2

include { GFFUTILS_DB } from './modules/rates.nf'
include { CREATE_GENE_LIST } from './modules/rates.nf'
include { SPLIT_GENE_LIST } from './modules/rates.nf'
include { RATE_CREATION } from './modules/rates.nf'
include { MERGE_RATES } from './modules/rates.nf'

include { BCFTOOLS_CSQ_FULL; BCFTOOLS_CSQ_FULL as DNM_BCFTOOLS_CSQ_FULL} from './modules/annotation.nf'
include { CADD; CADD as DNM_CADD } from './modules/annotation.nf'
include { GNOMAD; GNOMAD as DNM_GNOMAD } from './modules/annotation.nf'
include { CONSTRAINTS; CONSTRAINTS as DNM_CONSTRAINTS} from './modules/annotation.nf'
include { SHET; SHET as DNM_SHET } from './modules/annotation.nf'

include { GET_EXPECTED_COUNTS } from './modules/weights.nf'
include { MERGE_EXPECTED } from './modules/weights.nf'
include { GET_OBSERVED_COUNTS } from './modules/weights.nf'
include { MERGE_COUNTS } from './modules/weights.nf'
include { LOESS } from './modules/weights.nf'




workflow{ 


    // DEBUG and TEST only : avoid to create gffutils database every time
    if (params.containsKey("gff_db")) {
      gffutils_db_ch = Channel.fromPath(params.gff_db)
    }
    // Create gffutils database
    else {
      gff_ch = Channel.fromPath(params.gff)
      gffutils_db_ch = GFFUTILS_DB(gff_ch) 
    }

    // If a gene list is provided we use it
    if (params.containsKey("gene_list")) {
      gene_list_ch = Channel.fromPath(params.gene_list)
    }
    // Otherwise we create one from the gff file
    else {
      gene_list_ch = CREATE_GENE_LIST(gffutils_db_ch)
    }

    // Split gene list in smaller lists
    split_gene_list_ch = SPLIT_GENE_LIST(gene_list_ch, params.split_step)

    // Create rates files (one per split gene list)
    rate_creation_ch = RATE_CREATION(split_gene_list_ch.toSortedList().flatten(), gffutils_db_ch.first(), params.genome_fasta, params.mutation_rate_model)
    
    // Annotate rates file
    if (params.containsKey("annotation"))
    {
      rates_bcftools_csq_ch = BCFTOOLS_CSQ_FULL(rate_creation_ch, params.genome_fasta, params.genome_fasta + ".fai", params.gff, gffutils_db_ch.first() )
      rates_cadd_ch = CADD(rates_bcftools_csq_ch, params.annotation.cadd_file, params.annotation.cadd_file + ".tbi")
      rates_gnomad_ch = GNOMAD(rates_cadd_ch, params.annotation.gnomad_file, params.annotation.gnomad_file + ".tbi" )
      rates_constrained_ch = CONSTRAINTS(rates_gnomad_ch, params.annotation.gene_full_constraints, params.annotation.gene_region_constraints )
      rates_shet_ch = SHET(rates_constrained_ch, params.annotation.shet)
    }

    // Merge results
    // rates_merged_ch = MERGE_RATES(rates_shet_ch.collect())

    // Annotate DNM file
    if (params.containsKey("annotation"))
    {
      dnm_bcftools_csq_ch = DNM_BCFTOOLS_CSQ_FULL(params.dnm, params.genome_fasta, params.genome_fasta + ".fai", params.gff,  gffutils_db_ch.first() )
      dnm_cadd_ch = DNM_CADD(dnm_bcftools_csq_ch, params.annotation.cadd_file, params.annotation.cadd_file + ".tbi")
      dnm_gnomad_ch = DNM_GNOMAD(dnm_cadd_ch, params.annotation.gnomad_file, params.annotation.gnomad_file + ".tbi" )
      dnm_constrained_ch = DNM_CONSTRAINTS(dnm_gnomad_ch, params.annotation.gene_full_constraints, params.annotation.gene_region_constraints )
      dnm_shet_ch = DNM_SHET(dnm_constrained_ch, params.annotation.shet)
    }

    // Weights
    if (params.containsKey("weights"))
    {
      expected_ch = GET_EXPECTED_COUNTS(rates_shet_ch, params.weights.n_males, params.weights.n_females)
      expected_merged_ch = MERGE_EXPECTED(expected_ch.collect())

      observed_ch = GET_OBSERVED_COUNTS(dnm_shet_ch)

      merged_counts_ch = MERGE_COUNTS(expected_merged_ch, observed_ch)

      weights_ch = LOESS(merged_counts_ch)
    }




} 
