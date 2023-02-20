#!/usr/bin/env nextflow     
// bsub -R 'select[mem>5000] rusage[mem=5000]' -M 5000 -J denovowest -o /lustre/scratch123/hgi/mdt1/teams/hurles/ed11/DNW/out_rate_creation/after_migration/denovowest.%J.o -e /lustre/scratch123/hgi/mdt1/teams/hurles/ed11/DNW/out_rate_creation/after_migration/denovowest.%J.e "nextflow run denovowest.nf -resume -w /lustre/scratch123/hgi/mdt1/teams/hurles/ed11/DNW/out_rate_creation/after_migration"

nextflow.enable.dsl=2

include { GFFUTILS_DB } from './modules/rates.nf'
include { CREATE_GENE_LIST } from './modules/rates.nf'
include { SPLIT_GENE_LIST } from './modules/rates.nf'
include { RATE_CREATION } from './modules/rates.nf'
include { MERGE_RATES } from './modules/rates.nf'

include { RATES_TO_VCF } from './modules/annotation.nf'
include { BCFTOOLS_CSQ } from './modules/annotation.nf'
include { BCFTOOLS_CSQ_FULL } from './modules/annotation.nf'
include { VCF_TO_RATES } from './modules/annotation.nf'
include { CADD } from './modules/annotation.nf'
include { GNOMAD } from './modules/annotation.nf'
include { CONSTRAINTS } from './modules/annotation.nf'




workflow{ 


    // If the input provided is already a gffutils database file, no need to create it
    if (params.gff_file.endsWith(".db")) {
      gffutils_db_ch = Channel.fromPath(params.gff_file)
    }
    // If it is a gff, then create the database
    else {
      gff_ch = Channel.fromPath(params.gff_file)
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
    rates_bcftools_csq_ch = BCFTOOLS_CSQ_FULL(rate_creation_ch, params.genome_fasta, params.genome_fasta + ".fai", params.gff )
    rates_cadd_ch = CADD(rates_bcftools_csq_ch, params.cadd_file, params.cadd_file + ".tbi")
    rates_gnomad_ch = GNOMAD(rates_cadd_ch, params.gnomad_file, params.gnomad_file + ".tbi" )
    rates_constrained_ch = CONSTRAINTS(rates_gnomad_ch, params.gene_full_constraints, params.gene_region_constraints )

    // Merge results
    rates_merged_ch = MERGE_RATES(rates_constrained_ch.collect())

} 
