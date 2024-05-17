#!/usr/bin/env nextflow     
// bsub -R 'select[mem>5000] rusage[mem=5000]' -M 5000 -J denovowest -o /lustre/scratch123/hgi/mdt1/teams/hurles/ed11/DNW/out_rate_creation/after_migration/denovowest.%J.o -e /lustre/scratch123/hgi/mdt1/teams/hurles/ed11/DNW/out_rate_creation/after_migration/denovowest.%J.e "nextflow run denovowest.nf -resume -w /lustre/scratch123/hgi/mdt1/teams/hurles/ed11/DNW/out_rate_creation/after_migration"

nextflow.enable.dsl=2

include { GFFUTILS_DB } from './modules/rates.nf'
include { CREATE_GENE_LIST_FROM_GFF} from './modules/rates.nf'
include { CREATE_GENE_LIST_FROM_RATES} from './modules/rates.nf'
include { SPLIT_GENE_LIST } from './modules/rates.nf'
include { RATE_CREATION } from './modules/rates.nf'
include { MERGE_RATES } from './modules/rates.nf'
include { SPLIT_RATES } from './modules/rates.nf'

include { FILTER_DNM } from './modules/dnm.nf'
include { FILTER_DNM_GFF } from './modules/dnm.nf'
include { PUBLISH_DNM } from './modules/dnm.nf'


include { BCFTOOLS_CSQ_FULL; BCFTOOLS_CSQ_FULL as DNM_BCFTOOLS_CSQ_FULL} from './modules/annotation.nf'
include { CADD; CADD as DNM_CADD } from './modules/annotation.nf'
include { GNOMAD; GNOMAD as DNM_GNOMAD } from './modules/annotation.nf'
include { CONSTRAINTS; CONSTRAINTS as DNM_CONSTRAINTS} from './modules/annotation.nf'
include { SHET; SHET as DNM_SHET } from './modules/annotation.nf'
include { DBNSFP; DBNSFP as DNM_DBNSFP } from './modules/annotation.nf'


include { GET_EXPECTED_COUNTS } from './modules/weights.nf'
include { MERGE_EXPECTED } from './modules/weights.nf'
include { GET_OBSERVED_COUNTS } from './modules/weights.nf'
include { MERGE_COUNTS } from './modules/weights.nf'
include { LOESS } from './modules/weights.nf'

include { SIMULATION } from './modules/simulation.nf'
include { MERGE_SIMULATION } from './modules/simulation.nf'



workflow{ 

    /////////////////////////// 
    // DEFAULT_CONFIGURATION //
    ///////////////////////////
    

    if (!params.containsKey("create_rates")) {
      params.create_rates = true
    }


    if (!params.containsKey("annotate_rates")) {
      params.annotate_rates = true
    }

    if (!params.containsKey("annotate_dnm")) {
      params.annotate_dnm = true
    }
    if (params.annotate_dnm || params.annotate_rates) {

      if (!params.annotation.containsKey("annotate_bcftoolscsq")) {
        params.annotation.annotate_bcftoolscsq = true
      }
    }
    if (!params.containsKey("run_simulation")) {
      params.run_simulation = true
    }
    if (!params.containsKey("run_weights_creation")) {
      params.run_weights_creation = true
    }


    // Defines the number of genes to process in one batch
    if (!params.containsKey("split_step")) {
      params.split_step = 100
    }



    ///////////////////////////////////////
    // TODO : Check config file validity //
    ///////////////////////////////////////


    ///////////////////////////////// 
    // GFF UTILS DATABASE CREATION //
    /////////////////////////////////


    // DEBUG and TEST only : avoid to create gffutils database every time
    if (params.containsKey("gff_db")) {
      gffutils_db_ch = Channel.fromPath(params.gff_db)
    }
    // Create gffutils database
    else if (params.containsKey("gff")) {

      gff_ch = Channel.fromPath(params.gff)
      gffutils_db_ch = GFFUTILS_DB(gff_ch) 
    }
    

    // If a gene list is provided we use it
    if (params.containsKey("gene_list")) {
      gene_list_ch = Channel.fromPath(params.gene_list)
    }
    // Otherwise we create one from the gff file
    else {

      if (params.containsKey("rates"))
      {
        rates_ch =  Channel.fromPath(params.rates)
        gene_list_ch = CREATE_GENE_LIST_FROM_RATES(rates_ch)
      }
      else if (params.containsKey("gff_db") || params.containsKey("gff"))
      {
        gene_list_ch = CREATE_GENE_LIST_FROM_GFF(gffutils_db_ch)
      }
    }

    // Split gene list in smaller lists
    split_gene_list_ch = SPLIT_GENE_LIST(gene_list_ch, params.split_step)

    ///////////////////////// 
    // RATES FILE CREATION //
    /////////////////////////

    // If the user provide a rates file we split it in smaller rates files
    if (params.containsKey("rates")) {
      
      rates_ch = SPLIT_RATES(split_gene_list_ch.toSortedList().flatten(), file(params.rates))
    }
    // Otherwise we generate rates files from the GFF file
    else if (params.create_rates)
    {

      // Default mutation rate model is kmer trinucleotide model
      if (!params.containsKey("mutation_rate_model_type")) {
        params.mutation_rate_model_type = "kmer"
      }

      // Create rates files (one per split gene list)
      rates_ch = RATE_CREATION(split_gene_list_ch.toSortedList().flatten(), gffutils_db_ch.first(), file(params.genome_fasta), file(params.mutation_rate_model),  params.mutation_rate_model_type)

      // If the point was to generate a rates file with no annotation, we simply merge the unannotated individual rates file
      if (!(params.annotate_rates)) {
            rates_merged_ch = MERGE_RATES(rates_ch.collect())
      }
    }


    /////////////////////////// 
    // RATES FILE ANNOTATION //
    ///////////////////////////
    
    // Annotate rates file
    if (params.containsKey("annotation") and params.containsKey("annotate_rates") and (params.annotate_rates))
    {
      rates_annotated_ch = rates_ch
      
      if (params.annotation.containsKey("annotate_bcftoolscsq") and (params.annotation.annotate_bcftoolscsq)) {
        rates_annotated_ch = BCFTOOLS_CSQ_FULL(rates_annotated_ch, file(params.genome_fasta), file(params.genome_fasta + ".fai"),  file(params.gff), gffutils_db_ch.first(), "rates" )
      }
      
      if (params.annotation.containsKey("cadd_file")){
        rates_annotated_ch = CADD(rates_annotated_ch,  file(params.annotation.cadd_file),  file(params.annotation.cadd_file + ".tbi"), "rates")
      }
      if (params.annotation.containsKey("gene_full_constraints")){
        rates_annotated_ch = CONSTRAINTS(rates_annotated_ch,  file(params.annotation.gene_full_constraints),  file(params.annotation.gene_region_constraints), "rates")
      }
      
      if (params.annotation.containsKey("shet")){
        rates_annotated_ch = SHET(rates_annotated_ch, file(params.annotation.shet), "rates")
      }

      if (params.annotation.containsKey("gnomad_file")){
        rates_annotated_ch = GNOMAD(rates_annotated_ch,  file(params.annotation.gnomad_file),  file(params.annotation.gnomad_file + ".tbi"), "rates")
      }

      if (params.annotation.containsKey("dbnsfp_file")){
        rates_annotated_ch = DBNSFP(rates_annotated_ch,  file(params.annotation.dbnsfp_file),  file(params.annotation.dbnsfp_file + ".tbi"), file(params.annotation.dbnsfp_columns),  gffutils_db_ch.first(), "rates")
      }

      // Merge results
      rates_merged_ch = MERGE_RATES(rates_annotated_ch.collect())
    }
    // If the rates file is already annotated, nothing to do
    else if (params.containsKey("rates")){
        rates_annotated_ch = rates_ch
    }



    ///////////////////////// 
    // DNM FILE ANNOTATION //
    /////////////////////////

    if (params.containsKey("dnm")) {

      dnm_ch = Channel.fromPath(params.dnm)

      if ((params.containsKey("gff_db") || params.containsKey("gff"))) {
          dnm_ch = FILTER_DNM_GFF(dnm_ch, gene_list_ch, gffutils_db_ch)[0]
      }
      else {
           dnm_ch = FILTER_DNM(dnm_ch, gene_list_ch)[0]
      }

      // Annotate DNM file
      if (params.containsKey("annotation") and params.containsKey("annotate_dnm") and (params.annotate_dnm))
      {
        dnm_annotated_ch = dnm_ch

        if (params.annotation.containsKey("annotate_bcftoolscsq") and (params.annotation.annotate_bcftoolscsq)) {
          dnm_annotated_ch = DNM_BCFTOOLS_CSQ_FULL(dnm_annotated_ch, file(params.genome_fasta), file(params.genome_fasta + ".fai"), file(params.gff),  gffutils_db_ch.first(), "dnm" )
        }
        
        if (params.annotation.containsKey("cadd_file")){
          dnm_annotated_ch = DNM_CADD(dnm_annotated_ch,  file(params.annotation.cadd_file), file(params.annotation.cadd_file + ".tbi"), "dnm")
        }

        if (params.annotation.containsKey("gene_full_constraints")){
          dnm_annotated_ch = DNM_CONSTRAINTS(dnm_annotated_ch,  file(params.annotation.gene_full_constraints), file(params.annotation.gene_region_constraints), "dnm")
        }

        if (params.annotation.containsKey("shet")){
          dnm_annotated_ch = DNM_SHET(dnm_annotated_ch,  file(params.annotation.shet), "dnm")
        }

        if (params.annotation.containsKey("gnomad_file")){
          dnm_annotated_ch = DNM_GNOMAD(dnm_annotated_ch,  file(params.annotation.gnomad_file),  file(params.annotation.gnomad_file + ".tbi"), "dnm")
        }

        if (params.annotation.containsKey("dbnsfp_file")){
          dnm_annotated_ch = DNM_DBNSFP(dnm_annotated_ch,  file(params.annotation.dbnsfp_file),  file(params.annotation.dbnsfp_file + ".tbi"), file(params.annotation.dbnsfp_columns),  gffutils_db_ch.first(),"dnm")
        }

        PUBLISH_DNM(dnm_annotated_ch)
      }
      // If the DNM file is already annotated, nothing to do
      else
      {
        dnm_annotated_ch = dnm_ch
      }
    }

    ////////////////////// 
    // WEIGHTS CREATION //
    //////////////////////


    if (params.containsKey("weights"))
    {
      weights_ch = Channel.fromPath(params.weights)
    }
    else if (params.run_weights_creation)
    {
      expected_ch = GET_EXPECTED_COUNTS(rates_annotated_ch, params.nmales, params.nfemales)
      expected_merged_ch = MERGE_EXPECTED(expected_ch.collect())

      observed_ch = GET_OBSERVED_COUNTS(dnm_annotated_ch)

      merged_counts_ch = MERGE_COUNTS(expected_merged_ch, observed_ch).merged_counts_ch

      weights_ch = LOESS(merged_counts_ch).weights_ch
    }

    //////////////// 
    // SIMULATION //
    ////////////////

    if (params.run_simulation)
    {
      simulation_ch = dnm_annotated_ch.combine(rates_annotated_ch)
      simulation_ch = SIMULATION(simulation_ch, params.score, params.nmales, params.nfemales)
      MERGE_SIMULATION(simulation_ch.collect())
    }




} 
