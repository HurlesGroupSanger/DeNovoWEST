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
include { RATES_STATS } from './modules/rates.nf'
include { MERGE_RATES_STATS } from './modules/rates.nf'
include { SUMMARIZE_RATES_STATS } from './modules/rates.nf'

include { FILTER_DNM } from './modules/dnm.nf'
include { FILTER_DNM_GFF } from './modules/dnm.nf'
include { PUBLISH_DNM } from './modules/dnm.nf'

include { BCFTOOLS_CSQ_FULL; BCFTOOLS_CSQ_FULL as DNM_BCFTOOLS_CSQ_FULL} from './modules/annotation.nf'
include { CADD; CADD as DNM_CADD } from './modules/annotation.nf'
include { GNOMAD; GNOMAD as DNM_GNOMAD } from './modules/annotation.nf'
include { CONSTRAINTS; CONSTRAINTS as DNM_CONSTRAINTS} from './modules/annotation.nf'
include { SHET; SHET as DNM_SHET } from './modules/annotation.nf'
include { DBNSFP; DBNSFP as DNM_DBNSFP } from './modules/annotation.nf'
include { CUSTOM; CUSTOM as DNM_CUSTOM } from './modules/annotation.nf'
include { CUSTOM as CUSTOM_2; CUSTOM as DNM_CUSTOM_2 } from './modules/annotation.nf'
include { CUSTOM as CUSTOM_3; CUSTOM as DNM_CUSTOM_3 } from './modules/annotation.nf'
include { VCF; VCF as DNM_VCF } from './modules/annotation.nf'
include { VCF as VCF_2; VCF as DNM_VCF_2 } from './modules/annotation.nf'
include { VCF as VCF_3; VCF as DNM_VCF_3 } from './modules/annotation.nf'


include { SIMULATION; SIMULATION as SIMULATION_NS; SIMULATION as SIMULATION_MIS } from './modules/simulation.nf'
include { MERGE_SIMULATION; MERGE_SIMULATION as MERGE_SIMULATION_NS; MERGE_SIMULATION as MERGE_SIMULATION_MIS } from './modules/simulation.nf'

include { DENOVONEAR_LINEAR; DENOVONEAR_3D; PREPARE_DNM_DENOVONEAR; SPLIT_DNM; SPLIT_GENE_LIST_CLUSTERING;
 MERGE_CLUSTERING; MERGE_CLUSTERING as MERGE_CLUSTERING_LINEAR; MERGE_CLUSTERING as MERGE_CLUSTERING_3D; COMBINE_CLUSTERING; COMBINE_DNE_DNN;
 ADD_GENE_ID; ADD_GENE_ID as ADD_GENE_ID_3D; ADD_GENE_ID as ADD_GENE_ID_LINEAR} from './modules/denovonear.nf'



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

    if (!params.containsKey("run_clustering")) {
      params.run_clustering = true
    }

    // Defines the number of genes to process in one batch
    if (!params.containsKey("split_step")) {
      params.split_step = 100
    }

    // Export or not the weights annotated DNM and rates file
    if (!params.containsKey("export_weighted_dnmrates")) {
      params.export_weighted_dnmrates = false
    }


    // By default we run the non-synonymous test
    if (params.run_simulation) {
      if (!params.containsKey("runtype")) {
        params.runtype = "ns"
      }
    }

    // By default we run the 3D clustering test
    if (params.run_clustering) {
      if (!params.containsKey("clustering_runtype")) {
        params.clustering_runtype = "3D"
      }
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
        rates_annotated_ch = DBNSFP(rates_annotated_ch,  file(params.annotation.dbnsfp_file),  file(params.annotation.dbnsfp_file + ".tbi"), params.annotation.dbnsfp_columns,  gffutils_db_ch.first(), "rates")
      }

      if (params.annotation.containsKey("custom")){

        // TODO : dirty solution to run mutliple annotation from several custom files (max 3)
          // Did not find a better workaround
          cpt = 1
          params.annotation.custom.each { key, custom ->
              
              if (cpt == 1){
                rates_annotated_ch = CUSTOM(
                    rates_annotated_ch,
                    file(custom.path),
                    file(custom.path + ".tbi"),
                    custom.columns,
                    gffutils_db_ch.first(),
                    "rates"
                )
              }
              else if (cpt == 2){
                rates_annotated_ch = CUSTOM_2(
                    rates_annotated_ch,
                    file(custom.path),
                    file(custom.path + ".tbi"),
                    custom.columns,
                    gffutils_db_ch.first(),
                    "rates"
                )
              }
              else if (cpt == 3){
                rates_annotated_ch = CUSTOM_3(
                    rates_annotated_ch,
                    file(custom.path),
                    file(custom.path + ".tbi"),
                    custom.columns,
                    gffutils_db_ch.first(),
                    "rates"
                )
              }
              cpt = cpt + 1
          }
      }

      if (params.annotation.containsKey("vcf")){

          // TODO : dirty solution to run mutliple annotation from several VCF files (max 3)
          // Did not find a better workaround
          cpt = 1
          params.annotation.vcf.each { key, vcf ->
              
              if (cpt == 1){
                rates_annotated_ch = VCF(
                    rates_annotated_ch,
                    file(vcf.path),
                    file(vcf.path + ".tbi"),
                    vcf.columns,
                    gffutils_db_ch.first(),
                    "rates"
                )
              }
              else if (cpt == 2){
                rates_annotated_ch = VCF_2(
                    rates_annotated_ch,
                    file(vcf.path),
                    file(vcf.path + ".tbi"),
                    vcf.columns,
                    gffutils_db_ch.first(),
                    "rates"
                )
              }
              else if (cpt == 3){
                rates_annotated_ch = VCF_3(
                    rates_annotated_ch,
                    file(vcf.path),
                    file(vcf.path + ".tbi"),
                    vcf.columns,
                    gffutils_db_ch.first(),
                    "rates"
                )
              }
              cpt = cpt + 1
          }
      }

      // Merge results
      rates_merged_ch = MERGE_RATES(rates_annotated_ch.collect())

      // Create rates stats
      rates_stats_ch = RATES_STATS(rates_annotated_ch)
      rates_stats_merged_ch = MERGE_RATES_STATS(rates_stats_ch.collect())
      rates_stats_summary_ch = SUMMARIZE_RATES_STATS(rates_stats_merged_ch)

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
          dnm_annotated_ch = DNM_DBNSFP(dnm_annotated_ch,  file(params.annotation.dbnsfp_file),  file(params.annotation.dbnsfp_file + ".tbi"), params.annotation.dbnsfp_columns,  gffutils_db_ch.first(),"dnm")
        }

        if (params.annotation.containsKey("custom")){

          // TODO : dirty solution to run mutliple annotation from several custom files (max 3)
          // Did not find a better workaround
          cpt = 1
          params.annotation.custom.each { key, custom ->
              
              if (cpt == 1){
                dnm_annotated_ch = DNM_CUSTOM(
                    dnm_annotated_ch,
                    file(custom.path),
                    file(custom.path + ".tbi"),
                    custom.columns,
                    gffutils_db_ch.first(),
                    "dnm"
                )
              }
              else if (cpt == 2){
                dnm_annotated_ch = DNM_CUSTOM_2(
                    dnm_annotated_ch,
                    file(custom.path),
                    file(custom.path + ".tbi"),
                    custom.columns,
                    gffutils_db_ch.first(),
                    "dnm"
                )
              }
              else if (cpt == 3){
                dnm_annotated_ch = DNM_CUSTOM_3(
                    dnm_annotated_ch,
                    file(custom.path),
                    file(custom.path + ".tbi"),
                    custom.columns,
                    gffutils_db_ch.first(),
                    "dnm"
                )
              }
              cpt = cpt + 1
          }
        }

        if (params.annotation.containsKey("vcf")){

          // TODO : dirty solution to run mutliple annotation from several VCF files (max 3)
          // Did not find a better workaround
          cpt = 1
          params.annotation.vcf.each { key, vcf ->
              
              if (cpt == 1){
                dnm_annotated_ch = DNM_VCF(
                    dnm_annotated_ch,
                    file(vcf.path),
                    file(vcf.path + ".tbi"),
                    vcf.columns,
                    gffutils_db_ch.first(),
                    "rates"
                )
              }
              else if (cpt == 2){
                dnm_annotated_ch = DNM_VCF_2(
                    dnm_annotated_ch,
                    file(vcf.path),
                    file(vcf.path + ".tbi"),
                    vcf.columns,
                    gffutils_db_ch.first(),
                    "rates"
                )
              }
              else if (cpt == 3){
                dnm_annotated_ch = DNM_VCF_3(
                    dnm_annotated_ch,
                    file(vcf.path),
                    file(vcf.path + ".tbi"),
                    vcf.columns,
                    gffutils_db_ch.first(),
                    "rates"
                )
              }
              cpt = cpt + 1
          }
        }

        PUBLISH_DNM(dnm_annotated_ch)
      }
      // If the DNM file is already annotated, nothing to do
      else
      {
        dnm_annotated_ch = dnm_ch
      }
    }

    //////////////// 
    // SIMULATION //
    ////////////////

    if (params.run_simulation)
    {
      simulation_ch = dnm_annotated_ch.combine(rates_annotated_ch)

      if(params.runtype == "both") {

        SIMULATION_NS(simulation_ch, params.score, params.nmales, params.nfemales, "ns")
        MERGE_SIMULATION_NS(SIMULATION_NS.out.results.collect(), "ns")

        SIMULATION_MIS(simulation_ch, params.score, params.nmales, params.nfemales, "mis")
        MERGE_SIMULATION_MIS(SIMULATION_MIS.out.results.collect(), "mis")

      }
      else {

        SIMULATION(simulation_ch, params.score, params.nmales, params.nfemales, params.runtype)
        MERGE_SIMULATION(SIMULATION.out.results.collect(), params.runtype)

      }



    }


    //////////////// 
    // CLUSTERING //
    ////////////////

    if (params.run_clustering)
    {

      gene_list_clustering_ch = SPLIT_GENE_LIST_CLUSTERING(dnm_annotated_ch, gene_list_ch, params.split_step)
      dnm_split_ch = SPLIT_DNM(gene_list_clustering_ch.toSortedList().flatten(), dnm_annotated_ch.first())
      dnm_dnn_ch = PREPARE_DNM_DENOVONEAR(dnm_split_ch, gffutils_db_ch.first())

      if(params.clustering_runtype == "both") {
        DENOVONEAR_LINEAR(dnm_dnn_ch, file(params.gff), file(params.genome_fasta))
        MERGE_CLUSTERING_LINEAR(DENOVONEAR_LINEAR.out.results.collect(), "linear")
        ADD_GENE_ID_LINEAR(MERGE_CLUSTERING_LINEAR.out,  gffutils_db_ch.first(), "linear")

        DENOVONEAR_3D(dnm_dnn_ch, file(params.gff), file(params.genome_fasta), file(params.protein_structures))
        MERGE_CLUSTERING_3D(DENOVONEAR_3D.out.results.collect(), "3D")
        ADD_GENE_ID_3D(MERGE_CLUSTERING_3D.out,  gffutils_db_ch.first(), "3D")

        COMBINE_CLUSTERING(ADD_GENE_ID_LINEAR.out, ADD_GENE_ID_3D.out)

        // Combine enrichment and clustering results
        if (params.run_simulation) {
          if(params.runtype == "both") {
            COMBINE_DNE_DNN(MERGE_SIMULATION_NS.out, MERGE_SIMULATION_MIS.out,COMBINE_CLUSTERING.out)
          }
        }

      }

      if(params.clustering_runtype == "linear") {

        DENOVONEAR_LINEAR(dnm_dnn_ch, file(params.gff), file(params.genome_fasta))
        MERGE_CLUSTERING_LINEAR(DENOVONEAR_LINEAR.out.results.collect(), "linear")
        ADD_GENE_ID_LINEAR(MERGE_CLUSTERING_LINEAR.out,  gffutils_db_ch.first(), "linear")

        // Combine enrichment and clustering results
        if (params.run_simulation){
          if(params.runtype == "both") {
            COMBINE_DNE_DNN(MERGE_SIMULATION_NS.out, MERGE_SIMULATION_MIS.out,ADD_GENE_ID_LINEAR.out)
          }
        }

      }

      if(params.clustering_runtype == "3D") {

        DENOVONEAR_3D(dnm_dnn_ch, file(params.gff), file(params.genome_fasta),  file(params.protein_structures))
        MERGE_CLUSTERING_3D(DENOVONEAR_3D.out.results.collect(), "3D")
        ADD_GENE_ID_3D(MERGE_CLUSTERING_3D.out,  gffutils_db_ch.first(), "3D")


        // Combine enrichment and clustering results
        if (params.run_simulation) {
          if(params.runtype == "both") {
            COMBINE_DNE_DNN(MERGE_SIMULATION_NS.out, MERGE_SIMULATION_MIS.out,ADD_GENE_ID_3D.out)
          }
        }

      }



    }


}

import java.time.LocalDateTime
import java.time.format.DateTimeFormatter

workflow.onComplete {

    def outDir = params.outdir
    def configOutDir = "${outDir}/configs"

    def timestamp = LocalDateTime.now().format(DateTimeFormatter.ofPattern("yyyy-MM-dd_HH-mm"))
    
    // Create the directory if it doesn't exist
    new File(configOutDir).mkdirs()

    // Get the current git commit hash
    def commitHash = ""
    try {
        def command = "git rev-parse HEAD".execute()
        commitHash = command.text.trim()
    } catch (Exception e) {
        println "Unable to retrieve git commit hash: ${e.message}"
        commitHash = "Unknown"
    }

    // Copy the config file to the output directory
    workflow.configFiles.each{configFile ->
        try {
            def sourceFile = configFile.toString()
            def filename = sourceFile.split('/').last()
            def targetFile = "${configOutDir}/${filename}.${timestamp}"
              
            def command = ["cp", sourceFile, targetFile]
            def process = command.execute()
            process.waitFor()

            if (process.exitValue() != 0) {
                println "Failed to copy ${sourceFile}: ${process.err.text}"
            } else {
                def file = new File(targetFile)
                file.append("\n\n//Version used (git commit): ${commitHash}\n")
                println "Saved config file in ${targetFile}"
            }
        } catch (Exception e) {
            println "Error processing ${configFile}: ${e.message}"
        }
    }

}
