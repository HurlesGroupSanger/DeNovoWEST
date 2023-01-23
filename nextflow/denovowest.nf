#!/usr/bin/env nextflow     
// bsub -R 'select[mem>5000] rusage[mem=5000]' -M 5000 -J denovowest -o /lustre/scratch123/hgi/mdt1/teams/hurles/ed11/DNW/out_rate_creation/after_migration/denovowest.%J.o -e /lustre/scratch123/hgi/mdt1/teams/hurles/ed11/DNW/out_rate_creation/after_migration/denovowest.%J.e "nextflow run denovowest.nf -resume -w /lustre/scratch123/hgi/mdt1/teams/hurles/ed11/DNW/out_rate_creation/after_migration"

nextflow.enable.dsl=2

include { RATES_TO_VCF } from './modules/annotation.nf'
include { BCFTOOLS_CSQ } from './modules/annotation.nf'

/*
 * Build the gffutils database from GFF file
 */
process GFFUTILS_DB {

  input:
  path gff_file

  output:
  path "${gff_file}.db"

  script :
  """
  #!/usr/bin/env python
  import gffutils
  gff_db = gffutils.create_db("$gff_file", "${gff_file}.db", merge_strategy="create_unique")    
  """
}


/*
 * Create gene list from the GFF file
 */
process CREATE_GENE_LIST {

  input:
  path gff_db

  output:
  path "gene_list.tsv"

  script :
  """
  #!/usr/bin/env python
  import gffutils
  gff_db = gffutils.FeatureDB('$gff_db', keep_order=True)
  
  with open("gene_list.tsv", "w") as f :
    for gene in gff_db.all_features(featuretype="gene", order_by="start") :
        print(gene.attributes)
        f.write(str(gene.attributes['ID'][0]) + '\\n') 
  """
}


/*
 * Split gene list for rate creation parallelisation
 */
process SPLIT_GENE_LIST {

    input :
    path gene_list 
    val split_step

    output :
    path 'chunk_*'

    script :
    """
    split -l $split_step $gene_list chunk_
    """

}

/*
 * Create rate file
 */
process RATE_CREATION {

    input :
    path gene_list 
    path gff_db
    path fasta
    path mutation_rate_model

    output :
    path "${gene_list}_mutation/mutation_rates.tsv"

    script :
    """
    create_rates_file.py \
      --gff $gff_db \
      --fasta $fasta \
      --mutation_rate_model $mutation_rate_model \
      --gene_list $gene_list \
      --outdir ${gene_list}_mutation
    """
}

/*
 * Merge rate files 
 */
process MERGE_RATES {

  input : 
  path rates, stageAs : "mutation_rates_*.tsv"
  
  output : 
  path "merged_rates.tsv"

  script : 
  """
  # Concatenates the tsv files and keep only the first header
  awk 'FNR==1 && NR!=1{next;}{print}' $rates > merged_rates.tsv
  """

}



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

    // Merge rates files into one
    rate_merged_ch = MERGE_RATES(rate_creation_ch.collect())

    // Turn rates file into VCF file
    vcf_ch = RATES_TO_VCF(rate_creation_ch, params.genome_fasta + ".fai")

    // Bcftools csq
    bcftools_csq_ch = BCFTOOLS_CSQ(vcf_ch, params.gff, params.genome_fasta)

} 
