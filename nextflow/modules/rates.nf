/*
 * Build the gffutils database from GFF file
 */
process GFFUTILS_DB {

  publishDir "${params.outdir}/rates", mode: 'symlink', overwrite: true

  input:
  path gff_file

  output:
  path "${gff_file}.db", emit : gffutils_db_ch

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
process CREATE_GENE_LIST_FROM_GFF {

  publishDir "${params.outdir}/rates", mode: 'copy', overwrite: true

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
    for gene in gff_db.all_features(featuretype="gene") :
        f.write(str(gene.attributes['ID'][0]) + '\\n') 
  """
}

/*
 * Create gene list from a provided rates file
 */
process CREATE_GENE_LIST_FROM_RATES {

  publishDir "${params.outdir}/rates", mode: 'copy', overwrite: true

  input:
  path rates_ch

  output:
  path "gene_list.tsv"

  script :
  """
  zcat $rates_ch | tail -n +2 | cut -f1 | sort | uniq > gene_list.tsv
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

    publishDir "${params.outdir}/rates/split/", mode: 'symlink', overwrite: true

    input :
    path gene_list 
    path gff_db
    path fasta
    path mutation_rate_model
    val mutation_rate_model_type

    output :
    path "${gene_list}_mutation/mutation_rates.tsv"

    script :
    """
    create_rates_file.py \
      --gff $gff_db \
      --fasta $fasta \
      --mutation_rate_model $mutation_rate_model \
      --gene_list $gene_list \
      --outdir ${gene_list}_mutation \
      --model $mutation_rate_model_type
    """
}

/*
 * Merge rate files 
 */
process MERGE_RATES {

  publishDir "${params.outdir}/rates/", mode: 'copy', overwrite: true

  input : 
  path rates, stageAs : "mutation_rates_*.tsv"
  
  output : 
  path "merged_rates.tsv.gz"

  script : 
  """
  # Concatenates the tsv files and keep only the first header
  awk 'FNR==1 && NR!=1{next;}{print}' $rates > merged_rates.tsv
  gzip merged_rates.tsv
  """

}

/*
 * Split rate files 
 */
process SPLIT_RATES {

  input : 
  path gene_list
  path rates

  
  output : 
  path "${gene_list}_mutation_rates.tsv"

  script : 
  """
  # Extract header
  zcat $rates | head -n 1 > ${gene_list}_mutation_rates.tsv
  # Keep only the rows where gene is in gene_list
  zcat $rates | awk 'FNR==NR {values[\$1]; next} \$1 in values' $gene_list - >> ${gene_list}_mutation_rates.tsv
  """

}