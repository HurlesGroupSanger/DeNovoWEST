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
    for gene in gff_db.all_features(featuretype="gene") :
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

    publishDir "${params.outdir}/rates/split/", mode: 'symlink', overwrite: true

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