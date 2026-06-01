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
    val scaling_factor

    output :
    tuple (path "${gene_list}_mutation/mutation_rates.tsv"),  val(gene_list.name)

    script :
    def sf_arg = scaling_factor ? "--scaling_factor ${scaling_factor}" : ""
    """
    create_rates_file.py \
      --gff $gff_db \
      --fasta $fasta \
      --rates_model_path $mutation_rate_model \
      --gene_list $gene_list \
      --outdir ${gene_list}_mutation \
      --model $mutation_rate_model_type \
      ${sf_arg}

    """
}

/*
 * Merge rate files 
 */
process MERGE_RATES {

  beforeScript = params.useModules
    ? "module load $params.tabixModule"
    : ""

  publishDir "${params.outdir}/rates/", mode: 'copy', overwrite: true

  input : 
  path rates, stageAs : "mutation_rates_*.tsv"
  
  output : 
  path "merged_rates.tsv.gz"
  path "merged_rates.tsv.gz.tbi"

  script : 
  """
  # Concatenates the tsv files and keep only the first header
  awk 'FNR==1 && NR!=1{next;}{print}' $rates > merged_rates.tsv

  # Sort rates file based on chromosome, position and alternate allele
  (head -n 1 merged_rates.tsv && tail -n +2 merged_rates.tsv | sort -k2,2 -k3,3n -k5,5) > sorted_merged_rates.tsv
  
  # Compress and index rates file
  mv sorted_merged_rates.tsv merged_rates.tsv
  bgzip merged_rates.tsv
  tabix -S 1 -s 2 -b 3 -e 3 merged_rates.tsv.gz
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
  tuple path ("${gene_list}_mutation_rates.tsv"), val(gene_list)

  script : 
  """
  # Extract header
  zcat $rates | head -n 1 > ${gene_list}_mutation_rates.tsv
  # Keep only the rows where gene is in gene_list
  zcat $rates | awk 'FNR==NR {values[\$1]; next} \$1 in values' $gene_list - >> ${gene_list}_mutation_rates.tsv
  """

}


process FILTER_RATES_REGION {

  beforeScript = params.useModules
    ? "module load $params.tabixModule;module load $params.bedtoolsModule"
    : ""
  
  publishDir "${params.outdir}/rates", mode: 'copy', overwrite: true

  input :
  path rates
  path index
  path excluded_regions
  path fasta

  output :
  tuple (path "rates_kept.tsv.gz"), path ("rates_kept.tsv.gz.tbi")
  tuple (path "rates_discarded.tsv.gz"), path ("rates_discarded.tsv.gz.tbi")


  script :
  """

  # Resolve/create fasta index (.fai) if needed
  if [[ -f "${fasta}.fai" ]]; then
    fai="${fasta}.fai"
  elif [[ -f "${fasta}" ]]; then
    fai="${fasta}.fai"
    samtools faidx "${fasta}"
  else
    echo "ERROR: fasta file not found: ${fasta}" >&2
    exit 1
  fi

  # Make sure excluded regions are sorted
  bedtools sort -i $excluded_regions -g \$fai > excluded_regions_sorted.bed

  # Create the regions to keep bed file (complement of excluded regions)
  bedtools complement -i excluded_regions_sorted.bed -g \$fai > regions_to_keep.bed

  # Add header
  zcat  "$rates" | head -n 1 > rates_kept.tsv
  zcat  "$rates" | head -n 1 > rates_discarded.tsv

  # Build the list of DNMs to keep/discard based on regions
  tabix -S 1 -R regions_to_keep.bed "$rates" >> rates_kept.tsv
  tabix -S 1 -R excluded_regions_sorted.bed "$rates" >> rates_discarded.tsv

  # Compress and index rates file
  bgzip rates_kept.tsv
  tabix -S 1 -s 2 -b 3 -e 3 rates_kept.tsv.gz

  bgzip rates_discarded.tsv
  tabix -S 1 -s 2 -b 3 -e 3 rates_discarded.tsv.gz

  # Cleanup
  rm regions_to_keep.bed excluded_regions_sorted.bed

  """
}




/*
 * Generate a JSON file containing statistics about the rates 
 */
process RATES_STATS {

  input : 
  path rates

  
  output : 
  path "rates_stats.json"

  script : 
  """
  rates_stats.py $rates rates_stats.json
  """

}

process MERGE_RATES_STATS {

  publishDir "${params.outdir}/rates/", mode: 'copy', overwrite: true

  input : 
  path rates_stats, stageAs : "rates_stats_*.json"
  
  output : 
  path "rates_stats.json"

  script : 
  """
  jq -s 'reduce .[] as \$item ({}; . * \$item)' *.json > rates_stats.json
  #echo "toto" > rates_stats.json
  """

}


process SUMMARIZE_RATES_STATS {

  publishDir "${params.outdir}/rates/", mode: 'copy', overwrite: true

  input : 
  path rates_stats
  
  output : 
  path "summarized_stats"

  script : 
  """
  summarize_stats.py $rates_stats summarized_stats
  """

}