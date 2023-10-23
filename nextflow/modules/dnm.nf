
/*
 * Filter the DNM files to retain only the variants in gene found in gene_list
 */
process FILTER_DNM {


  publishDir "${params.outdir}/dnm", mode: 'copy', overwrite: true


  input:
  path dnm
  path gene_list


  output:
  path "dnm_kept.tsv"
  path "dnm_discarded.tsv"


  script :
  """
  filter_dnm.py ${dnm} ${gene_list} --output_kept_dnm dnm_kept.tsv  --output_discarded_dnm dnm_discarded.tsv
  """
}

/*
 * Filter the DNM files to retain only the variants in gene found in gene_list and found in CDS regions of the gff file
 */
process FILTER_DNM_GFF {


  publishDir "${params.outdir}/dnm", mode: 'copy', overwrite: true


  input:
  path dnm
  path gene_list
  path gff_db


  output:
  path "dnm_kept.tsv"
  path "dnm_discarded.tsv"


  script :
  """
  filter_dnm.py ${dnm} ${gene_list} --output_kept_dnm dnm_kept.tsv  --output_discarded_dnm dnm_discarded.tsv --gff ${gff_db}
  """
}


process PUBLISH_DNM {


  publishDir "${params.outdir}/dnm/", mode: 'copy', overwrite: true


  input:
  path dnm

  output :
  path "dnm_annotated.tsv"


  script :
  """
  cp ${dnm} dnm_annotated.tsv
  """
}
