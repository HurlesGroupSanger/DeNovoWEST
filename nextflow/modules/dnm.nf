/*
 * Filter the DNM files to retain only the variants in gene found in gene_list
 */
process FILTER_DNM {

  input:
  path dnm
  path gene_list

  output:
  path "dnm_filtered.tsv"

  script :
  """
  filter_dnm.py ${dnm} ${gene_list} --output dnm_filtered.tsv
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
