
/*
 * Filter the DNM files to retain only the variants in gene found in gene_list
 */
process FILTER_DNM {


  input:
  path dnm
  path gene_list


  output:
  tuple path ("dnms_kept_after_gene_filtering.tsv"), path ("dnms_discarded_after_gene_filtering.tsv")


  beforeScript "[ -v NF_TEST ] && export PYTHONPATH=$baseDir/../../../../../denovowest/;"

  script :
  """
  filter_dnm.py ${dnm} ${gene_list} --output_kept_dnm dnms_kept_after_gene_filtering.tsv  --output_discarded_dnm dnms_discarded_after_gene_filtering.tsv
  """
}

/*
 * Filter the DNM files to retain only the variants in gene found in gene_list and found in CDS regions of the gff file
 */
process FILTER_DNM_GFF {


  input:
  path dnm
  path gene_list
  path gff_db


  output:
  tuple path ("dnms_kept_after_cds_filtering.tsv"), path ("dnms_discarded_after_cds_filtering.tsv")

  script :
  """
  filter_dnm.py ${dnm} ${gene_list} --output_kept_dnm dnms_kept_after_cds_filtering.tsv  --output_discarded_dnm dnms_discarded_after_cds_filtering.tsv --gff ${gff_db}
  """
}

process FILTER_DNM_REGION {

  beforeScript = params.useModules
    ? "module load $params.tabixModule;module load $params.bedtoolsModule"
    : ""
  

  input :
  tuple path (kept_dnm), path(discarded_dnm)
  path excluded_regions
  path fasta

  output :
  tuple path ("dnms_kept_after_region_filtering.tsv"), path ("dnms_discarded_after_region_filtering.tsv")

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

  # Create the regions to keep bed file (complement of excluded regions)
  bedtools complement -i ${excluded_regions} -g \$fai > regions_to_keep.bed

  # Ensure DNM is bgzipped & indexed
  if [[ "${kept_dnm}" != *.gz ]]; then
    bgzip -c "${kept_dnm}" > "${kept_dnm}.gz"
    tabix -S 1 -s 2 -b 3 -e 3 "${kept_dnm}.gz"
    dnm_gz="${kept_dnm}.gz"
  else
    dnm_gz="${kept_dnm}"
  fi

  # Add header
  zcat  "\${dnm_gz}" | head -n 1 > dnms_kept_after_region_filtering.tsv
  zcat  "\${dnm_gz}" | head -n 1 > dnms_discarded_after_region_filtering.tsv

  # Build the list of DNMs to keep/discard based on regions
  tabix -R regions_to_keep.bed "\$dnm_gz" >> dnms_kept_after_region_filtering.tsv
  tabix -R ${excluded_regions} "\$dnm_gz" >> dnms_discarded_after_region_filtering.tsv

  # Add a reason column to discarded DNMs
  awk 'BEGIN {OFS="\t"} NR==1 {print \$0, "reason"} NR>1 {print \$0, "excluded_region"}' dnms_discarded_after_region_filtering.tsv > tmp && mv tmp dnms_discarded_after_region_filtering.tsv

  # Merge with previously discarded DNMs
  cat  "${discarded_dnm}" | tail -n +2 >> dnms_discarded_after_region_filtering.tsv

  """
}

process PUBLISH_FILTERED_DNM {


  publishDir "${params.outdir}/dnm/", mode: 'copy', overwrite: true


  input:
  tuple (path kept_dnm), path (discarded_dnm)

  output :
  path "dnms_kept.tsv"
  path "dnms_discarded.tsv"


  script :
  """
  cp ${kept_dnm} dnms_kept.tsv
  cp ${discarded_dnm} dnms_discarded.tsv
  """
}


process PUBLISH_ANNOTATED_DNM {


  publishDir "${params.outdir}/dnm/", mode: 'copy', overwrite: true


  input:
  path dnm

  output :
  path "dnms_annotated.tsv"


  script :
  """
  cp ${dnm} dnms_annotated.tsv
  """
}
