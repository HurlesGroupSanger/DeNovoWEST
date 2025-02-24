process SPLIT_GENE_LIST_CLUSTERING {

    input :
    path dnm
    path gene_list 
    val split_step

    output :
    path 'chunk_*'

    script :
    """
    grep -F -f <(cut -f1 $dnm) $gene_list > gene_list_dnm.tsv
    split -l $split_step gene_list_dnm.tsv chunk_
    """

}





process SPLIT_DNM {

  input : 
  path gene_list
  path dnm

  
  output : 
  path "${gene_list}_dnm.tsv"

  script : 
  """
  # Extract header
  cat $dnm | head -n 1 > ${gene_list}_dnm.tsv
  # Keep only the rows where gene is in gene_list
  cat $dnm | awk 'FNR==NR {values[\$1]; next} \$1 in values' $gene_list - >> ${gene_list}_dnm.tsv
  """

}


process PREPARE_DNM_DENOVONEAR {


	input :
    path dnm
	path gtf

    output :
    path "dnn_denovonear.tsv"

    script :
    """
    prepare_dnms.py $dnm \
		$gtf \
		dnn_denovonear.tsv
    """
}




process DENOVONEAR_LINEAR {

    label "process_long"

	input :
    path dnm
	path gtf
    path fasta

    output :
    path "linear_clustering_results.tsv", emit :results

    script :
    """
    denovonear cluster --in $dnm --gencode $gtf --fasta $fasta --out linear_clustering_results.tsv
    """
}

process DENOVONEAR_3D {

    label "process_long"

	input :
    path dnm
	path gtf
    path fasta
	path protein_structures

    output :
    path "3d_clustering_results.tsv", emit :results

    script :
    """
    denovonear cluster-structure \
		--in $dnm \
		--structures $protein_structures \
		--gencode $gtf \
		--fasta $fasta \
        --genome-build grch38 \
		--out 3d_clustering_results.tsv
    """
}


process MERGE_CLUSTERING {
    
	input :
    path clustering_results, stageAs : "clustering_results_*.tsv"
    val runtype

    output :
    path "clustering_results.tsv"

    script :
    """
    # Concatenates the tsv files and keep only the first header
    awk 'FNR==1 && NR!=1{next;}{print}' $clustering_results > clustering_results.tsv

    """
}


process ADD_GENE_ID {

	publishDir "${params.outdir}/clustering/$runtype", mode: 'copy', overwrite: true
    
	input :
    path clustering_results, stageAs : "clustering_results_symbol.tsv"
    path gff
    val runtype

    output :
    path "clustering_results.tsv"

    script :
    """
    gene_symbol_to_id.py $clustering_results $gff clustering_results.tsv

    """
}



process COMBINE_CLUSTERING {

    publishDir "${params.outdir}/clustering/combined", mode: 'copy', overwrite: true

    input :
    path linear_clustering_results, stageAs : "linear_clustering_results.tsv"
    path threed_clustering_results, stageAs : "3d_clustering_results.tsv"

    output :
    path "clustering_results.tsv"

    script :
    """
    combine_test_results.py $linear_clustering_results $threed_clustering_results clustering_results.tsv
    """
}

process COMBINE_DNE_DNN {

    publishDir "${params.outdir}/combined_tests", mode: 'copy', overwrite: true

    input :
    path enrichment_results_ns, stageAs : "enrichment_results_ns.tsv"
    path enrichment_results_mis, stageAs : "enrichment_results_mis.tsv"
    path clustering_results_mis

    output :
    path "combined_tests_results.tsv"

    script :
    """
    combine_dne_dnn.py $enrichment_results_ns $enrichment_results_mis $clustering_results_mis combined_tests_results.tsv
    """
}