process SIMULATION {

    
	memory "100G"

	input :
    tuple path(dnm), path(rates), path(weights) 
	val nmales
	val nfemales

    output :
    path "enrichment_results.tsv"

    script :
    """
    simulation.py $dnm $rates $weights --nmales $nmales --nfemales $nfemales
    """
}

process MERGE_SIMULATION {

	publishDir "${params.outdir}/simulation", mode: 'copy', overwrite: true
    
	memory "100G"

	input :
    path enrichment_results, stageAs : "enrichment_results_*.tsv"

    output :
    path "enrichment_results.tsv"

    script :
    """
    # Concatenates the tsv files and keep only the first header
    awk 'FNR==1 && NR!=1{next;}{print}' $enrichment_results > enrichment_results.tsv

    """
}