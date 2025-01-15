process SIMULATION {

    
	label "process_long"

	input :
    tuple path(dnm), path(rates)
    val column
	val nmales
	val nfemales
    val runtype

    output :
    path "enrichment_results.tsv", emit :results
    path("weighted_DNM.tsv"), optional: true, emit : dnm
    path("weighted_rates.tsv"), optional: true, emit : rates

    beforeScript "[ -v NF_TEST ] && export PYTHONPATH=$baseDir/../../../../../denovowest/simulation/;"

    script :
    """
    simulation.py $dnm $rates $column --nmales $nmales --nfemales $nfemales --runtype $runtype
    """
}

process MERGE_SIMULATION {

	publishDir "${params.outdir}/simulation/$runtype", mode: 'copy', overwrite: true
    
	input :
    path enrichment_results, stageAs : "enrichment_results_*.tsv"
    val runtype

    output :
    path "enrichment_results.tsv"

    script :
    """
    # Concatenates the tsv files and keep only the first header
    awk 'FNR==1 && NR!=1{next;}{print}' $enrichment_results > enrichment_results.tsv

    """
}