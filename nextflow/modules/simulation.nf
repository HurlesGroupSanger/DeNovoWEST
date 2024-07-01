process SIMULATION {

    
	label "process_long"

	input :
    tuple path(dnm), path(rates)
    val column
	val nmales
	val nfemales
    val export_weighted_dnmrates

    output :
    path "enrichment_results.tsv", emit :results
    path("weighted_DNM.tsv"), optional: true, emit : dnm
    path("weighted_rates.tsv"), optional: true, emit : rates

    beforeScript "[ -v NF_TEST ] && export PYTHONPATH=$baseDir/../../../../../denovowest/simulation/;"

    script :
    """
    if $export_weighted_dnmrates; then
        echo "Running DNW simulation and exporting the weight annotated files"
        simulation.py $dnm $rates $column --nmales $nmales --nfemales $nfemales --export_weighted_dnmrates
    else
        echo "Running DNW simulation"
        simulation.py $dnm $rates $column --nmales $nmales --nfemales $nfemales 
    fi
    """
}

process MERGE_SIMULATION {

	publishDir "${params.outdir}/simulation", mode: 'copy', overwrite: true
    
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