process SIMULATION {

    
	label "process_long"

    publishDir "${params.outdir}/simulation/$runtype/batches", mode: 'symlink', overwrite: true


	input :
    tuple path(dnm), path(rates), val(id)
    val column
	val nmales
	val nfemales
    val runtype
    val nsim
    val impute_missing_scores

    output :
    path "${id}_enrichment_results.tsv", emit :results
    path "${id}_simulation_logs.json", emit :logs



    beforeScript "[ -v NF_TEST ] && export PYTHONPATH=$baseDir/../../../../../denovowest/simulation/;"

    script :
    """
    if [  "$impute_missing_scores" = "true"  ]; then
        simulation.py $dnm $rates $column --nmales $nmales --nfemales $nfemales --runtype $runtype --nsim $nsim --impute-missing-scores --outfile ${id}_enrichment_results.tsv --debug
    else
        simulation.py $dnm $rates $column --nmales $nmales --nfemales $nfemales --runtype $runtype --nsim $nsim --outfile ${id}_enrichment_results.tsv --debug
    fi

    mv simulation_logs.json ${id}_simulation_logs.json

    """
}

process MERGE_SIMULATION {

	publishDir "${params.outdir}/simulation/$runtype", mode: 'copy', overwrite: true
    
	input :
    path enrichment_results
    path simulation_logs
    val runtype

    output :
    path "enrichment_results.tsv"
    path "simulation_logs.json"

    script :
    """
    # Concatenates the tsv files and keep only the first header
    awk 'FNR==1 && NR!=1{next;}{print}' $enrichment_results > enrichment_results.tsv

    jq -s 'add' *.json > simulation_logs.json

    """
}