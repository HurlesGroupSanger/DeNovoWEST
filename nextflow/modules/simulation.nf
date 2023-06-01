process SIMULATION {

	publishDir "${params.outdir}/simulation", mode: 'copy'
    
	memory "100G"

	input :
    path dnm 
    path rates
	path weights
	val nmales
	val nfemales

    output :
    path "enrichment_results.tsv"

    script :
    """
    simulation.py $dnm $rates $weights --nmales $nmales --nfemales $nfemales
    """
}