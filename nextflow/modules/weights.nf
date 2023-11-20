process GET_OBSERVED_COUNTS {

    input :
    path dnm_file 


    output :
    path "observed_counts.tsv"

    script :
    """
    create_count_tables.py get-observed-counts ${dnm_file}  observed_counts.tsv 
    """
}

process GET_EXPECTED_COUNTS {

    input :
    path rates_file 
    val n_males
    val n_females

    output :
    path "expected_counts.tsv"

    script :
    """
    create_count_tables.py get-expected-counts ${rates_file} ${n_males} ${n_females} expected_counts.tsv
    """
}

process MERGE_EXPECTED {
    input : 
    path expected_counts, stageAs : "expected_counts_*.tsv"

    output :
    path "merged_expected_counts.tsv"

    script :
    """
    #!/usr/bin/env python
    import sys
    import pandas as pd
    expected_counts_files = "${expected_counts}".split(" ")
    list_df = [pd.read_csv(filename, sep="\t") for filename in expected_counts_files]

    # Sum up the expected values on each gene
    s = pd.Series()
    for df in list_df :
        if s.empty :
            s = df.exp
        else :
            s = s + df.exp

    merged_df = list_df[0]
    merged_df.exp = s

    merged_df.to_csv("merged_expected_counts.tsv", sep="\t", index=False)

    """
}


process MERGE_COUNTS {

    publishDir "${params.outdir}/weights/", mode: 'copy', overwrite:true, include: 'merged_counts.tsv'
    publishDir "${params.outdir}/weights/", mode: 'copy', overwrite:true, include : 'enrichment_plots/*.png',  exlude : 'enrichment_plots'


    input:
    path expected_counts
    path observed_counts

    output:
    path "merged_counts.tsv", emit : merged_counts_ch 
    path "enrichment_plots"

    script :
    """
    create_count_tables.py merge-expected-observed ${expected_counts} ${observed_counts} merged_counts.tsv
    """
}


process LOESS {

    publishDir "${params.outdir}/weights/", mode: 'copy', overwrite: true, include : 'weights.tsv'
    publishDir "${params.outdir}/weights/", mode: 'copy', overwrite: true, include : 'enrichment_plots_loess/*.png', exlude : 'enrichment_plots_loess'


    conda "/software/team29/ed11/miniconda3/envs/rpy2"

    input:
    path merged_counts

    output :
    path "weights.tsv", emit : weights_ch
    path "enrichment_plots_loess"


    script :
    """
    fit_loess.py ${merged_counts} weights.tsv --outdir enrichment_plots_loess
    """

}