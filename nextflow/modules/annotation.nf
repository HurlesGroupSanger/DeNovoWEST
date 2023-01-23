process RATES_TO_VCF {

    input :
    path rate_file 
    path fasta_index

    output :
    path "all_positions.vcf"

    script :
    """
    rates_to_vcf.py $rate_file $fasta_index all_positions.vcf
    """
}

process BCFTOOLS_CSQ {

    input :
    path vcf 
    path gff
    path fasta

    output :
    path "all_positions.vcf"

    script :
    """
    bcftools csq \
        -l \
        -v 0 \
        -g $gff \
        -f $fasta \
        -o all_positions_bcftoolscsq.vcf.gz \
        $vcf
    """
}