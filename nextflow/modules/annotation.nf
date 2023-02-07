process RATES_TO_VCF {

    input :
    path rate_file 
    path fasta_index

    output :
    path "mutation_rates.vcf"

    script :
    """
    rates_to_vcf.py rates-to-vcf $rate_file $fasta_index mutation_rates.vcf
    """
}

process BCFTOOLS_CSQ {

    input :
    path vcf 
    path gff
    path fasta

    output :
    path "mutation_rates_bcftoolscsq.vcf.gz"

    script :
    """
    bcftools csq \
        -l \
        -v 0 \
        -g $gff \
        -f $fasta \
        -o mutation_rates_bcftoolscsq.vcf.gz \
        $vcf
    """
}

process VCF_TO_RATES {

    input :
    path vcf 
    path rates

    output :
    path "mutation_rates_bcftoolscsq.tsv"

    script :
    """
    rates_to_vcf.py vcf-to-rates $vcf $rates mutation_rates_bcftoolscsq.tsv
    """
}

process BCFTOOLS_CSQ_FULL {

    input :
    path rate_file 
    path fasta
    path fasta_index
    path gff

    output :
    path "mutation_rates_bcftoolscsq.tsv"

    script:
    """
    # Turn rates file (tabular) into VCF
    rates_to_vcf.py rates-to-vcf $rate_file $fasta_index mutation_rates.vcf

    # Call bcftools consequences
    bcftools csq \
    -l \
    -v 0 \
    -g $gff \
    -f $fasta \
    -o mutation_rates_bcftoolscsq.vcf.gz \
    mutation_rates.vcf

    # Turn back VCF into rates file
    rates_to_vcf.py vcf-to-rates mutation_rates_bcftoolscsq.vcf.gz $rate_file mutation_rates_bcftoolscsq.tsv
    """
}

process CADD {

    input :
    path rate_file 
    path cadd_file
    path cadd_file_index

    output :
    path "mutation_rates_cadd.tsv"

    script:
    """
    # Annotate CADD
    annotate_cadd.py  $rate_file $cadd_file mutation_rates_cadd.tsv
    """
}

process GNOMAD {

    input :
    path rate_file 
    path gnomad_file
    path gnomad_file_index

    output :
    path "mutation_rates_gnomad.tsv"

    script:
    """
    # Annotate CADD
    annotate_gnomad.py  $rate_file $gnomad_file mutation_rates_gnomad.tsv
    """
}