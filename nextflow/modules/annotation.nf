process BCFTOOLS_CSQ_FULL {

    input :
    path rates_or_dnm 
    path fasta
    path fasta_index
    path gff
    path gff_db
    val type

    output :
    path "${type}_bcftoolscsq.tsv"

    script:
    """
    # Turn rates file (tabular) into VCF
    rates_to_vcf.py rates-to-vcf $rates_or_dnm $fasta_index ${type}.vcf

    # Call bcftools consequences
    bcftools csq \
    -l \
    -v 0 \
    -g $gff \
    -f $fasta \
    -o ${type}_bcftoolscsq.vcf.gz \
   ${type}.vcf

    # Turn back VCF into rates file
    rates_to_vcf.py vcf-to-rates ${type}_bcftoolscsq.vcf.gz $rates_or_dnm $gff_db ${type}_bcftoolscsq.tsv
    """
}

process CADD {

    input :
    path rates_or_dnm 
    path cadd_file
    path cadd_file_index
    val type

    output :
    path "${type}_cadd.tsv"

    script:
    """
    # Annotate CADD
    annotate_cadd.py  $rates_or_dnm $cadd_file ${type}_cadd.tsv
    """
}

process GNOMAD {

    input :
    path rates_or_dnm 
    path gnomad_file
    path gnomad_file_index
    val type

    output :
    path "${type}_gnomad.tsv"

    script:
    """
    # Annotate CADD
    annotate_gnomad.py  $rates_or_dnm $gnomad_file ${type}_gnomad.tsv
    """
}

process CONSTRAINTS {

    input :
    path rates_or_dnm 
    path gene_full_constraints
    path gene_region_constraints
    val type


    output :
    path "${type}_constrained.tsv"

    script:
    """
    # Annotate constraints
    annotate_constraint.py  $rates_or_dnm $gene_full_constraints $gene_region_constraints ${type}_constrained.tsv
    """
}

process SHET {

    input :
    path rates_or_dnm 
    path shet
    val type

    output :
    path "${type}_shet.tsv"

    script:
    """
    # Annotate shet
    annotate_shet.py  $rates_or_dnm $shet ${type}_shet.tsv
    """
}

process DBNSFP {

    input :
    path rates_or_dnm 
    path dbnsfp
    path dbnsfp_index
    path dbnsfp_columns_to_extract
    path gffutils_db
    val type

    output :
    path "${type}_dbnsfp.tsv"

    script:
    """
    # Annotate dbnsfp
    annotate_dbnsfp.py --gff $gffutils_db $rates_or_dnm $dbnsfp $dbnsfp_columns_to_extract ${type}_dbnsfp.tsv
    """
}