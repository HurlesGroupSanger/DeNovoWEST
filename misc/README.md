
## gff2gff
This perl script can be found in the [bcftools github repository](https://github.com/samtools/bcftools/blob/develop/misc/gff2gff). It was used to create a MANE gff file compatible with bcftoolscsq.
    
    zcat MANE.GRCh38.v1.0.ensembl_genomic.gff.gz | ./gff2gff | gzip -c > MANE.GRCh38.v1.0.ensembl_genomic.bcftoolscsq.gff.gz
