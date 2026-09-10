#!/usr/bin/env bash
set -euo pipefail

# Usage: ./tsv_to_vcf.sh input.tsv.gz output.vcf.gz

infile="$1"
outfile="$2"

zcat "$infile" | awk -F'\t' -v OFS='\t' '
NR == 1 {
    # Find column indices by header name
    for (i = 1; i <= NF; i++) {
        col[$i] = i
    }
    if (!("chrom" in col) || !("pos" in col) || !("ref" in col) || \
        !("alt" in col) || !("prob" in col)) {
        print "ERROR: missing required column in header" > "/dev/stderr"
        exit 1
    }

    # Print VCF header
    print "##fileformat=VCFv4.3"
    print "##contig=<ID=1,length=248956422,assembly=GRCh38>"
    print "##contig=<ID=2,length=242193529,assembly=GRCh38>"
    print "##contig=<ID=3,length=198295559,assembly=GRCh38>"
    print "##contig=<ID=4,length=190214555,assembly=GRCh38>"
    print "##contig=<ID=5,length=181538259,assembly=GRCh38>"
    print "##contig=<ID=6,length=170805979,assembly=GRCh38>"
    print "##contig=<ID=7,length=159345973,assembly=GRCh38>"
    print "##contig=<ID=8,length=145138636,assembly=GRCh38>"
    print "##contig=<ID=9,length=138394717,assembly=GRCh38>"
    print "##contig=<ID=10,length=133797422,assembly=GRCh38>"
    print "##contig=<ID=11,length=135086622,assembly=GRCh38>"
    print "##contig=<ID=12,length=133275309,assembly=GRCh38>"
    print "##contig=<ID=13,length=114364328,assembly=GRCh38>"
    print "##contig=<ID=14,length=107043718,assembly=GRCh38>"
    print "##contig=<ID=15,length=101991189,assembly=GRCh38>"
    print "##contig=<ID=16,length=90338345,assembly=GRCh38>"
    print "##contig=<ID=17,length=83257441,assembly=GRCh38>"
    print "##contig=<ID=18,length=80373285,assembly=GRCh38>"
    print "##contig=<ID=19,length=58617616,assembly=GRCh38>"
    print "##contig=<ID=20,length=64444167,assembly=GRCh38>"
    print "##contig=<ID=21,length=46709983,assembly=GRCh38>"
    print "##contig=<ID=22,length=50818468,assembly=GRCh38>"
    print "##INFO=<ID=MR,Number=1,Type=Float,Description=\"Probability score from input TSV\">"
    print "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"
    next
}
{
    chrom = $(col["chrom"])
    gsub(/^chr/, "", chrom)          # strip chr prefix
    pos   = $(col["pos"])
    ref   = $(col["ref"])
    alt   = $(col["alt"])
    prob  = $(col["prob"])

    print chrom, pos, ".", ref, alt, ".", ".", "MR=" prob
}
' > "$outfile"

echo "Wrote $outfile"

# Compress and index the VCF file
bgzip ${outfile}
tabix -p vcf ${outfile}.gz