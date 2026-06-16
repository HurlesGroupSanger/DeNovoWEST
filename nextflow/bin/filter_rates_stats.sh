#!/bin/bash


zcat "$1" | awk -v key="$2" '
BEGIN { OFS="\t" }
{
    gene = $1
    pos  = $3
    # New position block: reset per-position dedup set (stays tiny, ~1-3 entries)
    if (pos != prev_pos) {
        delete pos_genes
        prev_pos = pos
    }
    # First occurrence of this gene at this position
    if (!(gene in pos_genes)) {
        pos_genes[gene] = 1
        count[gene]++
    }
}
END {
    print "{"
    print "  \"" key "\": {"
    n = 0
    for (gene in count) n++
    i = 0
    for (gene in count) {
        i++
        comma = (i < n) ? "," : ""
        print "    \"" gene "\": " count[gene] comma
    }
    print "  }"
    print "}"
}
'
