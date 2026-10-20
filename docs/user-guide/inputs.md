# Input files

All examples on this page come from the `examples/` directory of the repository.

!!! info "Conventions"
    - **Genome build**: DeNovoWEST 2 is developed and tested on **GRCh38**.
    - **Coordinates**: 1-based, as in VCF files, with indels left-anchored
      (e.g. `TC` → `T`). BED files are 0-based, as usual.
    - **Chromosome names**: the DNM file, the rates file, the excluded regions and the
      genome FASTA must all use the same naming as the GFF (`chr1` or `1`).
    - **Gene identifiers**: the `ID` of the gene in the GFF, version included
      (e.g. `ENSG00000186092.7`).

## DNM file

Tab-separated text file with a header, one row per DNM (and per gene, if a DNM overlaps
several genes).

| Column | Required | Description |
|--------|----------|-------------|
| `gene_id` | yes, **first column** | Gene identifier, as in the GFF. |
| `chrom` | yes, **second column** | Chromosome. |
| `pos` | yes, **third column** | Position (1-based). |
| `ref` | yes | Reference allele. |
| `alt` | yes | Alternative allele. |
| `consequence` | for the tests | Added by the workflow when `annotate_dnm` is set (see [Consequences](#consequences)). |
| score columns | for the enrichment test | Added by the workflow from the [annotation sources](#annotation-sources). |
| any other column | no | Kept as is (e.g. `gene_symbol`, sample identifier). |

```
gene_id	chrom	pos	ref	alt	varid	gene_symbol
ENSG00000186092.7	chr1	65518	A	T	var_1	OR4F5
ENSG00000186092.7	chr1	65549	G	T	var_2	OR4F5
```

If your DNMs are not yet assigned to genes, `misc/map_dnm_to_genes.py` assigns each DNM
to the genes of a GFF that it overlaps.

### Filtering

Before annotation, the workflow removes DNMs that:

- have no `gene_id`;
- are in a gene that is not in the gene list (`not_in_gene_list`);
- are on chromosome X or Y, if `exclude_sex_chromosomes = true` (`allosome`);
- do not overlap the coding sequence of their gene, extended by 8 bp on each side to
  include splice sites (`not_in_cds`, or `gene_not_in_gff`);
- start in an excluded region (`excluded_region`).

Removed DNMs are written to `dnm/dnms_discarded.tsv` with the reason in the `reason`
column (except those without `gene_id`).

## Rates file

Built by the workflow when `create_rates = true`. To reuse one, pass it with `rates`.
See [Mutation rates](../method/rates.md) for how it is built.

One row per possible single-nucleotide substitution in the coding sequence of each gene,
extended by 8 bp on each side.

| Column | Description |
|--------|-------------|
| `gene_id`, `chrom`, `pos` | **First three columns**, in this order. |
| `ref`, `alt` | Alleles. |
| `prob` | Probability that this substitution arises *de novo* on one haploid genome in one generation. |
| `gene_symbol` | Gene name, if present in the GFF (`gene_name` attribute). |
| `consequence`, `full_consequence_bcftoolscsq` | Added by [consequence annotation](#consequences). |
| score columns | Added by the [annotation sources](#annotation-sources). |

```
gene_id	chrom	pos	ref	alt	prob	consequence	full_consequence_bcftoolscsq	raw	score	transcript_in_dbnsfp	AlphaMissense_score
ENSG00000117139.19	chr1	202729039	C	A	1.164306866667408e-09	missense		2.672775	19.74	True	0.3881
ENSG00000117139.19	chr1	202729039	C	G	2.386829123982257e-09	missense		2.680624	19.77	True	0.3881
```

A rates file given with `rates` must be:

- sorted by `chrom`, `pos` and `alt`;
- compressed with `bgzip` and indexed with `tabix -S 1 -s 2 -b 3 -e 3`, with the index
  next to it (`<rates>.tbi`).

## GFF

A GFF3 file describing the genes to analyse, in the format expected by
[bcftools csq](https://samtools.github.io/bcftools/bcftools.html#csq): `gene` features,
transcripts as children of genes, and `CDS` features as children of transcripts.

Rates are built from **every transcript** of each gene in the GFF, so use a GFF with
one transcript per gene, typically the MANE Select transcripts.

The following attributes are used:

| Feature | Attribute | Use |
|---------|-----------|-----|
| gene | `ID` | Gene identifier (`gene_id` everywhere). |
| gene | `gene_name` (or `Name`) | Gene symbol. |
| gene | `gene_id` | Matching of VCF annotations by gene. |
| transcript | `transcript_id` | Matching of dbNSFP scores by transcript. |

```
chr1	HAVANA	gene	65419	71585	.	+	.	ID=ENSG00000186092.7;gene_id=ENSG00000186092.7;gene_type=protein_coding;gene_name=OR4F5;biotype=protein_coding;Name=OR4F5
```

GENCODE and Ensembl GFF3 files can be converted to the bcftools csq format with
`misc/gff2gff`, from the bcftools repository:

```bash
zcat gencode.v49.basic.annotation.gff3.gz | misc/gff2gff | gzip -c > gencode.v49.bcftoolscsq.gff3.gz
```

## Genome FASTA

The reference genome, with its index (`samtools faidx`) next to it as
`<genome_fasta>.fai`. Chromosome names must match the GFF.

## Gene list

Plain text file, one gene identifier per line, no header:

```
ENSG00000186092.7
```

Identifiers must match the GFF `ID`, version included.

## Excluded regions

A BED file of regions to exclude, e.g. simple repeats, segmental duplications and the
ENCODE blacklist. Variants in these regions are removed from the rates file and the DNM
file, so that both are restricted to the same part of the genome.

```
chr1	0	10000
chr1	207666	257666
```

## Sperm positive selection adjustment

A JSON file of multipliers applied to the mutation rates of genes or sites under positive
selection in the male germline. See
[Mutation rates](../method/rates.md#adjusting-for-positive-selection-in-sperm).

```json
{
  "genes": {
    "ENSG00000117139": { "missense": 4.88, "stop_gained": 27.92 }
  },
  "sites": {
    "chr12_112477719_A_G": {
      "chrom": "chr12", "pos": 112477719, "ref": "A", "alt": "G", "multiplier": 420.3
    }
  }
}
```

- `genes`: for each gene (with or without version), a multiplier per consequence. The
  consequence must match the `consequence` column of the rates file, so the rates file
  must be annotated.
- `sites`: a multiplier for a specific substitution. The key is free text; `chrom` must
  use the same naming as the rates file.

Both sections are optional, but at least one is required. A site in an adjusted gene gets
both multipliers.

## Annotation sources

All sources must be compressed with `bgzip` and indexed with `tabix` (`.tbi` next to the
file), and **use chromosome names without the `chr` prefix**. Annotations are matched on
chromosome, position, reference and alternative allele.

If an annotation column has the same name as an existing column, it is renamed
`<column>_<source>` (e.g. `score_cadd`). Missing values are written as `.`.

### Consequences

With `annotation.annotate_bcftoolscsq = true` (default), each variant is annotated with
`bcftools csq` using the GFF. Two columns are added:

- `consequence`: the most severe consequence of the variant on its own gene
  (e.g. `missense`, `stop_gained`, `splice_donor`, `synonymous`);
- `full_consequence_bcftoolscsq`: the full combined consequence when it differs
  (e.g. `missense&splice_region`).

Consequences use the bcftools csq vocabulary. If you provide your own `consequence`
column, it must use the same terms.

### CADD

The CADD whole-genome SNV file, e.g. `whole_genome_SNVs.tsv.gz` (CADD v1.7, GRCh38),
with six columns: chromosome, position, reference, alternative, raw score, PHRED score.

```
#chrom	pos	ref	alt	raw	score
1	65419	C	A	1.159795	11.89
```

Adds the columns `raw` and `score`.

### dbNSFP

The genome-wide [dbNSFP](https://www.dbnsfp.org/) file for GRCh38 (tested with versions
5.2a and 5.3.1a). The header must include `Ensembl_geneid` and `Ensembl_transcriptid`.

When dbNSFP has several values for a variant (one per transcript), the maximum over the
transcripts present in the GFF is kept, or else the maximum over the gene's transcripts.
Adds the requested columns and `transcript_in_dbnsfp`.

### Custom TSV

Any tab-separated file whose first four columns are chromosome, position, reference and
alternative allele, with a header line (starting with `#` or not). Adds the requested
columns, or all columns after the fourth.

### VCF

Any VCF file with annotations in the INFO column, e.g. popEVE or COSMIC. Adds the
requested INFO fields, or all of them. Only the first alternative allele of each record
is used.

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	65565	.	A	G	.	.	gene=OR4F5;popEVE=-4.754;EVE=10.471;ESM1v=-4.841
```

## O/E ratio inputs

| Parameter | Format |
|-----------|--------|
| `oe_ratios.gene_categories` | TSV with a header, `gene_id` in the first column (as in the rates file) and a `category` column. Genes absent from the file are grouped as `not_in_user_gene_list`. |
| `oe_ratios.gnomad_constraints` | The gnomAD v4.1 `constraint_metrics.tsv` file. |
