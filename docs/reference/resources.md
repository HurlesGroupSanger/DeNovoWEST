# Resources

## Included in the repository

| File | Description | Used with |
|------|-------------|-----------|
| `resources/rates_model_dnn.tsv` | Trinucleotide mutation rate model of [Samocha *et al.* 2014](https://doi.org/10.1038/ng.3050), from [denovonear](https://github.com/jeremymcrae/denovonear/blob/master/denovonear/data/rates.txt). | `rates_creation.path` with `model = "kmer"` |
| `resources/MANE.GRCh38.v1.0.ensembl_genomic_mapping.tsv` | Versioned and unversioned Ensembl gene identifiers, symbols and HGNC identifiers of MANE v1.0 genes. | Converting gene identifiers in your own files |
| `examples/` | Small inputs and expected outputs for each step, on the *OR4F5* gene. | Tests and the [quick start](../getting-started/quickstart.md) |

## External resources

| Resource | Use | Source |
|----------|-----|--------|
| GRCh38 reference genome | `genome_fasta` | e.g. [GRCh38 analysis set](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/) |
| GENCODE / MANE gene annotation | `gff` | [GENCODE](https://www.gencodegenes.org/human/), [MANE](https://www.ncbi.nlm.nih.gov/refseq/MANE/) |
| Roulette mutation rates | `rates_creation.path`, `clustering.rates` | [Roulette](https://github.com/vseplyarskiy/Roulette) |
| CADD scores | `annotation.cadd_file` | [CADD downloads](https://cadd.gs.washington.edu/download) |
| dbNSFP | `annotation.dbnsfp` | [dbNSFP](https://www.dbnsfp.org/) |
| popEVE, EVE, ESM1v | `annotation.vcf` | [EVE / popEVE](https://evemodel.org/) |
| AlphaFold human proteome | `clustering.protein_structures` | [AlphaFold DB](https://alphafold.ebi.ac.uk/download) |
| gnomAD constraint metrics | `oe_ratios.gnomad_constraints` | [gnomAD downloads](https://gnomad.broadinstitute.org/downloads) |
| Gene2Phenotype panels | `denovowest report --g2p` | [Gene2Phenotype](https://www.ebi.ac.uk/gene2phenotype/download) |

## Data from the 2020 publication

The rates files used in the 2020 publication (DeNovoWEST v1, keyed by gene symbol) are
available on the Sanger FTP:

- GRCh37: <ftp://ftp.sanger.ac.uk/pub/project/ddd/rates>
- GRCh38 (same trinucleotide model): <ftp://ftp.sanger.ac.uk/pub/project/ddd/B38_rates/B38_merged_rates.gz>

They are not compatible with DeNovoWEST 2, which uses Ensembl gene identifiers and
bcftools consequences: build a new rates file with the workflow instead.

The DNMs, weights and results of the publication are in the `input/` directory of the
repository, described in `input/README.md`.
