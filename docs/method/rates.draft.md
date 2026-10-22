# Mutation rates

The rates file is the null model of DeNovoWEST: for every possible single-nucleotide
substitution in each gene, the probability that it arises *de novo* in one haploid genome
in one generation. The enrichment test multiplies these probabilities by 2 × (number of
individuals) to get expected numbers of DNMs in the cohort.

## Which sites are included

For each gene of the GFF, the workflow takes the coding sequence (CDS) of each of its
transcripts, extends each exon by **8 bp** on each side to capture splice sites and most of
the splice region, and lists the three possible substitutions at every position.

Use a GFF with one transcript per gene (e.g. MANE Select): with several transcripts, the
union of their CDS is used.

## Mutation rate models

Set the model with `rates_creation.model` and its data with `rates_creation.path`.

### `kmer` (default)

A sequence context model: the rate of a substitution depends on the bases around it.
`path` is a table with the columns `from`, `to` and `mu_snp`:

```
from	to	mu_snp
AAA	ACA	1.6681119889870802e-09
AAA	AGA	2.9280578256688802e-09
```

The repository includes `resources/rates_model_dnn.tsv`, the trinucleotide model of
[Samocha *et al.* 2014](https://doi.org/10.1038/ng.3050) as distributed with
denovonear. Larger contexts (5-mers, 7-mers…) work too: the context size is taken
from the length of `from`. The genome FASTA is used to read the context, on the strand
of the transcript.

### `roulette` and `carlson`

Pre-computed, genome-wide rates from
[Roulette](https://github.com/vseplyarskiy/Roulette), distributed as one VCF per
chromosome. `path` is the directory holding these VCFs. Each file name must start with
the chromosome followed by `_` and end with `all.vcf.gz` (e.g. `1_all.vcf.gz`); files must
be indexed and use chromosome names without `chr`.

- `roulette` uses the Roulette rate (`INFO/MR`);
- `carlson` uses the rate of the [Carlson *et al.* 2018](https://doi.org/10.1038/s41467-018-05936-5)
  model provided in the same files (`INFO/MC`).

Roulette files do not cover chromosomes X and Y, so genes on sex chromosomes get no rates
with these models.

#### Scaling factor

Roulette and Carlson rates are relative: they must be multiplied by a scaling factor,
`rates_creation.scaling_factor`, to become the per-generation, **per haploid genome**
probabilities that DeNovoWEST expects. There are two ways to choose it.

**Default factors.** The [Roulette authors](https://github.com/vseplyarskiy/Roulette/tree/main/adding_mutation_rate)
give factors that convert the raw rates into per-generation rates per **diploid**
genome. Divide them by 2:

| Model | Roulette's factor (diploid) | `scaling_factor` (haploid) |
|-------|-----------------------------|----------------------------|
| `roulette` (`INFO/MR`) | 1.015 × 10<sup>-7</sup> | `5.075e-8` |
| `carlson` (`INFO/MC`) | 2.086 × 10<sup>-9</sup> | `1.043e-9` |

**Factor calibrated on your cohort (recommended).** The number of DNMs a cohort
yields depends on parental age, sequencing and variant calling, so a default factor can
over- or under-estimate the expected number of DNMs. Roulette's authors recommend
calibrating the rates on putatively neutral variants: for exome data, the synonymous DNMs.
The factor is then chosen so that the expected number of synonymous DNMs equals the
observed number:

```
scaling_factor = n_syn_observed / (2 × (n_males + n_females) × Σ MR over synonymous sites)
```

where the sum runs over all possible synonymous substitutions of the analysed genes.
Roulette's `add_scaled_rates.py --denovo 1` computes the same ratio without the
2 × (n<sub>males</sub> + n<sub>females</sub>) term, which DeNovoWEST applies itself.

In practice:

1. build and annotate a rates file with `scaling_factor = 1`, so that `prob` holds the raw rates;
2. sum `prob` over the rows whose `consequence` is `synonymous`, and count the synonymous
   SNVs in `dnm/dnms_annotated.tsv`, after the same filtering (gene list, excluded regions);
3. compute the factor with the formula above and rebuild the rates file with it.

Then run the [O/E ratio QC](../user-guide/configuration.md#observed-expected-ratios-oe_ratios):
the observed/expected ratio of synonymous DNMs should be close to 1.

!!! note
    Calibrating on synonymous DNMs absorbs any genome-wide excess or deficit of DNMs in
    your cohort. Genes are then tested for an excess relative to the cohort's own
    background rate.

## Filtering excluded regions

With `excluded_regions`, substitutions in the regions of the BED file are removed from the
rates file, and DNMs in the same regions are removed from the DNM file. The number of
positions removed per gene is written to `rates/stats/filtered_rates_stats.json`<!--, and the
[report](../user-guide/report.md) flags genes that lost more than 10% of their sites-->.

## Adjusting for positive selection in sperm

A few mutations give spermatogonial stem cells a growth advantage, so that they become
more frequent in sperm as men age. Genes carrying these mutations (e.g. *FGFR2*, *FGFR3*,
*HRAS*, *PTPN11*) show more DNMs than predicted by their mutation rate, which the
enrichment test would otherwise attribute to disease.

With `sperm_positive_selection_adjustment`, the workflow multiplies the rates of these
genes (per consequence) or sites by the factors given in a JSON file
(see [Input files](../user-guide/inputs.md#sperm-positive-selection-adjustment)). The total
rate of each adjusted gene before and after adjustment is written to
`rates/stats/adjusted_rates.json`<!--, and the [report](../user-guide/report.md) flags these
genes-->.

## VCF version of the rates

After annotation, filtering and adjustment, the workflow also writes the rates as a VCF
(`rates/rates_adjusted.vcf.gz`, rate in `INFO/MR`). Pass it to `clustering.rates` so that
the clustering test uses the same mutation rates as the enrichment test.

## Statistics

`rates/stats/rates_stats.json` summarises the rates file: number of sites and total rate
per gene and consequence, and number of sites without a score for each score column.<!-- The
[report](../user-guide/report.md) uses it to flag genes with many missing scores.-->
