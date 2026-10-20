# Glossary

**All-coding test**
:   The enrichment test run on all coding and splice variants. See [Run types](../method/enrichment.md#run-types).

**CADD**
:   Combined Annotation Dependent Depletion, a score of the deleteriousness of variants
    across the genome ([Kircher *et al.* 2014](https://doi.org/10.1038/ng.2892)). The PHRED-scaled
    score is added as the `score` column.

**CEP**
:   Computational effect predictor: any tool that scores how damaging a variant is
    (CADD, AlphaMissense, REVEL, popEVE…). DeNovoWEST 2 uses CEP scores as variant weights.

**Clustering test (DNN)**
:   Test of whether missense DNMs of a gene are closer to each other than expected, run
    with denovonear. See [Clustering test](../method/clustering.md).

**Consequence**
:   The predicted effect of a variant on the protein (missense, stop gained, splice donor…),
    annotated with bcftools csq.

**DNM**
:   *De novo* mutation: a variant present in a child and absent from both parents.

**Enrichment test (DNE)**
:   Test of whether the summed score of the DNMs of a gene is higher than expected. See
    [Enrichment test](../method/enrichment.md).

**Gene score**
:   Sum of the scores of the DNMs observed in a gene.

**Mutation rate model**
:   A model giving the probability that each possible substitution arises *de novo* in one
    generation: trinucleotide context ([Samocha *et al.* 2014](https://doi.org/10.1038/ng.3050)),
    [Roulette](https://github.com/vseplyarskiy/Roulette) or Carlson. See [Mutation rates](../method/rates.md).

**O/E ratio**
:   Ratio of the observed to the expected number of DNMs. Close to 1 for synonymous variants
    when the mutation rate model is well calibrated.

**PPV weights**
:   Positive predictive values used as variant weights in DeNovoWEST v1: the probability that a
    DNM of a given class and CADD score is pathogenic, estimated from the enrichment of DNMs
    in the DDD cohort. Replaced by CEP scores in v2.

**PTV**
:   Protein-truncating variant: stop gained, frameshift, splice donor and acceptor.

**Rates file**
:   The table of all possible substitutions in the analysed genes, with their mutation
    probability and annotations. See [Input files](../user-guide/inputs.md#rates-file).

**s<sub>het</sub>**
:   The selection coefficient against heterozygous protein-truncating variants in a gene
    ([Cassa *et al.* 2017](https://doi.org/10.1038/ng.3831)). Used in v1 to stratify PPV weights.
