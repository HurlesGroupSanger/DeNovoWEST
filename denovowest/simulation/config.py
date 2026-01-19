# config.py
from dataclasses import dataclass

from denovowest.utils.params import RUNTYPE_ALL_CODING


@dataclass
class Config:
    # Cohort
    nmales: int = 0  # Required
    nfemales: int = 0  # Required
    gene_list: str = ""

    # Variant processing
    impute_missing_scores: bool = False
    runtype: str = RUNTYPE_ALL_CODING

    # Inframe variant scoring
    inframe_missense_ratio: float = 0.03
    inframe_score_quantile: float = 0.5  # 0.5 = median

    # Frameshift variant scoring
    frameshift_nonsense_ratio: float = 1.3
    frameshift_score_quantile: float = 0.6  # 0.6 = 60th percentile

    # Simulation control
    nsim: int = 10**7
    pvalcap: float = 0.01
    jobs: int = 1

    # Output / debugging
    debug: bool = False
    outdir: str = "./"
    outfile: str = "enrichment_results.tsv"
