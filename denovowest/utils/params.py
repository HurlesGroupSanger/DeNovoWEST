from enum import Enum

#############################
# SHARED BY SEVERAL MODULES #
#############################

# Number of adjacent nucleotides to consider when retrieving CDS
CDS_OFFSET = 50


# Bcftools consequences ordered by severity
CONSEQUENCES_SEVERITIES = {
    "intergenic": 1,
    "feature_truncation": 2,
    "feature_elongation": 2,
    "regulatory": 3,
    "TF_binding_site": 4,
    "TFBS": 4,
    "downstream": 5,
    "upstream": 5,
    "non_coding_transcript": 6,
    "non_coding": 6,
    "intron": 7,
    "NMD_transcript": 7,
    "non_coding_transcript_exon": 8,
    "5_prime_utr": 9,
    "3_prime_utr": 9,
    "coding_sequence": 10,
    "mature_miRNA": 10,
    "stop_retained": 11,
    "start_retained": 11,
    "synonymous": 11,
    "incomplete_terminal_codon": 12,
    "splice_region": 13,
    "missense": 14,
    "inframe": 14,
    "inframe_insertion": 14,
    "inframe_deletion": 14,
    "protein_altering": 14,
    "transcript_amplification": 15,
    "exon_loss": 16,
    "disruptive": 17,
    "start_lost": 18,
    "stop_lost": 18,
    "stop_gained": 18,
    "frameshift": 18,
    "splice_acceptor": 19,
    "splice_donor": 19,
    "transcript_ablation": 20,
}


##############
# RATES #
##############


# Per generation mutation rate scaling factors taken from https://github.com/vseplyarskiy/Roulette/tree/main/adding_mutation_rate
# Divided by 2 because the provided rate are per diploid genomes
# TODO : Because roulette scaling factor is cohort dependent, we should provide it as a parameter rather than hardcoded value
# ROULETTE_SCALING_FACTOR = 1.015e-7 / 2  # This is the scaling factor recommended on Roulette's github repo
ROULETTE_SCALING_FACTOR = (
    5.418997e-06 / 2
)  # This reflects the number of patients in our cohort ~1/185000 applied to the scaled roulette rates based on synonymous background
CARLSON_SCALING_FACTOR = 2.086e-9 / 2


##############
# SIMULATION #
##############

# We use meta category in the simulation (e.g. start lost are assimilated to missense)
# Part of the mapping is legacy from DNW v1 where consequences were using mixing two format
CONSEQUENCES_MAPPING = {
    "frameshift_variant": "frameshift",
    "frameshift": "frameshift",
    "inframe_insertion": "inframe",
    "inframe_deletion": "inframe",
    "inframe": "inframe",
    "missense_variant": "missense",
    "missense": "missense",
    "stop_gained": "nonsense",
    "nonsense": "nonsense",
    "splice_acceptor_variant": "splice_lof",
    "splice_donor_variant": "splice_lof",
    "splice_acceptor": "splice_lof",
    "splice_donor": "splice_lof",
    "splice_lof": "splice_lof",
    "splice_region_variant": "splice_region",
    "splice_region": "splice_region",
    "start_lost": "start_lost",
    "stop_lost": "stop_lost",
    "stop_retained": "stop_retained",
    "synonymous": "synonymous",
    "synonymous_variant": "synonymous",
}


# Default maximum number of expected mutation to test for using the simulation approach
DEFAULT_MAX_NB_MUTATIONS_SIM = 250

# Threshold to stop/skip simulation when poisson probabilities are extremely low
STOP_SKIP_SIMULATION_THRESHOLD = 10**-12


# Simulation run type
class RunType(str, Enum):
    ALL_CODING = "all-coding"
    MISSENSE = "mis"
    SYNONYMOUS = "syn"
    PTV = "ptv"
    PROTEIN_ALTERING = "protein-altering"


class ConsequenceGroups(list, Enum):
    SYNONYMOUS = ["synonymous"]
    MISSENSE = ["missense"]
    PTV = ["stop_gained", "frameshift", "splice_acceptor", "splice_donor"]
    PROTEIN_ALTERING = [
        "missense",
        "inframe_insertion",
        "inframe_deletion",
        "stop_gained",
        "frameshift",
        "splice_acceptor",
        "splice_donor",
        "start_lost",
        "stop_lost",
    ]


# Some indels are located in canonical splice site or splice region
# We estimate their mutation rate based on the indel/snv ratio observed in the 180k cohort
RATIO_SPLICE_SITE_INDEL_SNV = 0.16
RATIO_SPLICE_REGION_INDEL_SNV = 0.1
