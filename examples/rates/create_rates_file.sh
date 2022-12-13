#!/bin/bash
#BSUB -q long 
#BSUB -G "ddd-grp"
#BSUB -M50000
#BSUB -R "select[mem>50000] rusage[mem=50000]"
#BSUB -J create_rates
#BSUB -o /lustre/scratch119/realdata/mdt2/teams/hurles/users/ed11/DNW/out_rate_creation/full_fasta_unzip/stdout.log
#BSUB -e /lustre/scratch119/realdata/mdt2/teams/hurles/users/ed11/DNW/out_rate_creation/full_fasta_unzip/stderr.log
source /software/team29/ed11/miniconda3/etc/profile.d/conda.sh
conda activate denovowest
python \
    /lustre/scratch119/realdata/mdt2/teams/hurles/users/ed11/DNW/DeNovoWEST/denovowest/rates/create_rates_file.py \
    --config /lustre/scratch119/realdata/mdt2/teams/hurles/users/ed11/DNW/DeNovoWEST/examples/rates/rate_creation_config.yaml
conda deactivate


