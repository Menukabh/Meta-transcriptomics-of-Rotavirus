#!/bin/bash
#SBATCH --account=PAS3107
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=5GB
#SBATCH --output=slurm-fastp_output-%j.out

set -euo pipefail

## Define constant
CONTAINER=oras://community.wave.seqera.io/library/fastp:1.3.2--916946baf992e235

## Define variable
R1=$1
R2=$2

## Run fastp
apptainer exec "$CONTAINER" fastp -i in.R1.fq.gz -I in.R2.fq.gz \
      -o out.R1.fq.gz -O out.R2.fq.gz \
      --disable_adapter_trimming \
      --disable_quality_filtering \
      --disable_length_filtering \
      --trim_poly_g

echo 'Done with the fastp script'
date