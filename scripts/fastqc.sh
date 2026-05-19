#!/bin/bash
#SBATCH --account=PAS3107
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=5GB
#SBATCH --output=slurm-fastqc_output-%j.out

set -euo pipefail

## Container for fastqc
Container=oras://community.wave.seqera.io/library/fastqc:0.12.1--104d26ddd9519960

## Define variable
outdir=$1
fastq=$2

## Run Fastqc
#apptainer exec $Container fastqc --help
apptainer exec $Container fastqc --outdir "$outdir" \
    --threads 4 \
    "$fastq"

echo 'Done with script'
date