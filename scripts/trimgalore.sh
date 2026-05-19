#!/bin/bash
#SBATCH --account=PAS3107
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=24
#SBATCH --output=slurm-trimgalore-%j.out

set -euo pipefail

## Define constants: container of trim-galore
CONTAINER=oras://community.wave.seqera.io/library/trim-galore:0.6.10--bc38c9238980c80e

## Define placeholder variable
outdir=$1
R1=$2
R2=$3

## Run trimgalore: Run fastqc afgter runnin trimgalore
#apptainer exec $CONTAINER trim_galore --help
apptainer exec $CONTAINER trim_galore \
    --paired \
    --fastqc \
    --output_dir "$outdir" \
    --2colour 20 \
    --max_n 0 \
    --cores 24 \
    "$R1" "$R2"


echo "Done with Script"
date

