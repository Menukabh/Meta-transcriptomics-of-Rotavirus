#!/bin/bash
#SBATCH --account=PAS3107
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=16GB
#SBATCH --output=slurm-spades_output-%j.out

## Bash strict settings
set -euo pipefail

## Load container and load the either virus_odb10 database
CONTAINER=oras://community.wave.seqera.io/library/spades:4.2.0--3313822b80929818

## Arguments
R1=$1
R2=$2
singleton=$3
outdir=$4

# Run metaspades
#apptainer exec "$CONTAINER" metaspades.py --help
apptainer exec "$CONTAINER" metaspades.py \
    -1 "$R1" \
    -2 "$R2" \
    -s "$singleton" \
    -o "$outdir"

echo "Done with the metaspades"
date