#!/bin/bash
#SBATCH --account=PAS3107
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=16GB
#SBATCH --output=slurm-seqkit_output-%j.out

## Bash strict settings
set -euo pipefail

# Define constant
container=oras://community.wave.seqera.io/library/seqkit:2.13.0--205358a3675c7775

# Define variable
contigs=$1
outfile=$2

## Run seqkit to filter contigs less than 300bp
#apptainer exec $container seqkit seq --help
apptainer exec $container seqkit seq \
--min-len 300 "$contigs" \
> "$outfile"

echo "Done with seqkit"
date