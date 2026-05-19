#!/bin/bash
#SBATCH --account=PAS3107
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=5GB
#SBATCH --output=slurm-multiqc_output-%j.out

set -euo pipefail

## Define constant
Container=oras://community.wave.seqera.io/library/multiqc:1.33--e3576ddf588fa00d

## Define variable
fastqc=$1
outdir=$2

# Run multiqc
## apptainer exec $Container multiqc --help
apptainer exec $Container multiqc "$fastqc" \
    --outdir "$outdir"

echo 'Done with the script'
date