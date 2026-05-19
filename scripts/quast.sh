#!/bin/bash
#SBATCH --account=PAS3107
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=16GB
#SBATCH --output=slurm-quast_output-%j.out

set -euo pipefail

# Define constant
container=oras://community.wave.seqera.io/library/quast:5.3.0--bfd4c029fde7e696

# Define variable
fasta=$1
outdir=$2

## Run quast
#apptainer exec $container quast.py --help
apptainer exec $container quast.py \
    "$fasta" \
    -o "$outdir"

echo "Done with Quast"
date