#!/bin/bash
#SBATCH --account=PAS3107
#SBATCH --time=15:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --output=slurm-bowtie_%j.output

# Bash strict settings
set -euo pipefail

# Define constant
CONTAINER=oras://community.wave.seqera.io/library/bowtie2:2.5.5--21b0835eb76ba3c0

## Define variables
echo "Started script for bowtie2 to build an index"
date

## Run Bowties2
#apptainer exec "$CONTAINER" bowtie2-build --help
apptainer exec "$CONTAINER" bowtie2-build menuka_metatrans/results/blast_NCBI_rota/rota_genomes.fna \
  menuka_metatrans/results/bowtie2_index_rota/rota

echo "Successfully ran bowtie2 to build an index"
date