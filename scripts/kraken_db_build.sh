#!/bin/bash
#SBATCH --account=PAS3107
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --output=slurm-krakendb-%j.out

# Define container and db
kraken_cont=oras://community.wave.seqera.io/library/kraken2:2.17.1--1738c34504f3fb18
DB=/fs/scratch/PAS0471/menuka/kraken_host_db_metatrans

# Build the db after adding pig and cell genome sequence
apptainer exec "$kraken_cont" kraken2-build \
  --build \
  --threads 8 \
  --db "$DB"

echo "Build the custom db of kraken"
date