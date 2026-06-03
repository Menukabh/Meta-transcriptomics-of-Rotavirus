#!/bin/bash
#SBATCH --account=PAS3107
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=150GB
#SBATCH --output=slurm-krakendb-%j.out

# Define container and db
kraken_cont=oras://community.wave.seqera.io/library/kraken2:2.17.1--1738c34504f3fb18
kraken_db="/fs/scratch/PAS0471/menuka/kraken_stnd_db"

# build the db after adding pig and cell genome sequence
apptainer exec $kraken_cont kraken2-build \
  --build \
  --threads 24 \
  --db $kraken_db

echo "Added sequences to the standard db"
date