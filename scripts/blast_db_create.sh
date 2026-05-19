#!/bin/bash
#SBATCH --account=PAS3107
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --output=slurm-blastdb_create_%j.output

# bash strict settings
set -euo pipefail

# Load blast module
module load blast-plus/2.16.0

# Define variables
QUERY=$1

echo "Start to create blast db"
date

# Create blast database of Rotavirus
# makeblastdb -help
makeblastdb -in "$QUERY" \
 -dbtype nucl \
 -parse_seqids \
 -out menuka_metatrans/results/rota_db/NCBI_rota_genome

echo "Blast database created at menuka_metatrans/results/rota_db"
date