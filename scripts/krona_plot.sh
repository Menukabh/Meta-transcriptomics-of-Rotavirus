#!/bin/bash
#SBATCH --account=PAS3107
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=16GB
#SBATCH --output=slurm-krona_plot_output-%j.out

## Bash strict settings
set -euo pipefail

# Load conda env
module load miniconda3/24.1.2-py310
conda activate /fs/ess/PAS0471/condaenv_database_Menuka/krona

# Define variable
input=$1
output=$2

## Run kronaplot to create plots of the kraken output
#https://telatin.github.io/microbiome-bioinformatics/Kraken-to-Krona/
ktImportTaxonomy -t 5 -m 3 "$input" -o "$output"

echo "Done with kronaplot"
date