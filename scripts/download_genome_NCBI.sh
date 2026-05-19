#!/bin/bash
#SBATCH --account=PAS3107
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=12G
#SBATCH --output=slurm-rotavirus_download_%j.output

## Bash strict settings
set -euo pipefail

# Load conda env
module load miniconda3/24.1.2-py310
conda activate /fs/ess/PAS0471/condaenv_database_Menuka/ncbi_env/

echo "Start to download rotavirus genome from NCBI"
date

## Download the dataset from NCBI of the rotavirus genome
datasets download virus genome taxon 10912 \
--include genome \
--filename menuka_metatrans/data/NCBI_rota_genome/rotavirus_data.zip

echo "Done with script"
date