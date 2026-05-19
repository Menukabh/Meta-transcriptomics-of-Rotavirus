#!/bin/bash
#SBATCH --account=PAS3107
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=36
#SBATCH --mem=120GB
#SBATCH --output=slurm-kraken_database_create-%j.out

## Define variable
kraken_stan=$1
output_directory=$2
echo "The output directory is $output_directory"

##Download Kraken standard database
wget "$kraken_stan" -P "$output_directory"
#tar -xvf "$output_directory"/k2_standard_20251015.tar.gz -C "$output_directory"
tar -xvf "$output_directory"/k2_gtdb_genome_reps_20250609.tar.gz -C "$output_directory"

echo "Done with script"
date