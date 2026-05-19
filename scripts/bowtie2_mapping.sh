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
samtools=oras://community.wave.seqera.io/library/samtools:1.23.1--5cb989b890127f7a

## Define variables
R1=$1
R2=$2
outfile=$3

echo "R1 file: $R1"
echo "R2 file: $R2"
echo "Output file: $outfile"

echo "Started script for bowtie2 for mapping raw reads to viral genome"
date

## Run Bowties2
apptainer exec "$CONTAINER" bowtie2 \
--time \
--threads 4 \
-x menuka_metatrans/results/bowtie2_index_rota/rota \
-1 "$R1" \
-2 "$R2" \
| apptainer exec "$samtools" samtools sort -@ 4 -o "$outfile"

echo "Successfully Ran bowtie2 for mapping raw reads to viral genome"
date