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
contigs=$1
outfile=$2

echo "Contigs file: $contigs"
echo "Output BAM file: $outfile"

echo "Started bowtie2 mapping contigs to viral genome"
date

## Run Bowties2
apptainer exec "$CONTAINER" bowtie2 --help
apptainer exec "$CONTAINER" bowtie2 \
  --time \
  --threads 4 \
  --local \
  -f \
  -x menuka_metatrans/results/bowtie2_index_rota/rota \
  -U "$contigs" \
| apptainer exec "$samtools" samtools sort -@ 4 -o "$outfile"

apptainer exec "$samtools" samtools index "$outfile"

echo "Successfully ran bowtie2 for mapping contigs to viral genome"
date