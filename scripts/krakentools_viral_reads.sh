#!/bin/bash
#SBATCH --account=PAS3107
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=36GB
#SBATCH --output=slurm-kraken_-%j.out

set -euo pipefail

# load the conda env
module load miniconda3/24.1.2-py310
conda activate /users/PAS0471/menuka/.conda/envs/krakentools

# Define variable
outfile=$1
report=$2
Read1=$3
Read2=$4
viral_reads1=$5
viral_reads2=$6

# Run krakentools: --fastq-output: To get the output in fastq format, the default output is FASTA
#--include-children: includes reads classified to any descendant of 10239/viral reads
#KrakenTools/extract_kraken_reads.py --help
KrakenTools/extract_kraken_reads.py \
  -k "$outfile" \
  -r "$report" \
  -s1 "$Read1" \
  -s2 "$Read2" \
  -t 10239 \
  --include-children \
  --fastq-output \
  -o "$viral_reads1" \
  -o2 "$viral_reads2"

echo 'Done with the script krakentool'
date
