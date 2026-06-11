#!/bin/bash
#SBATCH --account=PAS3107
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=16GB
#SBATCH --output=slurm-samtool_convert_fastq-%j.out

## Bash strict settings
set -euo pipefail

# Define constant
container=oras://community.wave.seqera.io/library/samtools:1.21--84c9d77c3901e90b

# Define variable
bam_infile=$1
outfile1=$2
outfile2=$3
outfile3=$4

## Run samtools to convert bam file to fastqc file
# samtool documentation: https://www.htslib.org/doc/samtools.html
# flag option https://www.htslib.org/doc/samtools-flags.html
apptainer exec $container samtools fastq "$bam_infile" \
  -1 "$outfile1" \
  -2 "$outfile2" \
  -0 /dev/null \
  -s "$outfile3"

echo "Done with samtools"
date