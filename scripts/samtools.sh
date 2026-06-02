#!/bin/bash
#SBATCH --account=PAS3107
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=16GB
#SBATCH --output=slurm-samtool_output-%j.out

## Bash strict settings
set -euo pipefail

# Define constant
container=oras://community.wave.seqera.io/library/samtools:1.21--84c9d77c3901e90b

# Define variable
infile=$1
outfile=$2

## Run samtools to filter out the unampped reads
# samtool documentation: https://www.htslib.org/doc/samtools.html
#apptainer exec $container samtools view --help
# flag option https://www.htslib.org/doc/samtools-flags.html
apptainer exec $container samtools view  -@ 8 -h -F 4 -b "$infile" > "$outfile"

echo "Done with samtools"
date