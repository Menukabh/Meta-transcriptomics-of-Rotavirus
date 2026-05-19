#!/bin/bash
#SBATCH --account=PAS3107
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=40GB
#SBATCH --output=slurm-trimmomatic_output-%j.out

set -euo pipefail

## Define constant
CONTAINER=oras://community.wave.seqera.io/library/trimmomatic:0.40--7b5b7590373e6fc4

## Define variable
R1=$1
R2=$2
outfile1_pair=$3
outfile1_unpair=$4
outfile2_pair=$5
outfile2_unpair=$6

# Run trimmomatic - remove Nextera adapters
# sliding window of 4 removes bases quality below 20 in each 4 bases

# apptainer exec "$CONTAINER" trimmomatic PE --help

#apptainer exec "$CONTAINER" trimmomatic PE \
    #-threads 24 \
    #$R1 $R2 \
    #$outfile1_pair $outfile1_unpair \
    #$outfile2_pair $outfile2_unpair \
    #ILLUMINACLIP:/fs/ess/PAS3107/menuka_metatrans/scripts/NexteraPE-PE.fa:2:30:10 \
    #LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

# For multiple samples: Got an error java.lang.OutOfMemoryError: Java heap space
# use larger java heap and reduce the thread
apptainer exec --env JAVA_TOOL_OPTIONS="-Xmx32g" "$CONTAINER" trimmomatic PE \
    -threads 8 \
    "$R1" "$R2" \
    "$outfile1_pair" "$outfile1_unpair" \
    "$outfile2_pair" "$outfile2_unpair" \
    ILLUMINACLIP:/fs/ess/PAS3107/menuka_metatrans/scripts/NexteraPE-PE.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

echo 'Done with the trimmomatic script'
date