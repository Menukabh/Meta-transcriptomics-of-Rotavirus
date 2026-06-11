#!/bin/bash
#SBATCH --account=PAS3107
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=150GB
#SBATCH --output=slurm-kraken_output-%j.out

##Define constant
kraken_cont=oras://community.wave.seqera.io/library/kraken2:2.17.1--1738c34504f3fb18

## Create directory
#mkdir -p menuka_metatrans/results/kraken_bacteria

## Define variables
#kraken_db="/fs/scratch/PAS0471/menuka/kraken_stnd_db"
kraken_db=/fs/scratch/PAS0471/menuka/kraken_host_db_metatrans
echo "Using kraken2 standard database: $kraken_db"
outfile=$1
report=$2
Read1=$3
Read2=$4

echo "The output file is $outfile"
echo "The report file is $report"
echo "The forward read is $Read1"
echo "The reverse read is $Read2"

## Run kraken2 to separate viral reads
apptainer exec $kraken_cont kraken2 --db "$kraken_db" --threads 12 \
        --output "$outfile" --report "$report" \
        --paired "$Read1" "$Read2"

echo "Done with script"
date