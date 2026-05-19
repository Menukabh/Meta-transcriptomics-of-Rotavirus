#!/bin/bash
#SBATCH --account=PAS3107
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=12G
#SBATCH --output=slurm-blast_%j.output

## Bash strict settings
set -euo pipefail

# Load container
cont=oras://community.wave.seqera.io/library/blast:2.17.0--3d1eb1104ccfd59c

# Define variables
input_file=$1
db_name=$2
out_file=$3

## Blast the viral contigs against the NCBI ref seq to identify virus - set percentage identity to 80%
apptainer exec "$cont" blastn -query "$input_file" \
    -db "$db_name" \
    -evalue 1e-10 \
    -perc_identity 80 \
    -outfmt "6 qseqid saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
    -num_threads 12 \
    -out "$out_file"
   
echo "Done with script"
date