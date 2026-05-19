#!/bin/bash
#SBATCH --account=PAS3107
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=16GB
#SBATCH --output=slurm-busco_output-%j.out

## Bash strict settings
set -euo pipefail

## Load container and load the either virus_odb10 database
CONTAINER=oras://community.wave.seqera.io/library/busco:6.0.0--7def4b2c35a1aed1
LINEAGE_DIR=/fs/ess/PAS0471/menuka/Erin_data/busco_downloads/lineages/hemiptera_odb10

## Arguments
input=$1
outdir=$2

# Run BUSCO - list the lists of the available datasets
#apptainer exec "$CONTAINER" busco --help
apptainer exec "$CONTAINER" busco --list-datasets
apptainer exec "$CONTAINER" busco \
    -m genome \
    -l  "$LINEAGE_DIR" \
    --offline \
    --download_path /fs/ess/PAS0471/menuka/Erin_data \
    -i "$input" \
    -o "$outdir" \
    -c 30

echo "Done with script"
date