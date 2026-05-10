#!/bin/bash
#SBATCH --job-name=checkm2_all
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=80G
#SBATCH --output=/data/habib/metagenome/biosurfactant/logs/checkm2_all_%j.log

source ~/.bashrc
conda activate checkm2_new

PROJECT=/data/habib/metagenome

for s in C3 C6 H2 H6; do
    echo "=== CheckM2: $s $(date) ==="
    checkm2 predict \
        --input  $PROJECT/biosurfactant/$s/09_binning/metawrap_refined/bins \
        --output-directory $PROJECT/biosurfactant/$s/10_checkm2 \
        --extension .fa \
        --threads 16 \
        --force
    echo "=== Done: $s $(date) ==="
done
