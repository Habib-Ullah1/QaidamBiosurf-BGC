#!/bin/bash
#SBATCH --job-name=deseq2
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/data/habib/metagenome/biosurfactant/logs/deseq2_%j.log

source ~/.bashrc
conda activate r_deseq2

echo "=== DESeq2 seasonal analysis $(date) ==="
Rscript /data/habib/metagenome/biosurfactant/scripts/deseq2_seasonal.R
echo "=== Job done $(date) ==="
