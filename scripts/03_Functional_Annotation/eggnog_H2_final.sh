#!/bin/bash
#SBATCH --job-name=eggnog_H2
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=64G
#SBATCH --time=99:00:00
#SBATCH --output=/data/habib/metagenome/biosurfactant/logs/eggnog_H2_final_%j.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate eggnog_env

PROJECT=/data/habib/metagenome
EGGNOG_DB=/data/habib/metagenome/biosurfactant/databases/eggnog

echo "=== eggNOG FINAL RESUME: H2 $(date) ==="
echo "Current annotations: $(grep -v '^#' /data/habib/metagenome/biosurfactant/H2/05_eggnog/H2_eggnog.emapper.annotations 2>/dev/null | wc -l)"

emapper.py \
    -i /data/habib/metagenome/biosurfactant/H2/01_prodigal/H2_proteins.faa \
    --output H2_eggnog \
    --output_dir /data/habib/metagenome/biosurfactant/H2/05_eggnog \
    --data_dir $EGGNOG_DB \
    --itype proteins \
    --cpu 16 \
    --dmnd_db $EGGNOG_DB/eggnog_proteins.dmnd \
    --resume \
    2> /data/habib/metagenome/biosurfactant/logs/eggnog_H2_final.log

echo "Final annotations: $(grep -v '^#' /data/habib/metagenome/biosurfactant/H2/05_eggnog/H2_eggnog.emapper.annotations 2>/dev/null | wc -l)"
echo "=== Done: H2 $(date) ==="
