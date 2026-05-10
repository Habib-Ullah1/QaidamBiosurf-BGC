#!/bin/bash
#SBATCH --job-name=eggnog_H6
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=64G
#SBATCH --output=/data/habib/metagenome/biosurfactant/logs/eggnog_H6_%j.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate eggnog_env

PROJECT=/data/habib/metagenome
EGGNOG_DB=$PROJECT/biosurfactant/databases/eggnog

echo "=== eggNOG: H6 $(date) ==="
emapper.py \
    -i $PROJECT/biosurfactant/H6/01_prodigal/H6_proteins.faa \
    --output H6_eggnog \
    --output_dir $PROJECT/biosurfactant/H6/05_eggnog \
    --data_dir $EGGNOG_DB \
    --itype proteins \
    --cpu 16 \
    --dmnd_db $EGGNOG_DB/eggnog_proteins.dmnd \
    --override \
    2> $PROJECT/biosurfactant/logs/eggnog_H6_run.log

echo "Annotations: $(grep -v '^#' $PROJECT/biosurfactant/H6/05_eggnog/H6_eggnog.emapper.annotations 2>/dev/null | wc -l)"
echo "=== Done: H6 $(date) ==="
