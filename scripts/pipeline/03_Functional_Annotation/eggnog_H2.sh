#!/bin/bash
#SBATCH --job-name=eggnog_H2
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=64G
#SBATCH --output=/data/habib/metagenome/biosurfactant/logs/eggnog_H2_%j.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate eggnog_env

PROJECT=/data/habib/metagenome
EGGNOG_DB=$PROJECT/biosurfactant/databases/eggnog

echo "=== eggNOG: H2 $(date) ==="
emapper.py \
    -i $PROJECT/biosurfactant/H2/01_prodigal/H2_proteins.faa \
    --output H2_eggnog \
    --output_dir $PROJECT/biosurfactant/H2/05_eggnog \
    --data_dir $EGGNOG_DB \
    --itype proteins \
    --cpu 16 \
    --dmnd_db $EGGNOG_DB/eggnog_proteins.dmnd \
    --override \
    2> $PROJECT/biosurfactant/logs/eggnog_H2_run.log

echo "Annotations: $(grep -v '^#' $PROJECT/biosurfactant/H2/05_eggnog/H2_eggnog.emapper.annotations 2>/dev/null | wc -l)"
echo "=== Done: H2 $(date) ==="
