#!/bin/bash
#SBATCH --job-name=bigscape_all4
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --output=/data/habib/metagenome/biosurfactant/logs/bigscape_all4_%j.log

source ~/.bashrc
conda activate bigscape_env

PROJECT=/data/habib/metagenome
PFAM=$PROJECT/biosurfactant/databases/pfam
BGCDIR=$PROJECT/biosurfactant/bigscape_input_3samples

echo "=== BiG-SCAPE all 4 samples $(date) ==="
echo "Input BGCs: $(ls $BGCDIR/*.gbk | wc -l)"

bigscape \
    --inputdir $BGCDIR \
    --outputdir $PROJECT/biosurfactant/bigscape_output_all4 \
    --pfam_dir $PFAM \
    --cores 16 \
    --cutoffs 0.3 0.5 0.7 \
    --mix \
    --include_gbk_str region \
    --verbose

echo "=== BiG-SCAPE done $(date) ==="
ls $PROJECT/biosurfactant/bigscape_output_all4/
