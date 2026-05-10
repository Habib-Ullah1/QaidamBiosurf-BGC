#!/bin/bash
#SBATCH --job-name=bigscape_3s
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --output=/data/habib/metagenome/biosurfactant/logs/bigscape_3s_%j.log

source ~/.bashrc
conda activate bigscape_env

PROJECT=/data/habib/metagenome
PFAM=$PROJECT/biosurfactant/databases/pfam
BGCDIR=$PROJECT/biosurfactant/bigscape_input_3samples

echo "=== Verifying Pfam files ==="
ls -lh $PFAM/*.h3* 2>/dev/null || echo "ERROR: Pfam .h3 files missing at $PFAM"

echo "=== Verifying BGC input files ==="
ls $BGCDIR/*.gbk 2>/dev/null | wc -l

echo "=== BiG-SCAPE starting $(date) ==="
bigscape \
    --inputdir $BGCDIR \
    --outputdir $PROJECT/biosurfactant/bigscape_output_3samples \
    --pfam_dir $PFAM \
    --cores 16 \
    --cutoffs 0.3 0.5 0.7 \
    --mix \
    --include_gbk_str region \
    --verbose

echo "=== BiG-SCAPE done $(date) ==="
ls $PROJECT/biosurfactant/bigscape_output_3samples/
