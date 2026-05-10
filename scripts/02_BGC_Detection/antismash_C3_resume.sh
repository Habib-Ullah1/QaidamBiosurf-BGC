#!/bin/bash
#SBATCH --job-name=antismash_C3
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=120G
#SBATCH --output=/data/habib/metagenome/biosurfactant/logs/antismash_C3_resume_%j.log

source ~/.bashrc
conda activate antismash

PROJECT=/data/habib/metagenome
IN=$PROJECT/biosurfactant/C3/00_contigs_filtered/contigs_500bp.fa
OUT=$PROJECT/biosurfactant/C3/02_antismash
LOG=$PROJECT/biosurfactant/logs/antismash_C3.log

echo "=== antiSMASH C3 RESUME $(date) ==="

antismash \
    --taxon bacteria \
    --cpus 16 \
    --minlength 1000 \
    --genefinding-tool prodigal-m \
    --cb-general \
    --cb-subclusters \
    --cb-knownclusters \
    --clusterhmmer \
    --tigrfam \
    --asf \
    --pfam2go \
    --smcog-trees \
    --rre \
    --output-dir $OUT \
    $IN \
    >> $LOG 2>&1

echo "=== antiSMASH C3 done $(date) ==="
REGIONS=$(ls $OUT/*.region*.gbk 2>/dev/null | wc -l)
echo "BGC regions produced: $REGIONS"
