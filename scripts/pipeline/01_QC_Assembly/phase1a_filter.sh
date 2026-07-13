#!/bin/bash
#SBATCH --job-name=phase1a_filter
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/data/habib/metagenome/biosurfactant/logs/phase1a_filter_%j.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biosurfactant_base

PROJECT=/data/habib/metagenome
SAMPLES="C3 C6 H2 H6"

for s in $SAMPLES; do
  echo "=== Filtering contigs for sample: $s ==="
  IN=/data/habib/metagenome/analysis/MEDUSA/sample_${s}/Pipeline/data/assembled/${s}_assembly/final.contigs.fa
  OUT=$PROJECT/biosurfactant/$s/00_contigs_filtered/contigs_500bp.fa
  LOG=$PROJECT/biosurfactant/logs/seqkit_filter_$s.log

  BEFORE=$(grep -c '>' $IN)

  seqkit seq \
    --min-len 500 \
    --out-file $OUT \
    $IN 2> $LOG

  AFTER=$(grep -c '>' $OUT)
  echo "  Before: $BEFORE contigs"
  echo "  After : $AFTER contigs retained"
  echo "  Log   : $LOG"
done

echo "Phase 1a COMPLETE"
