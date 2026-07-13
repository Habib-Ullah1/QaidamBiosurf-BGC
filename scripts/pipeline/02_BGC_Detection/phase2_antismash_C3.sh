#!/bin/bash
#SBATCH --job-name=antismash_C3
#SBATCH --time=99:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --output=/data/habib/metagenome/biosurfactant/logs/phase2_antismash_C3_%j.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate antismash

PROJECT=/data/habib/metagenome
s=C3

echo "=== antiSMASH BGC detection for sample: $s ==="
echo "Started: $(date)"

IN=$PROJECT/biosurfactant/$s/00_contigs_filtered/contigs_500bp.fa
OUTDIR=$PROJECT/biosurfactant/$s/02_antismash

antismash \
  --cpus 8 \
  --minlength 1000 \
  --genefinding-tool prodigal-m \
  --asf \
  --rre \
  --tigrfam \
  --clusterhmmer \
  --output-dir $OUTDIR \
  --output-basename $s \
  --logfile $PROJECT/biosurfactant/logs/antismash_${s}.log \
  $IN

echo "Finished: $(date)"
echo "Phase 2 $s COMPLETE"
