#!/bin/bash
#SBATCH --job-name=antismash_1500
#SBATCH --array=1-4
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/data/habib/metagenome/biosurfactant/logs/antismash_1500_%a.out
#SBATCH --error=/data/habib/metagenome/biosurfactant/logs/antismash_1500_%a.err

samples=(C3 C6 H2 H6)
s=${samples[$SLURM_ARRAY_TASK_ID-1]}

source ~/miniconda3/etc/profile.d/conda.sh
conda activate antismash

INPUT=/data/habib/metagenome/biosurfactant/$s/00_contigs_filtered/contigs_1500bp.fa
OUTDIR=/data/habib/metagenome/biosurfactant/$s/02_antismash_1500bp

rm -rf $OUTDIR
mkdir -p $OUTDIR

antismash \
  --cpus 16 \
  --taxon bacteria \
  --genefinding-tool prodigal-m \
  --output-dir $OUTDIR \
  $INPUT

echo "$s antismash 1500bp done"
