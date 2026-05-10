#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate biosurfactant_base

PROJECT=/data/habib/metagenome
SAMPLES="C3 C6 H2 H6"

for s in $SAMPLES; do
  echo "=== Prodigal ORF prediction for sample: $s ==="
  IN=$PROJECT/biosurfactant/$s/00_contigs_filtered/contigs_500bp.fa
  OUTDIR=$PROJECT/biosurfactant/$s/01_prodigal

  prodigal \
    -i $IN \
    -a $OUTDIR/${s}_proteins.faa \
    -d $OUTDIR/${s}_genes.fna \
    -f gff \
    -o $OUTDIR/${s}_prodigal.gff \
    -p meta \
    -q \
    2> $PROJECT/biosurfactant/logs/prodigal_${s}.log

  COUNT=$(grep -c '>' $OUTDIR/${s}_proteins.faa)
  echo "  ORFs predicted: $COUNT"
done

echo "Phase 1b COMPLETE"
