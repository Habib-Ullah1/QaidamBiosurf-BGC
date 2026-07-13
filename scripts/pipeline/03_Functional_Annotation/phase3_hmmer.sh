#!/bin/bash
#SBATCH --job-name=phase3_hmmer
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=/data/habib/metagenome/biosurfactant/logs/phase3_hmmer_%j.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biosurfactant_base

PROJECT=/data/habib/metagenome
PFAM=/data/db/antismashDB/pfam/35.0/Pfam-A.hmm
TIGRFAM=/data/db/antismashDB/tigrfam/TIGRFam.hmm
SAMPLES="C3 C6 H2 H6"

for s in $SAMPLES; do
  echo "=== HMMER for sample: $s ==="
  echo "Started: $(date)"

  PROTEINS=$PROJECT/biosurfactant/$s/01_prodigal/${s}_proteins.faa
  OUTDIR=$PROJECT/biosurfactant/$s/03_hmmer

  # Search against Pfam-A
  echo "  Running Pfam search..."
  hmmsearch \
    --cpu 16 \
    --domtblout $OUTDIR/${s}_pfam.domtblout \
    -E 1e-5 \
    --domE 1e-5 \
    $PFAM \
    $PROTEINS \
    > $OUTDIR/${s}_pfam.out \
    2> $PROJECT/biosurfactant/logs/hmmer_pfam_${s}.log

  PFAM_HITS=$(grep -v "^#" $OUTDIR/${s}_pfam.domtblout | wc -l)
  echo "  Pfam hits: $PFAM_HITS"

  # Search against TIGRFAM
  echo "  Running TIGRFAM search..."
  hmmsearch \
    --cpu 16 \
    --domtblout $OUTDIR/${s}_tigrfam.domtblout \
    -E 1e-5 \
    --domE 1e-5 \
    $TIGRFAM \
    $PROTEINS \
    > $OUTDIR/${s}_tigrfam.out \
    2> $PROJECT/biosurfactant/logs/hmmer_tigrfam_${s}.log

  TIGR_HITS=$(grep -v "^#" $OUTDIR/${s}_tigrfam.domtblout | wc -l)
  echo "  TIGRFAM hits: $TIGR_HITS"

  echo "Finished: $(date)"
  echo ""
done

echo "Phase 3 COMPLETE"
EOF

head -3 $PROJECT/biosurfactant/scripts/phase3_hmmer.sh
