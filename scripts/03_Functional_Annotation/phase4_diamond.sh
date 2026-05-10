#!/bin/bash
#SBATCH --job-name=phase4_diamond
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=/data/habib/metagenome/biosurfactant/logs/phase4_diamond_%j.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biosurfactant_base

PROJECT=/data/habib/metagenome
DBDIR=$PROJECT/biosurfactant/databases/diamond_refs
SAMPLES="C3 C6 H2 H6"

echo "Building DIAMOND database..."
diamond makedb --in $DBDIR/biosurfactant_uniprot.fasta --db $DBDIR/biosurfactant_refs --threads 16

for s in $SAMPLES; do
  mkdir -p $PROJECT/biosurfactant/$s/04_diamond
  PROTEINS=$PROJECT/biosurfactant/$s/01_prodigal/${s}_proteins.faa
  OUTDIR=$PROJECT/biosurfactant/$s/04_diamond
  echo "=== DIAMOND blastp: $s ==="
  diamond blastp --db $DBDIR/biosurfactant_refs --query $PROTEINS --out $OUTDIR/${s}_diamond.tsv --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle --evalue 1e-5 --id 30 --query-cover 50 --threads 16 --sensitive 2> $PROJECT/biosurfactant/logs/diamond_${s}.log
  echo "  Hits: $(wc -l < $OUTDIR/${s}_diamond.tsv)"
done

echo "Phase 4 COMPLETE"
