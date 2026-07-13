#!/bin/bash
#SBATCH --job-name=coverm_finish
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/data/habib/metagenome/biosurfactant/16_mag_abundance/coverm_finish.log

set -euo pipefail
B=/data/habib/metagenome/biosurfactant
READS=/data/habib/metagenome/ncbi_upload
MAGDIR=$B/16_mag_abundance
BAMS=$MAGDIR/bams
REF=$MAGDIR/all_MAGs_concat.fna
T=${SLURM_CPUS_PER_TASK:-16}

source ~/miniconda3/etc/profile.d/conda.sh
conda activate coverm_env

echo "$(date) — Step 1: Map H6 to concat reference"
if [ ! -f "$BAMS/H6.bam" ]; then
  minimap2 -ax sr -t "$T" "$REF.mmi" \
    "$READS/H6.R1.fq.gz" "$READS/H6.R2.fq.gz" \
    | samtools sort -@ 4 -m 2G -o "$BAMS/H6.bam" -
  samtools index -@ "$T" "$BAMS/H6.bam"
  echo "$(date) — H6 BAM done: $(ls -lh "$BAMS/H6.bam" | awk '{print $5}')"
else
  echo "$(date) — H6 BAM already exists, skipping"
fi

echo "$(date) — Step 2: Run coverm genome on all 4 BAMs"
coverm genome \
  --bam-files "$BAMS/C3.bam" "$BAMS/C6.bam" "$BAMS/H2.bam" "$BAMS/H6.bam" \
  --genome-fasta-directory "$MAGDIR/all_MAGs_renamed" \
  -x fa \
  -m relative_abundance mean covered_fraction \
  --min-covered-fraction 0 \
  -t "$T" \
  -o "$MAGDIR/MAG_relabund.tsv"

echo "$(date) — Done. Output:"
wc -l "$MAGDIR/MAG_relabund.tsv"
head -5 "$MAGDIR/MAG_relabund.tsv"
echo
echo "=== Top 10 MAGs by relative abundance (sample C3) ==="
head -1 "$MAGDIR/MAG_relabund.tsv"
tail -n+2 "$MAGDIR/MAG_relabund.tsv" | sort -t$'\t' -k2 -rn | head -10
