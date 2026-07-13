#!/bin/bash
#SBATCH --job-name=phase11_C6
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --output=/data/habib/metagenome/biosurfactant/logs/phase11_C6_%j.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate gtdbtk_env

export GTDBTK_DATA_PATH=/data/db/gtdbtk_R214/release214

s=C6
BINS=/data/habib/metagenome/biosurfactant/C6/09_binning/metawrap_refined/metawrap_50_10_bins
OUTDIR=/data/habib/metagenome/biosurfactant/C6/11_gtdbtk

echo "=== Phase 11 GTDB-Tk: C6 ==="
echo "Started: $(date)"
mkdir -p $OUTDIR

gtdbtk classify_wf \
  --genome_dir $BINS \
  --out_dir $OUTDIR \
  --extension fa \
  --cpus 16 \
  --skip_ani_screen \
  2> /data/habib/metagenome/biosurfactant/logs/gtdbtk_C6.log

echo "--- GTDB-Tk Summary ---"
cat $OUTDIR/gtdbtk.bac120.summary.tsv 2>/dev/null | head -20
cat $OUTDIR/gtdbtk.ar53.summary.tsv 2>/dev/null | head -10

echo "=== Phase 11 C6 COMPLETE ==="
echo "Finished: $(date)"
