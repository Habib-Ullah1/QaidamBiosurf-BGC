#!/bin/bash
#SBATCH --job-name=phase10_C6
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/data/habib/metagenome/biosurfactant/logs/phase10_C6_%j.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate checkm2_env

s=C6
BINS=/data/habib/metagenome/biosurfactant/C6/09_binning/metawrap_refined/metawrap_50_10_bins
OUTDIR=/data/habib/metagenome/biosurfactant/C6/10_checkm2
DB=/data/habib/metagenome/biosurfactant/databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd

echo "=== Phase 10 CheckM2: C6 ==="
echo "Started: $(date)"
mkdir -p $OUTDIR

checkm2 predict \
  --input $BINS \
  --output-directory $OUTDIR \
  --database_path $DB \
  --extension fa \
  --threads 16 \
  2> /data/habib/metagenome/biosurfactant/logs/checkm2_C6.log

echo "--- CheckM2 Summary ---"
cat $OUTDIR/quality_report.tsv | head -20
echo "High quality MAGs (>90% complete, <5% contam):"
awk -F'\t' 'NR>1 && $2>90 && $3<5' $OUTDIR/quality_report.tsv | wc -l
echo "Medium quality MAGs (>50% complete, <10% contam):"
awk -F'\t' 'NR>1 && $2>50 && $3<10' $OUTDIR/quality_report.tsv | wc -l

echo "=== Phase 10 C6 COMPLETE ==="
echo "Finished: $(date)"
