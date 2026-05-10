#!/bin/bash
#SBATCH --job-name=phase9d_H6
#SBATCH --time=99:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --output=/data/habib/metagenome/biosurfactant/logs/phase9d_H6_%j.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate metawrap_env

s=H6
BINDIR=/data/habib/metagenome/biosurfactant/H6/09_binning
LOGDIR=/data/habib/metagenome/biosurfactant/logs

echo "=== Phase 9d MetaWRAP bin refinement: H6 ==="
echo "Started: $(date)"

rm -rf $BINDIR/metawrap_refined
mkdir -p $BINDIR/metawrap_refined

metawrap bin_refinement \
  -o $BINDIR/metawrap_refined \
  -t 16 -m 120 \
  -c 30 -x 15 \
  -A $BINDIR/metabat2 \
  -B $BINDIR/maxbin2_clean \
  -C $BINDIR/concoct/bins \
  2> $LOGDIR/metawrap_H6.log

echo "Total bins: $(ls $BINDIR/metawrap_refined/metawrap_30_15_bins/*.fa 2>/dev/null | wc -l)"
cat $BINDIR/metawrap_refined/metawrap_30_15_bins.stats 2>/dev/null | head -20

echo "=== Phase 9d H6 COMPLETE ==="
echo "Finished: $(date)"
