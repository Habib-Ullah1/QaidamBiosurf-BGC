#!/bin/bash
#SBATCH -J precheck
#SBATCH -p cu
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time 00:20:00
#SBATCH -o /data/habib/metagenome/biosurfactant/logs/%x-%j.out
#SBATCH -e /data/habib/metagenome/biosurfactant/logs/%x-%j.err

source ~/.bashrc
cd /data/habib/metagenome/biosurfactant || exit 1
mkdir -p logs

echo "===== A. is featureCounts available? (check the few likely envs only, not all) ====="
for e in biosurfactant_base subread featurecounts metawrap_env quast_env; do
  p=$(conda run -n "$e" bash -lc 'command -v featureCounts' 2>/dev/null)
  [ -n "$p" ] && echo "  FOUND in $e -> $p"
done
echo "  (if nothing printed above, featureCounts is not installed -> we use samtools bedcov fallback)"
command -v samtools && conda run -n biosurfactant_base bash -lc 'command -v samtools'

echo; echo "===== B. Prodigal headers carry coordinates? (needed to build gene intervals) ====="
grep -m1 '>' C3/01_prodigal/*.faa

echo; echo "===== C. BAM paths + are they indexed? ====="
for S in C3 C6 H2 H6; do
  b=$(ls $S/09_binning/bam/*.bam 2>/dev/null | head -1)
  echo "$S: BAM=$b  index=$(ls ${b}.bai 2>/dev/null || echo MISSING)"
done

echo; echo "===== D. quick BAM sanity (first one) ====="
B=$(ls C3/09_binning/bam/*.bam 2>/dev/null | head -1)
conda run -n biosurfactant_base bash -lc "samtools quickcheck '$B' && echo 'BAM ok' || echo 'BAM problem'"
echo "DONE"
