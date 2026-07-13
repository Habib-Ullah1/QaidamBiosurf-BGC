#!/bin/bash
#SBATCH -J find_tools
#SBATCH -p cu
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time 00:20:00
#SBATCH -o /data/habib/metagenome/biosurfactant/logs/%x-%j.out
#SBATCH -e /data/habib/metagenome/biosurfactant/logs/%x-%j.err

source ~/.bashrc
cd /data/habib/metagenome/biosurfactant || exit 1

echo "===== which envs have samtools / featureCounts / bedtools? ====="
for e in $(conda env list | awk 'NR>2 && $1!~"#"{print $1}'); do
  out=$(conda run -n "$e" bash -lc 'for t in samtools featureCounts bedtools; do p=$(command -v $t 2>/dev/null); [ -n "$p" ] && echo "  $t"; done' 2>/dev/null)
  [ -n "$out" ] && { echo "ENV $e:"; echo "$out"; }
done
echo "===== also check the metawrap/binning env explicitly (it sorted the BAMs) ====="
for e in metawrap_env metawrap binning_env coverm_env; do
  conda run -n "$e" bash -lc 'command -v samtools featureCounts 2>/dev/null' 2>/dev/null | sed "s/^/  $e: /"
done
echo DONE
