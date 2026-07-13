#!/bin/bash
#SBATCH -J locate_tools
#SBATCH -p cu
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time 00:20:00
#SBATCH -o /data/habib/metagenome/biosurfactant/logs/%x-%j.out
#SBATCH -e /data/habib/metagenome/biosurfactant/logs/%x-%j.err

export LC_ALL=C.UTF-8; export LANG=C.UTF-8
cd /data/habib/metagenome/biosurfactant || exit 1

echo "===== 1. environment modules? ====="
module avail 2>&1 | grep -iE "samtools|subread|bedtools" || echo "  (no module system / none found)"

echo; echo "===== 2. direct filesystem search in conda envs (no activation needed) ====="
for exe in samtools featureCounts bedtools; do
  echo "--- $exe ---"
  find /home/habib/miniconda3 -maxdepth 4 -type f -name "$exe" 2>/dev/null
done

echo; echo "===== 3. anywhere else in home / data (capped depth) ====="
find /home/habib -maxdepth 5 -type f \( -name samtools -o -name featureCounts \) 2>/dev/null | grep -v miniconda | head

echo; echo "===== 4. what env did the binning scripts actually use? (they ran samtools) ====="
grep -rhiE "conda activate|module load|samtools" scripts/06_MAG_Binning/*.sh 2>/dev/null | grep -iE "activate|module|samtools" | sort -u | head

echo DONE
