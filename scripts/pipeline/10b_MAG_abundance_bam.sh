#!/bin/bash
#SBATCH --job-name=coverm_bam
#SBATCH --partition=cu
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --output=/data/habib/metagenome/biosurfactant/logs/coverm_bam_%j.log
set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda activate coverm_env
cd /data/habib/metagenome/biosurfactant

R=/data/habib/metagenome/ncbi_upload
OUT=16_mag_abundance
REF=$OUT/all_MAGs_concat.fna
BAMDIR=$OUT/bams
export TMPDIR=$OUT/tmp2
mkdir -p "$BAMDIR" "$TMPDIR"

[ -s "$REF" ]     || cat $OUT/all_MAGs_renamed/*.fa > "$REF"
[ -s "$REF.mmi" ] || minimap2 -x sr -d "$REF.mmi" "$REF"

map_one () {
  local name=$1 r1=$2 r2=$3
  if [ -s "$BAMDIR/$name.bam" ]; then echo "[$name] exists, skip"; return; fi
  echo "[$name] mapping $(date)"
  minimap2 -ax sr -t 16 "$REF.mmi" "$r1" "$r2" \
    | samtools sort -@ 8 -m 2G -T "$TMPDIR/$name.st" -o "$BAMDIR/$name.bam" -
  samtools index -@ 4 "$BAMDIR/$name.bam"
  echo "[$name] done $(date)"
}
map_one C3 $R/C3_1.fq.gz  $R/C3_2.fq.gz
map_one C6 $R/C6.R1.fq.gz $R/C6.R2.fq.gz
map_one H2 $R/H2.R1.fq.gz $R/H2.R2.fq.gz
map_one H6 $R/H6.R1.fq.gz $R/H6.R2.fq.gz

coverm genome --bam-files $BAMDIR/*.bam \
  --genome-fasta-directory $OUT/all_MAGs_renamed -x fa \
  -m relative_abundance mean covered_fraction --min-covered-fraction 10 \
  -t 16 -o $OUT/MAG_relabund.tsv
echo "EXIT=$?"; wc -l $OUT/MAG_relabund.tsv
