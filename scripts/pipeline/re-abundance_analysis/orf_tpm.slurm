#!/bin/bash
#SBATCH -J orf_tpm
#SBATCH -p cu
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --time 06:00:00
#SBATCH -o /data/habib/metagenome/biosurfactant/logs/%x-%j.out
#SBATCH -e /data/habib/metagenome/biosurfactant/logs/%x-%j.err

export LC_ALL=C.UTF-8; export LANG=C.UTF-8
source ~/.bashrc
conda activate counts
cd /data/habib/metagenome/biosurfactant || exit 1
mkdir -p logs integration/functional/tpm
OUT=integration/functional/tpm
THREADS=28

for S in C3 C6 H2 H6; do
  FAA=$(ls $S/01_prodigal/*.faa | head -1)
  BAM=$S/09_binning/bam/${S}_sorted.bam

  # 1) SAF from Prodigal headers:  ">k141_123_4 # start # end # strand # ..."
  awk '/^>/{
        id=substr($1,2);
        split($0,a," # "); s=a[2]; e=a[3]; st=(a[4]==1?"+":"-");
        contig=id; sub(/_[0-9]+$/,"",contig);
        print id"\t"contig"\t"s"\t"e"\t"st
      }' "$FAA" > $OUT/${S}.saf
  echo "$S: $(wc -l < $OUT/${S}.saf) genes in SAF"

  # 2) count read PAIRS per ORF against the full-assembly BAM
  featureCounts -F SAF -a $OUT/${S}.saf -o $OUT/${S}.counts.tsv \
     -T $THREADS -p --countReadPairs -M --primary "$BAM" 2> $OUT/${S}.fc.log

  # 3) counts -> TPM within the sample (col6=Length, last col=count)
  awk 'BEGIN{FS=OFS="\t"}
       /^#/||/^Geneid/{next}
       {len=$6; c=$NF; if(len>0){rpk=c/(len/1000); R[$1]=rpk; s+=rpk}}
       END{for(g in R) print g, (s>0?R[g]/s*1e6:0)}' $OUT/${S}.counts.tsv > $OUT/${S}.tpm.tsv
  echo "$S TPM: $(wc -l < $OUT/${S}.tpm.tsv) genes"
done
echo "ALL DONE"
