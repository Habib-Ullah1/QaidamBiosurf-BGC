#!/bin/bash
#SBATCH --job-name=phase9b_C3
#SBATCH --time=99:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --output=/data/habib/metagenome/biosurfactant/logs/phase9b_C3_%j.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate metabat2_env

s=C3
CONTIGS=/data/habib/metagenome/biosurfactant/C3/00_contigs_filtered/contigs_500bp.fa
BAM=/data/habib/metagenome/biosurfactant/C3/09_binning/bam/C3_sorted.bam
BINDIR=/data/habib/metagenome/biosurfactant/C3/09_binning
LOGDIR=/data/habib/metagenome/biosurfactant/logs

echo "=== Phase 9b Binning: C3 ==="
echo "Started: $(date)"

mkdir -p $BINDIR/metabat2 $BINDIR/maxbin2 $BINDIR/depth

# Step 1 - Generate depth file
echo "Generating depth file..."
jgi_summarize_bam_contig_depths \
  --outputDepth $BINDIR/depth/${s}_depth.txt \
  $BAM
echo "Depth done: $(date)"

# Step 2 - MetaBAT2
echo "Running MetaBAT2..."
metabat2 \
  -i $CONTIGS \
  -a $BINDIR/depth/${s}_depth.txt \
  -o $BINDIR/metabat2/${s}_bin \
  -t 16 -m 1500 --unbinned \
  2> $LOGDIR/metabat2_C3.log
echo "MetaBAT2 bins: $(ls $BINDIR/metabat2/*.fa 2>/dev/null | wc -l)"
echo "MetaBAT2 done: $(date)"

# Step 3 - MaxBin2
echo "Running MaxBin2..."
awk '{print $1"\t"$3}' $BINDIR/depth/${s}_depth.txt | \
  tail -n +2 > $BINDIR/depth/${s}_abundance.txt

run_MaxBin.pl \
  -contig $CONTIGS \
  -abund $BINDIR/depth/${s}_abundance.txt \
  -out $BINDIR/maxbin2/${s}_bin \
  -thread 16 \
  2> $LOGDIR/maxbin2_C3.log
echo "MaxBin2 bins: $(ls $BINDIR/maxbin2/*.fasta 2>/dev/null | wc -l)"
echo "MaxBin2 done: $(date)"

echo "=== Phase 9b C3 COMPLETE ==="
echo "Finished: $(date)"
