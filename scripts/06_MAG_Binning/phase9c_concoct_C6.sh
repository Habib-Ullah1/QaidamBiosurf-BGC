#!/bin/bash
#SBATCH --job-name=phase9c_C6
#SBATCH --time=99:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --output=/data/habib/metagenome/biosurfactant/logs/phase9c_C6_%j.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate concoct_env

s=C6
CONTIGS=/data/habib/metagenome/biosurfactant/C6/00_contigs_filtered/contigs_500bp.fa
BAM=/data/habib/metagenome/biosurfactant/C6/09_binning/bam/C6_sorted.bam
BINDIR=/data/habib/metagenome/biosurfactant/C6/09_binning
LOGDIR=/data/habib/metagenome/biosurfactant/logs

echo "=== Phase 9c CONCOCT: C6 ==="
echo "Started: $(date)"

mkdir -p $BINDIR/concoct/cut $BINDIR/concoct/output

# Step 1 - Cut contigs into chunks (CONCOCT works on 10kb chunks)
echo "Cutting contigs into 10kb chunks..."
cut_up_fasta.py $CONTIGS \
  -c 10000 -o 0 --merge_last \
  -b $BINDIR/concoct/cut/${s}_contigs_10k.bed \
  > $BINDIR/concoct/cut/${s}_contigs_10k.fa

# Step 2 - Generate coverage table
echo "Generating CONCOCT coverage table..."
concoct_coverage_table.py \
  $BINDIR/concoct/cut/${s}_contigs_10k.bed \
  $BAM \
  > $BINDIR/concoct/cut/C6_coverage_table.tsv

# Step 3 - Run CONCOCT
echo "Running CONCOCT clustering..."
concoct \
  --composition_file $BINDIR/concoct/cut/${s}_contigs_10k.fa \
  --coverage_file $BINDIR/concoct/cut/${s}_coverage_table.tsv \
  -b $BINDIR/concoct/output/ \
  --threads 16 \
  2> $LOGDIR/concoct_C6.log

# Step 4 - Merge clustering results back to original contigs
echo "Merging sub-contig clustering..."
merge_cutup_clustering.py \
  $BINDIR/concoct/output/clustering_gt1000.csv \
  > $BINDIR/concoct/output/clustering_merged.csv

# Step 5 - Extract bins as FASTA files
echo "Extracting bins..."
mkdir -p $BINDIR/concoct/bins
extract_fasta_bins.py $CONTIGS \
  $BINDIR/concoct/output/clustering_merged.csv \
  --output_path $BINDIR/concoct/bins/

echo "CONCOCT bins: $(ls $BINDIR/concoct/bins/*.fa 2>/dev/null | wc -l)"
echo "=== Phase 9c C6 COMPLETE ==="
echo "Finished: $(date)"
