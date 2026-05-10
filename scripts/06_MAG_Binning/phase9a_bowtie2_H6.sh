#!/bin/bash
#SBATCH --job-name=phase9a_H6
#SBATCH --time=99:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/data/habib/metagenome/biosurfactant/logs/phase9a_H6_%j.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biosurfactant_base

s=H6
CONTIGS=/data/habib/metagenome/biosurfactant/H6/00_contigs_filtered/contigs_500bp.fa
R1=/data/habib/metagenome/analysis/MEDUSA/sample_H6/Pipeline/data/trimmed/H6.R1_trim.fastq
R2=/data/habib/metagenome/analysis/MEDUSA/sample_H6/Pipeline/data/trimmed/H6.R2_trim.fastq
OUTDIR=/data/habib/metagenome/biosurfactant/H6/09_binning
LOGDIR=/data/habib/metagenome/biosurfactant/logs

echo "=== Phase 9a Bowtie2: H6 ==="
echo "Started: $(date)"

mkdir -p $OUTDIR/bowtie2_index $OUTDIR/bam

# Step 1 - Build Bowtie2 index
echo "Building index..."
bowtie2-build --threads 16 $CONTIGS $OUTDIR/bowtie2_index/${s}_index \
  2> $LOGDIR/bowtie2_index_H6.log
echo "Index done: $(date)"

# Step 2 - Map reads, sort and index BAM
echo "Mapping reads..."
bowtie2 --threads 16 \
  -x $OUTDIR/bowtie2_index/${s}_index \
  -1 $R1 -2 $R2 \
  --no-unal -q \
  2> $LOGDIR/bowtie2_map_H6.log | \
  samtools sort -@ 16 -o $OUTDIR/bam/${s}_sorted.bam
echo "Mapping done: $(date)"

# Step 3 - Index BAM
samtools index $OUTDIR/bam/${s}_sorted.bam

# Step 4 - Mapping stats
echo "--- Mapping stats ---"
samtools flagstat $OUTDIR/bam/${s}_sorted.bam

echo "=== Phase 9a H6 COMPLETE ==="
echo "Finished: $(date)"
