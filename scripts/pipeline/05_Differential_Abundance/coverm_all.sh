#!/bin/bash
#SBATCH --job-name=coverm_all
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=60G
#SBATCH --output=/data/habib/metagenome/biosurfactant/logs/coverm_%j.log

source ~/.bashrc
conda activate biosurfactant_base

PROJECT=/data/habib/metagenome
INTEG=$PROJECT/biosurfactant/integration
mkdir -p $INTEG

echo "=== CoverM: collecting BAM files $(date) ==="

# Collect all BAM files from the binning directories
BAMFILES=""
for s in C3 C6 H2 H6; do
    for bam in $PROJECT/biosurfactant/$s/09_binning/bam/*.bam; do
        [ -f "$bam" ] && BAMFILES="$BAMFILES $bam"
    done
done

echo "BAM files found: $(echo $BAMFILES | wc -w)"

# Run CoverM on contigs for all samples
for s in C3 C6 H2 H6; do
    echo "=== CoverM contig coverage: $s $(date) ==="
    mkdir -p $PROJECT/biosurfactant/$s/13_coverage/coverM

    coverm contig \
        --bam-files $PROJECT/biosurfactant/$s/09_binning/bam/*.bam \
        --methods mean trimmed_mean covered_fraction \
        --min-read-percent-identity 95 \
        --min-read-aligned-percent 75 \
        --threads 16 \
        --output-file $PROJECT/biosurfactant/$s/13_coverage/coverM/${s}_contig_coverage.tsv

    echo "=== Done: $s $(date) ==="
done

echo "=== All CoverM done $(date) ==="
