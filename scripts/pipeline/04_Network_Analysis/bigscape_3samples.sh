#!/bin/bash
#SBATCH --job-name=bigscape_3samples
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=80G
#SBATCH --output=/data/habib/metagenome/biosurfactant/logs/bigscape_3samples_%j.log

source ~/.bashrc
conda activate bigscape_env

PROJECT=/data/habib/metagenome

# Collect BGK files from C6, H2, H6 only (C3 antiSMASH still running)
BGCDIR=$PROJECT/biosurfactant/bigscape_input_3samples
mkdir -p $BGCDIR

echo "=== Collecting BGC files from C6, H2, H6 $(date) ==="
for s in C6 H2 H6; do
    ASDIR=$PROJECT/biosurfactant/$s/02_antismash
    N=$(ls $ASDIR/*.region*.gbk 2>/dev/null | wc -l)
    echo "  $s: $N BGC region files"
    for f in $ASDIR/*.region*.gbk; do
        [ -f "$f" ] && cp "$f" "$BGCDIR/${s}_$(basename $f)"
    done
done

TOTAL=$(ls $BGCDIR/*.gbk 2>/dev/null | wc -l)
echo "Total BGC files collected: $TOTAL"

echo "=== BiG-SCAPE starting $(date) ==="
bigscape \
    --inputdir $BGCDIR \
    --outputdir $PROJECT/biosurfactant/bigscape_output_3samples \
    --pfam_dir $PROJECT/biosurfactant/databases/pfam \
    --cores 16 \
    --cutoffs 0.3 0.5 0.7 \
    --mix \
    --include_gbk_str region \
    --verbose

echo "=== BiG-SCAPE done $(date) ==="
echo "GCF families produced:"
ls $PROJECT/biosurfactant/bigscape_output_3samples/ 2>/dev/null
