#!/bin/bash
#SBATCH --job-name=bigscape_mibig
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --output=/data/habib/metagenome/biosurfactant/logs/bigscape_mibig_%j.log

source ~/.bashrc
conda activate bigscape_env

PROJECT=/data/habib/metagenome
PFAM=$PROJECT/biosurfactant/databases/pfam
BGCDIR=$PROJECT/biosurfactant/bigscape_input_3samples
OUTDIR=$PROJECT/biosurfactant/bigscape_output_all4_mibig2

echo "=== BiG-SCAPE all4 + MIBiG  $(date) ==="
echo "Input BGCs: $(ls $BGCDIR/*.gbk | wc -l)"
echo "Output: $OUTDIR"

bigscape \
    --inputdir $BGCDIR \
    --outputdir $OUTDIR \
    --pfam_dir $PFAM \
    --cores 16 \
    --cutoffs 0.3 0.5 0.7 \
    --mix \
    --mibig \
    --include_gbk_str region \
    --clans-off \
    --verbose

echo "=== done $(date) ==="
echo "--- MIBiG refs now present? ---"
grep -c "^BGC0" $OUTDIR/network_files/*/Network_Annotations_Full.tsv 2>/dev/null
