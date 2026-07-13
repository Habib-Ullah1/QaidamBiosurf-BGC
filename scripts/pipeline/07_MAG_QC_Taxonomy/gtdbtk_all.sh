#!/bin/bash
#SBATCH --job-name=gtdbtk_all
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=200G
#SBATCH --output=/data/habib/metagenome/biosurfactant/logs/gtdbtk_all_%j.log

source ~/.bashrc
conda activate gtdbtk_env

PROJECT=/data/habib/metagenome

for s in C3 C6 H2 H6; do
    echo "=== Filtering bins: $s $(date) ==="

    BINS=$PROJECT/biosurfactant/$s/09_binning/metawrap_refined/bins
    FILTERED=$PROJECT/biosurfactant/$s/11_gtdbtk/filtered_bins
    OUTDIR=$PROJECT/biosurfactant/$s/11_gtdbtk/output
    mkdir -p $FILTERED $OUTDIR

    # Clear previous filtered bins
    rm -f $FILTERED/*.fa

    python3 << INNEREOF
import csv, os, shutil
qreport = "$PROJECT/biosurfactant/$s/10_checkm2/quality_report.tsv"
bins_dir = "$BINS"
out_dir = "$FILTERED"
kept = 0
with open(qreport) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        comp = float(row["Completeness"])
        cont = float(row["Contamination"])
        if comp >= 50.0 and cont <= 10.0:
            src = os.path.join(bins_dir, row["Name"] + ".fa")
            dst = os.path.join(out_dir, row["Name"] + ".fa")
            if os.path.exists(src):
                shutil.copy2(src, dst)
                kept += 1
print(f"$s: {kept} bins passed filter (comp>=50, cont<=10)")
INNEREOF

    N=$(ls $FILTERED/*.fa 2>/dev/null | wc -l)
    echo "  Bins for GTDB-Tk: $N"

    if [ "$N" -gt 0 ]; then
        echo "=== GTDB-Tk classify_wf: $s $(date) ==="
        gtdbtk classify_wf \
            --genome_dir $FILTERED \
            --out_dir $OUTDIR \
            --extension fa \
            --cpus 16 \
            --skip_ani_screen
        echo "=== GTDB-Tk done: $s $(date) ==="
    else
        echo "  WARNING: No bins passed filter for $s — skipping"
    fi
done

echo "=== All GTDB-Tk complete $(date) ==="
