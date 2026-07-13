#!/bin/bash
#SBATCH -J gtdb_check                # job name
#SBATCH -p cu                        # partition (corrected!)
#SBATCH -N 1                         # one node
#SBATCH -n 1                         # one task (same as --ntasks=1)
#SBATCH -c 2                         # 2 CPU cores
#SBATCH --mem=4G                     # memory
#SBATCH -t 00:30:00                  # time (30 min; adjust if needed)
#SBATCH -o logs/%x_%j.out            # stdout
#SBATCH -e logs/%x_%j.err            # stderr

source ~/.bashrc
cd /data/habib/metagenome/biosurfactant || exit 1
mkdir -p logs

echo "===== A. GTDB-Tk reference database present & usable? ====="
echo "GTDBTK_DATA_PATH = ${GTDBTK_DATA_PATH:-<unset>}"
find /data/habib -maxdepth 4 -type d -iname "*gtdb*release*" 2>/dev/null | head
find /data/habib -maxdepth 5 \( -iname "*gtdb*.dmnd" -o -iname "bac120*.tree" -o -iname "*ar53*.tree" \) 2>/dev/null | head
find /home/habib -maxdepth 3 -type d -iname "*gtdb*" 2>/dev/null | head

echo; echo "===== B. your MAG GTDB-Tk results (free taxonomy for binned; reference for unbinned) ====="
find . -maxdepth 4 \( -iname "*bac120.summary.tsv" -o -iname "*ar53.summary.tsv" -o -iname "gtdbtk*summary.tsv" \) 2>/dev/null
f=$(find . -maxdepth 4 -iname "*summary.tsv" -path "*gtdb*" 2>/dev/null | head -1); echo ">> $f"
[ -n "$f" ] && head -2 "$f" | cut -f1,2

echo; echo "===== C. classifiers available for the unbinned residual ====="
for e in $(conda env list | awk 'NR>2 && $1!~"#"{print $1}'); do
  echo -n "env $e: "; conda run -n "$e" bash -lc 'command -v mmseqs diamond gtdbtk 2>/dev/null | tr "\n" " "' 2>/dev/null; echo
done

echo; echo "===== D. join files + inputs ====="
wc -l integration/functional/all_bgc_products.tsv integration/bgc_mag_mapping/contig_to_bin.tsv
echo "MAG fastas: $(ls 16_mag_abundance/all_MAGs_renamed/*.fa 2>/dev/null | wc -l)"
echo "coverage files: $(ls */13_coverage/coverM/*_contig_coverage.tsv 2>/dev/null | wc -l)"
echo "DONE"
