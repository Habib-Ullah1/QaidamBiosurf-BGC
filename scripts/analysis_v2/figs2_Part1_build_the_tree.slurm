#!/bin/bash
#SBATCH -J figS2_tree
#SBATCH -p cu
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --time 02:00:00
#SBATCH -o /data/habib/metagenome/biosurfactant/logs/%x-%j.out
#SBATCH -e /data/habib/metagenome/biosurfactant/logs/%x-%j.err
export LC_ALL=C.UTF-8; export LANG=C.UTF-8
source ~/.bashrc; conda activate iqtree_env
cd /data/habib/metagenome/biosurfactant/analysis_v2_community_abundance/figS2 || exit 1
REF=/data/habib/metagenome/biosurfactant/databases/napdos2_refs

# clean up the failed run's oversized files
rm -f combined_aln.faa combined.faa figS2_Cdomain.* combined_clean.faa

# combine: 20 query C-domains + 5 biosurfactant refs + 40 curated MIBiG refs = 65 tips
cat query_Cdomains.faa \
    $REF/biosurfactant_Cdomains_refs.faa \
    $REF/mibig_condensation_refs.faa > combined_clean.faa
N=$(grep -c '>' combined_clean.faa); echo "combined tips: $N"
[ "$N" -gt 400 ] && { echo "unexpectedly large ($N) — stop"; exit 1; }

# align (fast at this size) + ML tree with fixed model
mafft --thread 28 --maxiterate 1000 --localpair combined_clean.faa > combined_aln.faa
iqtree2 -s combined_aln.faa -m LG+F+G -B 1000 -T 28 -pre figS2_Cdomain --redo
echo "DONE: $(grep -c '>' combined_aln.faa) tips in tree"
