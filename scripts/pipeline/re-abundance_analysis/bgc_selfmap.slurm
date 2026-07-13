#!/bin/bash
#SBATCH -J bgc_selfmap
#SBATCH -p cu
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --time 04:00:00
#SBATCH -o /data/habib/metagenome/biosurfactant/logs/%x-%j.out
#SBATCH -e /data/habib/metagenome/biosurfactant/logs/%x-%j.err

source ~/.bashrc
conda activate biosurfactant_base
cd /data/habib/metagenome/biosurfactant || exit 1
mkdir -p logs integration/bgc_taxonomy
OUT=integration/bgc_taxonomy
THREADS=28

# 1) MAG reference proteins — pulled DIRECTLY from existing per-sample prodigal .faa
#    (no re-running prodigal). contig_to_bin gives sample|contig -> bin; tag each protein header with its MAG.
python3 - <<'PY'
import collections, glob
base="/data/habib/metagenome/biosurfactant"
c2b={}
for l in open(f"{base}/integration/bgc_mag_mapping/contig_to_bin.tsv"):
    p=l.rstrip("\n").split("\t")
    if len(p)>=3: c2b[(p[0],p[1])]=f"{p[0]}_{p[2]}"     # -> MAG name e.g. C3_concoct_111
n=0
with open(f"{base}/integration/bgc_taxonomy/mag_proteins.faa","w") as o:
    for S in ["C3","C6","H2","H6"]:
        faa=glob.glob(f"{base}/{S}/01_prodigal/*.faa")[0]
        keep=False; mag=None
        for line in open(faa):
            if line[0]==">":
                orf=line[1:].split()[0]; contig=orf.rsplit("_",1)[0]
                mag=c2b.get((S,contig)); keep=mag is not None
                if keep: o.write(f">{mag}|{S}|{orf}\n"); n+=1
            elif keep: o.write(line)
print(f"MAG reference proteins written: {n}")
PY
# fail loudly if empty, instead of cascading
if [ ! -s $OUT/mag_proteins.faa ]; then echo "ERROR: mag_proteins.faa is empty — check contig_to_bin / prodigal paths"; exit 1; fi
diamond makedb --in $OUT/mag_proteins.faa -d $OUT/mag_ref --threads $THREADS

# 2) unbinned BGC-contig proteins (this half already worked)
python3 - <<'PY'
import collections, glob
base="/data/habib/metagenome/biosurfactant"
want=collections.defaultdict(set)
for l in open(f"{base}/integration/functional/all_bgc_products.tsv"):
    x=l.rstrip("\n").split("\t")
    if (x[3] if len(x)>3 else "UNBINNED")=="UNBINNED": want[x[0]].add(x[1])
with open(f"{base}/integration/bgc_taxonomy/unbinned_bgc_proteins.faa","w") as o:
    for S in ["C3","C6","H2","H6"]:
        faa=glob.glob(f"{base}/{S}/01_prodigal/*.faa")[0]
        keep=False
        for line in open(faa):
            if line[0]==">":
                orf=line[1:].split()[0]; contig=orf.rsplit("_",1)[0]
                keep=contig in want[S]
                if keep: o.write(f">{S}|{orf}\n")
            elif keep: o.write(line)
print("unbinned BGC proteins written")
PY

# 3) DIAMOND: unbinned BGC proteins vs MAG reference
diamond blastp -q $OUT/unbinned_bgc_proteins.faa -d $OUT/mag_ref \
  --outfmt 6 qseqid sseqid pident bitscore --max-target-seqs 1 --evalue 1e-5 \
  --id 50 --query-cover 50 --threads $THREADS -o $OUT/unbinned_vs_mags.tsv
echo "DONE: $(wc -l < $OUT/unbinned_vs_mags.tsv) hits"
