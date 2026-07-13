#!/bin/bash
#SBATCH --job-name=sandpuma_all
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=60G
#SBATCH --time=24:00:00
#SBATCH --output=/data/habib/metagenome/biosurfactant/logs/sandpuma_%j.log

source ~/.bashrc
conda activate antismash

PROJECT=/data/habib/metagenome

for s in C3 C6 H2 H6; do
    echo "=== Extracting A-domains: $s $(date) ==="
    ASDIR=$PROJECT/biosurfactant/$s/02_antismash
    OUTDIR=$PROJECT/biosurfactant/$s/06_sandpuma
    mkdir -p $OUTDIR

    # Extract A-domain sequences from antiSMASH genbank files
    python3 << INNEREOF
import glob, os
from Bio import SeqIO

asdir = "$ASDIR"
outdir = "$OUTDIR"
a_seqs = []

for gbk in sorted(glob.glob(f"{asdir}/*.region*.gbk")):
    for rec in SeqIO.parse(gbk, "genbank"):
        for feat in rec.features:
            if feat.type != "aSDomain":
                continue
            domain = feat.qualifiers.get("aSDomain", [""])[0]
            if "AMP-binding" not in domain and "A-OX" not in domain:
                continue
            locus = feat.qualifiers.get("locus_tag", ["unknown"])[0]
            domain_id = feat.qualifiers.get("domain_id", [locus])[0]
            try:
                seq = feat.extract(rec.seq).translate(to_stop=True)
                if len(seq) > 50:
                    a_seqs.append(f">{domain_id}\n{seq}")
            except Exception as e:
                pass

outfa = f"{outdir}/${s}_Adomains.faa"
with open(outfa, "w") as fh:
    fh.write("\n".join(a_seqs) + "\n")
print(f"$s: {len(a_seqs)} A-domain sequences -> {outfa}")
INNEREOF

    # Count extracted sequences
    N=$(grep -c "^>" $OUTDIR/${s}_Adomains.faa 2>/dev/null || echo 0)
    echo "  A-domains extracted: $N"

    # Run SANDPUMA if sequences exist and environment available
    if [ "$N" -gt 0 ]; then
        if conda activate sandpuma_env 2>/dev/null; then
            echo "  Running SANDPUMA on $N sequences..."
            sandpuma.py \
                --input $OUTDIR/${s}_Adomains.faa \
                --output $OUTDIR/${s}_sandpuma_results \
                --threads 16 \
                2> $PROJECT/biosurfactant/logs/sandpuma_${s}.log || \
            echo "  SANDPUMA failed — check log"
            conda activate antismash
        else
            echo "  sandpuma_env not available — using antiSMASH ASF substrate predictions"
            # Extract substrate predictions from antiSMASH output directly
            grep -h "substrate" $ASDIR/*.gbk 2>/dev/null | \
                grep "specificity" | \
                sort | uniq -c | sort -rn | head -30 \
                > $OUTDIR/${s}_substrate_from_antismash.txt
            echo "  Substrate predictions extracted from antiSMASH ASF output"
        fi
    fi
    echo "=== Done: $s $(date) ==="
done

echo "=== All SANDPUMA complete $(date) ==="
