#!/bin/bash
#SBATCH --job-name=gene_linkage
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=40G
#SBATCH --output=/data/habib/metagenome/biosurfactant/logs/gene_linkage_%j.log

source ~/.bashrc
conda activate biosurfactant_base

PROJECT=/data/habib/metagenome
INTEG=$PROJECT/biosurfactant/integration
mkdir -p $INTEG

python3 << INNEREOF
import csv, os, glob

PROJECT = "/data/habib/metagenome"
INTEG = f"{PROJECT}/biosurfactant/integration"
samples = ["C3", "C6", "H2", "H6"]

# Load GTDB-Tk taxonomy for all samples
taxonomy = {}  # bin_id -> taxonomy string
for s in samples:
    for f in [
        f"{PROJECT}/biosurfactant/{s}/11_gtdbtk/output/gtdbtk.bac120.summary.tsv",
        f"{PROJECT}/biosurfactant/{s}/11_gtdbtk/output/gtdbtk.ar53.summary.tsv"
    ]:
        if os.path.exists(f):
            with open(f) as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    taxonomy[f"{s}|{row['user_genome']}"] = row.get('classification', 'Unclassified')

print(f"Loaded taxonomy for {len(taxonomy)} MAGs")

# Load CheckM2 quality for all bins
quality = {}  # sample|bin -> (completeness, contamination)
for s in samples:
    f = f"{PROJECT}/biosurfactant/{s}/10_checkm2/quality_report.tsv"
    if os.path.exists(f):
        with open(f) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                quality[f"{s}|{row['Name']}"] = (row['Completeness'], row['Contamination'])

print(f"Loaded quality for {len(quality)} bins")

# Build contig -> bin mapping from metawrap refined bins
contig_to_bin = {}  # sample|contig -> bin_id
for s in samples:
    bins_dir = f"{PROJECT}/biosurfactant/{s}/09_binning/metawrap_refined/bins"
    for fa in glob.glob(f"{bins_dir}/*.fa"):
        bin_id = os.path.basename(fa).replace(".fa", "")
        with open(fa) as fh:
            for line in fh:
                if line.startswith(">"):
                    contig = line[1:].strip().split()[0]
                    contig_to_bin[f"{s}|{contig}"] = bin_id

print(f"Mapped {len(contig_to_bin)} contigs to bins")

# Load DIAMOND hits and link to bins + taxonomy
outrows = []
for s in samples:
    diamond_f = f"{PROJECT}/biosurfactant/{s}/04_diamond/{s}_diamond.tsv"
    if not os.path.exists(diamond_f):
        print(f"WARNING: DIAMOND file missing for {s}")
        continue
    with open(diamond_f) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            orf_id = parts[0]
            # Extract contig ID from ORF ID (Prodigal format: contig_orfnum)
            contig = "_".join(orf_id.split("_")[:-1])
            key = f"{s}|{contig}"
            bin_id = contig_to_bin.get(key, "unbinned")
            tax_key = f"{s}|{bin_id}"
            taxon = taxonomy.get(tax_key, "Unclassified")
            comp, cont = quality.get(tax_key, ("NA", "NA"))
            # Extract phylum from GTDB taxonomy string
            phylum = "Unknown"
            if "p__" in taxon:
                try:
                    phylum = taxon.split("p__")[1].split(";")[0]
                except:
                    pass
            outrows.append({
                "sample": s,
                "orf_id": orf_id,
                "contig_id": contig,
                "bin_id": bin_id,
                "gtdbtk_taxon": taxon,
                "phylum": phylum,
                "bin_quality": f"comp={comp},cont={cont}"
            })

outfile = f"{INTEG}/gene_to_organism_linkage.tsv"
fields = ["sample","orf_id","contig_id","bin_id","gtdbtk_taxon","phylum","bin_quality"]
with open(outfile, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
    w.writeheader()
    w.writerows(outrows)

print(f"Written {len(outrows)} gene-organism links to {outfile}")
INNEREOF
