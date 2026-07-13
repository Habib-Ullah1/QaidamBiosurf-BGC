#!/usr/bin/env python3
import csv, glob, os, sys

B="/data/habib/metagenome/biosurfactant"
OUT=f"{B}/integration/bgc_mag_mapping"
NET=glob.glob(f"{B}/bigscape_output_all4_mibig2/network_files/*/mix/mix_c0.50.network")[0]
CUTOFF=0.50

def is_mibig(x): return x.startswith("BGC0")
def is_ours(x):  return x.split("_")[0] in ("C3","C6","H2","H6")

# ---- collect all our BGCs and which have a MIBiG edge < cutoff ----
ours=set(); mibig_linked=set()
with open(NET) as f:
    r=csv.reader(f,delimiter='\t'); next(r)
    for row in r:
        a,b=row[0],row[1]
        try: d=float(row[2])
        except: continue
        for x in (a,b):
            if is_ours(x): ours.add(x)
        if d<CUTOFF:
            if is_ours(a) and is_mibig(b): mibig_linked.add(a)
            if is_ours(b) and is_mibig(a): mibig_linked.add(b)

# BGCs only appear in .network if they have >=1 edge; singletons (no edge) are absent.
# Pull the FULL BGC list from antiSMASH region mapping so singletons count as novel.
allbgc={}  # region_id -> (sample, contig)
with open(f"{OUT}/bgc_regions.tsv") as f:
    for row in csv.reader(f,delimiter='\t'):
        s,contig,region,prod=row[0],row[1],row[2],row[3]
        rid=f"{s}_{region}"
        allbgc[rid]=(s,contig,prod)

total=len(allbgc)
linked=sum(1 for b in allbgc if b in mibig_linked)
novel=total-linked
print("="*70)
print("BGC NOVELTY vs MIBiG 3.1 (BiG-SCAPE network, distance < 0.50)")
print("="*70)
print(f"  Total BGCs:            {total}")
print(f"  MIBiG-similar (<0.50): {linked}  ({100*linked/total:.1f}%)")
print(f"  Novel (no MIBiG link): {novel}  ({100*novel/total:.1f}%)")

# ---- novelty per BGC class ----
print("\n--- novelty by BGC product class ---")
cls_tot={}; cls_novel={}
for rid,(s,contig,prod) in allbgc.items():
    cls_tot[prod]=cls_tot.get(prod,0)+1
    if rid not in mibig_linked: cls_novel[prod]=cls_novel.get(prod,0)+1
print(f"  {'class':<22}{'total':>6}{'novel':>6}{'%novel':>8}")
for c in sorted(cls_tot,key=lambda x:-cls_tot[x]):
    n=cls_novel.get(c,0); print(f"  {c:<22}{cls_tot[c]:>6}{n:>6}{100*n/cls_tot[c]:>7.0f}%")

# ---- novelty by seasonal cohort (join BGC->MAG->season) ----
# load BGC->bin
binmap={}
with open(f"{OUT}/bgc_mag_taxonomy_v2.tsv") as f:
    for row in csv.reader(f,delimiter='\t'):
        s,contig,region,prod,bn,tax=row[0],row[1],row[2],row[3],row[4],row[5]
        binmap[f"{s}_{region}"]=f"{s}_{bn}" if bn!="unbinned" else "unbinned"
# load MAG seasonal class
season={}
with open(f"{OUT}/mag_master_table.tsv") as f:
    rd=csv.DictReader(f,delimiter='\t')
    for m in rd:
        wc=max(float(m['H2cov']),float(m['H6cov'])); sc=max(float(m['C3cov']),float(m['C6cov']))
        if wc>0.7 and sc>0.7: season[m['MAG']]="year-round"
        elif wc>0.7 and sc<0.3: season[m['MAG']]="winter"
        elif wc>0.7: season[m['MAG']]="winter-lean"
        elif sc>0.7 and wc<0.3: season[m['MAG']]="summer"
        elif sc>0.7: season[m['MAG']]="summer-lean"
        else: season[m['MAG']]="sporadic"

print("\n--- novelty by seasonal cohort (binned BGCs only) ---")
coh_tot={}; coh_novel={}
WINTER={"winter","winter-lean"}
for rid,(s,contig,prod) in allbgc.items():
    b=binmap.get(rid,"unbinned")
    if b=="unbinned": continue
    sea=season.get(b,"?")
    grp = "winter-cohort" if sea in WINTER else ("summer-cohort" if sea in ("summer","summer-lean") else "year-round")
    coh_tot[grp]=coh_tot.get(grp,0)+1
    if rid not in mibig_linked: coh_novel[grp]=coh_novel.get(grp,0)+1
print(f"  {'cohort':<16}{'total':>6}{'novel':>6}{'%novel':>8}")
for g in ("winter-cohort","year-round","summer-cohort"):
    if g in coh_tot:
        n=coh_novel.get(g,0); print(f"  {g:<16}{coh_tot[g]:>6}{n:>6}{100*n/coh_tot[g]:>7.0f}%")

# ---- the mechanistic cut: NRPS novelty in winter cohort ----
print("\n--- KEY: NRPS/lipopeptide-class BGCs in winter cohort ---")
NRPSish=lambda p: any(k in p for k in ("NRPS","NAPAA","RiPP","lasso","thio"))
wt=wn=0
for rid,(s,contig,prod) in allbgc.items():
    b=binmap.get(rid,"unbinned")
    if b=="unbinned": continue
    if season.get(b,"?") not in WINTER: continue
    if not NRPSish(prod): continue
    wt+=1
    if rid not in mibig_linked: wn+=1
print(f"  winter-cohort NRPS-type BGCs: {wt}, novel: {wn} ({100*wn/wt:.0f}%)" if wt else "  (none)")

with open(f"{OUT}/novelty_summary.txt","w") as o:
    o.write(f"total={total} linked={linked} novel={novel} pct_novel={100*novel/total:.1f}\n")
print(f"\nwritten: {OUT}/novelty_summary.txt")
