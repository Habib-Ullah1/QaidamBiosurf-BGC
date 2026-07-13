#!/usr/bin/env python3
import csv, re, statistics as st
from math import sqrt
from collections import defaultdict

B="/data/habib/metagenome/biosurfactant"
F=f"{B}/16S_data/level-6.csv"
ECON=f"{B}/integration/bgc_mag_mapping/genome_economics.tsv"

# explicit aliases for GTDB splits / synonyms between MAG genus and 16S genus
ALIAS={"Halomonas_B":["Halomonas_B","Halomonas","Billgrantia"],
       "Moraxella_A":["Moraxella_A","Moraxella"]}

# ---- load 16S: rows = samples, cols = genera (full taxonomy header) ----
rows=list(csv.reader(open(F)))
header=rows[0]
gcols=defaultdict(list)
for j,h in enumerate(header):
    if ";g__" in h:
        g=h.split(";g__")[-1].strip()
        if g and g!="__":
            gcols[g].append(j)
            gcols[re.sub(r'_[A-Z]$','',g)].append(j)   # suffix-stripped index too
samples=[r for r in rows[1:] if r and r[0] and r[0] not in ("Type","Season")]
N=len(samples)

def gvec(genus):
    cands={genus, re.sub(r'_[A-Z]$','',genus)}
    cands |= set(ALIAS.get(genus,[]))
    cols=set()
    for c in cands: cols|=set(gcols.get(c,[]))
    if not cols: return None
    v=[]
    for r in samples:
        s=0.0
        for j in cols:
            try: s+=float(r[j])
            except: pass
        v.append(s)
    return v

def levins_norm(v):                  # 1=single sample, ->1 spread evenly
    t=sum(v)
    if t<=0: return None
    p=[x/t for x in v]
    Bv=1.0/sum(pi*pi for pi in p)
    return (Bv-1)/(len(v)-1)

def spearman(x,y):
    n=len(x)
    def rk(v):
        s=sorted(range(n),key=lambda i:v[i]); r=[0]*n; i=0
        while i<n:
            j=i
            while j+1<n and v[s[j+1]]==v[s[i]]: j+=1
            for k in range(i,j+1): r[s[k]]=(i+j)/2+1
            i=j+1
        return r
    rx,ry=rk(x),rk(y); mx=sum(rx)/n; my=sum(ry)/n
    num=sum((rx[i]-mx)*(ry[i]-my) for i in range(n))
    den=sqrt(sum((rx[i]-mx)**2 for i in range(n))*sum((ry[i]-my)**2 for i in range(n)))
    rho=num/den if den else 0
    t=rho*sqrt((n-2)/(1-rho*rho)) if abs(rho)<1 and n>2 else 0
    return rho,t,n

# ---- match MAGs to 16S, compute niche metrics ----
res=[]
for m in csv.DictReader(open(ECON),delimiter='\t'):
    g=m['genus']
    if g in ("?","unknown","unclassified"): continue
    v=gvec(g)
    rec=dict(MAG=m['MAG'],genus=g,nBGC=int(m['nBGC']),
             gsize=float(m['est_gsize_mb']),season=m['season'],matched=False)
    if v and sum(v)>0:
        rec.update(matched=True, occ=sum(1 for x in v if x>0),
                   nicheB=levins_norm(v), mean16S=sum(v)/N, max16S=max(v))
    res.append(rec)

M=[r for r in res if r['matched']]
print("="*74)
print(f"NICHE-RESERVOIR TEST: is biosynthetic richness located in the rare/narrow niche?")
print("="*74)
print(f"MAGs matched to a 16S genus: {len(M)}/{len(res)}")
print(f"(unmatched are mostly GTDB placeholder lineages absent from 16S naming)\n")
if len(M)<6:
    print("  WARNING: matched n too small for correlation — interpret as anecdotal only.")

nb=[r['nBGC'] for r in M]
print("Correlations (BGC richness vs niche metrics):")
for lab,key,exp in [("occupancy (# of 17 samples)","occ","NEGATIVE if reservoir is rare"),
                    ("niche breadth (Levins, 0=specialist 1=generalist)","nicheB","NEGATIVE if specialists"),
                    ("mean 16S abundance","mean16S","NEGATIVE if rare")]:
    rho,t,n=spearman(nb,[r[key] for r in M])
    sig="**" if abs(t)>2.05 else "  "
    print(f"  vs {lab:<48} rho={rho:+.3f} t={t:+.2f} {sig}  [{exp}]")
print(f"\n  (n={len(M)}; |t|>~2.1 ~ p<0.05)")

print("\n--- matched MAGs, ranked by BGC richness ---")
print(f"  {'genus':<16}{'nBGC':>5}{'occ/17':>7}{'nicheB':>8}{'mean16S':>10}{'season':>12}")
for r in sorted(M,key=lambda x:-x['nBGC']):
    print(f"  {r['genus']:<16}{r['nBGC']:>5}{r['occ']:>7}{r['nicheB']:>8.2f}{r['mean16S']:>10.1f}{r['season']:>12}")

# quadrant: are the BGC-rich also low-occupancy (rare-biosphere reservoir)?
rich=[r for r in M if r['nBGC']>=st.median(nb)]
poor=[r for r in M if r['nBGC']<st.median(nb)]
if rich and poor:
    print(f"\n  median occupancy:  BGC-rich={st.median(r['occ'] for r in rich):.0f}/17   BGC-poor={st.median(r['occ'] for r in poor):.0f}/17")
    print(f"  median nicheB:     BGC-rich={st.median(r['nicheB'] for r in rich):.2f}      BGC-poor={st.median(r['nicheB'] for r in poor):.2f}")

with open(f"{B}/integration/bgc_mag_mapping/niche_reservoir.tsv","w") as o:
    o.write("MAG\tgenus\tnBGC\tgsize_mb\tseason\toccupancy\tnicheBreadth\tmean16S\tmax16S\n")
    for r in M:
        o.write(f"{r['MAG']}\t{r['genus']}\t{r['nBGC']}\t{r['gsize']:.2f}\t{r['season']}\t"
                f"{r['occ']}\t{r['nicheB']:.3f}\t{r['mean16S']:.2f}\t{r['max16S']:.2f}\n")
print(f"\nwritten: niche_reservoir.tsv")
