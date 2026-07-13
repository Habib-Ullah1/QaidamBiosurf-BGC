#!/usr/bin/env python3
import csv, sys
from math import sqrt

B="/data/habib/metagenome/biosurfactant"
OUT=f"{B}/integration/bgc_mag_mapping"
COV=f"{B}/16_mag_abundance/MAG_relabund.tsv"

# ---- Load CoverM (explicit column indices, verified) ----
# 0:Genome 1:C3rel 2:C3mean 3:C3cov 4:C6rel 5:C6mean 6:C6cov
# 7:H2rel 8:H2mean 9:H2cov 10:H6rel 11:H6mean 12:H6cov
cov={}
with open(COV) as f:
    r=csv.reader(f,delimiter='\t'); next(r)
    for row in r:
        if row[0] in ("unmapped","unbinned"): continue
        g=row[0]
        cov[g]=dict(
            C3rel=float(row[1]),  C3mean=float(row[2]),  C3cov=float(row[3]),
            C6rel=float(row[4]),  C6mean=float(row[5]),  C6cov=float(row[6]),
            H2rel=float(row[7]),  H2mean=float(row[8]),  H2cov=float(row[9]),
            H6rel=float(row[10]), H6mean=float(row[11]), H6cov=float(row[12]),
        )

# ---- Load BGC richness per MAG ----
bgc={}
with open(f"{OUT}/mag_bgc_summary.tsv") as f:
    for row in csv.reader(f,delimiter='\t'):
        # mag, phylum, genus, totalBGC, nclass, types
        bgc[row[0]]=dict(phylum=row[1],genus=row[2],nbgc=int(row[3]),ncls=int(row[4]))

# Every MAG gets a BGC count (0 if no BGCs)
for g in cov:
    if g not in bgc:
        bgc[g]=dict(phylum="?",genus="?",nbgc=0,ncls=0)

mags=sorted(cov.keys())

print("="*78)
print("TEST 2: Is 'winter-only' true absence or assembly artifact?")
print("="*78)
print("Logic: classify each MAG by summer vs winter covered_fraction AND mean depth.")
print("  - winter present  = max(H2cov,H6cov) > 0.7")
print("  - summer present  = max(C3cov,C6cov) > 0.7")
print("  - 'true absence' in summer requires summer reads DID map elsewhere")
print("    (we have summer data; low cov + low depth = genuine scarcity, not missing data)")
print()
hdr=f"{'MAG':<22}{'genus':<16}{'BGC':>4} | {'C3cov':>6}{'C6cov':>6}{'H2cov':>6}{'H6cov':>6} | {'C3dep':>7}{'C6dep':>7}{'H2dep':>7}{'H6dep':>7} | class"
print(hdr); print("-"*len(hdr))

cats={}
rows=[]
for g in mags:
    c=cov[g]
    wcov=max(c['H2cov'],c['H6cov']); scov=max(c['C3cov'],c['C6cov'])
    wdep=max(c['H2mean'],c['H6mean']); sdep=max(c['C3mean'],c['C6mean'])
    if wcov>0.7 and scov>0.7: cat="year-round"
    elif wcov>0.7 and scov<0.3: cat="WINTER-true"   # clearly present winter, clearly absent summer
    elif wcov>0.7 and scov<=0.7: cat="winter-lean"  # winter present, summer marginal
    elif scov>0.7 and wcov<0.3: cat="summer-true"
    elif scov>0.7: cat="summer-lean"
    else: cat="sporadic"
    cats[cat]=cats.get(cat,0)+1
    rows.append((c,g,cat,wcov,scov))

# sort: winter-true first then by winter depth
order={"WINTER-true":0,"winter-lean":1,"year-round":2,"summer-lean":3,"summer-true":4,"sporadic":5}
for c,g,cat,wcov,scov in sorted(rows,key=lambda x:(order[x[2]],-max(x[0]['H2mean'],x[0]['H6mean']))):
    print(f"{g:<22}{bgc[g]['genus']:<16}{bgc[g]['nbgc']:>4} | "
          f"{c['C3cov']:>6.2f}{c['C6cov']:>6.2f}{c['H2cov']:>6.2f}{c['H6cov']:>6.2f} | "
          f"{c['C3mean']:>7.1f}{c['C6mean']:>7.1f}{c['H2mean']:>7.1f}{c['H6mean']:>7.1f} | {cat}")

print()
print("Category counts:")
for k in sorted(cats,key=lambda x:order[x]):
    print(f"  {k:<14} {cats[k]:>2}")
print()
print("INTERPRETATION GUIDE:")
print("  WINTER-true MAGs with summer mean-depth NEAR ZERO = genuine winter specialists")
print("  (summer reads were sequenced & mapped to OTHER genomes, so absence is real).")
print("  If a WINTER-true MAG shows summer mean-depth > a few X, treat with caution.")

print()
print("="*78)
print("TEST 1: Spearman correlation — BGC richness vs winter abundance (all 44 MAGs)")
print("="*78)

def spearman(x,y):
    n=len(x)
    def ranks(v):
        s=sorted(range(n),key=lambda i:v[i])
        r=[0]*n; i=0
        while i<n:
            j=i
            while j+1<n and v[s[j+1]]==v[s[i]]: j+=1
            avg=(i+j)/2.0+1
            for k in range(i,j+1): r[s[k]]=avg
            i=j+1
        return r
    rx,ry=ranks(x),ranks(y)
    mx=sum(rx)/n; my=sum(ry)/n
    num=sum((rx[i]-mx)*(ry[i]-my) for i in range(n))
    den=sqrt(sum((rx[i]-mx)**2 for i in range(n))*sum((ry[i]-my)**2 for i in range(n)))
    rho=num/den if den else 0
    # t-approx p-value
    if abs(rho)<1:
        t=rho*sqrt((n-2)/(1-rho*rho))
    else:
        t=float('inf')
    return rho,t,n

nbgc=[bgc[g]['nbgc'] for g in mags]
wab=[(cov[g]['H2rel']+cov[g]['H6rel'])/2 for g in mags]
sab=[(cov[g]['C3rel']+cov[g]['C6rel'])/2 for g in mags]
ncls=[bgc[g]['ncls'] for g in mags]

for label,yv in [("winter abundance",wab),("summer abundance",sab),("total abundance",[wab[i]+sab[i] for i in range(len(mags))])]:
    rho,t,n=spearman(nbgc,yv)
    print(f"  BGC count vs {label:<18}: rho={rho:+.3f}  t={t:+.2f}  n={n}")
rho,t,n=spearman(ncls,wab)
print(f"  BGC CLASS count vs winter abundance : rho={rho:+.3f}  t={t:+.2f}  n={n}")

print()
print("  Critical |t| at n=44 (df=42): ~2.02 (p=0.05).  |t|<2.02 => NOT significant.")
print()
print("INTERPRETATION:")
print("  rho near 0 / non-sig => richness and abundance are DECOUPLED")
print("  (supports the 'two strategies' framing, NOT 'dominant=richest').")

# ---- write clean table for figures ----
with open(f"{OUT}/mag_master_table.tsv","w") as o:
    o.write("MAG\tphylum\tgenus\tnBGC\tnClass\tC3rel\tC6rel\tH2rel\tH6rel\tC3cov\tC6cov\tH2cov\tH6cov\twinterAbund\tsummerAbund\n")
    for g in mags:
        c=cov[g]
        o.write(f"{g}\t{bgc[g]['phylum']}\t{bgc[g]['genus']}\t{bgc[g]['nbgc']}\t{bgc[g]['ncls']}\t"
                f"{c['C3rel']:.4f}\t{c['C6rel']:.4f}\t{c['H2rel']:.4f}\t{c['H6rel']:.4f}\t"
                f"{c['C3cov']:.4f}\t{c['C6cov']:.4f}\t{c['H2cov']:.4f}\t{c['H6cov']:.4f}\t"
                f"{(c['H2rel']+c['H6rel'])/2:.4f}\t{(c['C3rel']+c['C6rel'])/2:.4f}\n")
print()
print(f"Master table written: {OUT}/mag_master_table.tsv (44 MAGs, clean, for figures)")
