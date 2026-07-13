#!/usr/bin/env python3
import csv
from math import sqrt

B="/data/habib/metagenome/biosurfactant"
F=f"{B}/16S_data/level-6.csv"
MASTER=f"{B}/integration/bgc_mag_mapping/mag_master_table.tsv"

# ---- genus aliases: GTDB-Tk (MAG) name -> possible 16S names ----
# Billgrantia = GTDB split of Halomonas; add others as detected
ALIAS={
 "Halomonas_B":["Halomonas_B","Halomonas","Billgrantia"],
 "Moraxella_A":["Moraxella_A","Moraxella"],
 "Gillisia":["Gillisia"],
 "Loktanella":["Loktanella"],
 "Roseovarius":["Roseovarius"],
 "Salinibacter":["Salinibacter"],
 "Haloarcula":["Haloarcula"],
 "Haloferax":["Haloferax"],
 "Persicimonas":["Persicimonas"],
 "Ralstonia":["Ralstonia"],
}

# ---- load 16S: rows=samples, cols=genera ----
with open(F) as f:
    rows=list(csv.reader(f))
header=rows[0]
# map column index -> genus (last g__ token); skip index/Type/Season
genus_col={}
for j,h in enumerate(header):
    if ";g__" in h:
        g=h.split(";g__")[-1]
        if g and g!="__": genus_col[j]=g
samples={}
for r in rows[1:]:
    sid=r[0]
    if not sid: continue
    samples[sid]=r
# season by sample id prefix: C*=summer, H*=winter
def season_of(sid): return "summer" if sid.startswith("C") else "winter"

# build genus -> {sample: relabund}
def genus_abund(gnames):
    cols=[j for j,g in genus_col.items() if g in gnames]
    out={}
    for sid,r in samples.items():
        v=0.0
        for j in cols:
            try: v+=float(r[j])
            except: pass
        out[sid]=v
    return out,cols

# Mann-Whitney U (summer vs winter), returns U, direction, mean_s, mean_w
def mwu(vals):
    s=[v for sid,v in vals.items() if season_of(sid)=="summer"]
    w=[v for sid,v in vals.items() if season_of(sid)=="winter"]
    ns,nw=len(s),len(w)
    # rank all
    allv=[(v,'s') for v in s]+[(v,'w') for v in w]
    allv.sort(key=lambda x:x[0])
    ranks={}
    i=0; n=len(allv)
    rankvals=[0]*n
    while i<n:
        j=i
        while j+1<n and allv[j+1][0]==allv[i][0]: j+=1
        avg=(i+j)/2.0+1
        for k in range(i,j+1): rankvals[k]=avg
        i=j+1
    Rs=sum(rankvals[k] for k in range(n) if allv[k][1]=='s')
    Us=Rs-ns*(ns+1)/2.0
    Uw=ns*nw-Us
    U=min(Us,Uw)
    ms=sum(s)/ns if ns else 0; mw=sum(w)/nw if nw else 0
    direction="winter" if mw>ms else ("summer" if ms>mw else "tie")
    # normal approx p (n large enough-ish)
    mu=ns*nw/2.0; sig=sqrt(ns*nw*(ns+nw+1)/12.0)
    z=(U-mu)/sig if sig else 0
    return U,z,direction,ms,mw,ns,nw

# ---- load MAG master table ----
mags=[]
with open(MASTER) as f:
    rd=csv.DictReader(f,delimiter='\t')
    for row in rd: mags.append(row)

# classify MAG season from covered fraction
def mag_season(m):
    wc=max(float(m['H2cov']),float(m['H6cov']))
    sc=max(float(m['C3cov']),float(m['C6cov']))
    if wc>0.7 and sc>0.7: return "year-round"
    if wc>0.7 and sc<0.3: return "winter"
    if wc>0.7: return "winter-lean"
    if sc>0.7 and wc<0.3: return "summer"
    if sc>0.7: return "summer-lean"
    return "sporadic"

print("="*108)
print("MG vs 16S CONCORDANCE — do BGC-bearing MAG genera show matching seasonal direction in independent 16S?")
print("="*108)
print(f"{'MAG':<20}{'genus':<16}{'nBGC':>5}{'MAGseason':>12} | {'16S_dir':>8}{'16Smean_S':>10}{'16Smean_W':>10}{'MWz':>7}{'found?':>8}  concord")
print("-"*108)

# dedup genera (some MAGs share genus); report per MAG but compute 16S once per genus
seen16s={}
nconc=0; ntest=0
for m in sorted(mags,key=lambda x:-int(x['nBGC'])):
    g=m['genus']
    if g in ("?","unknown","unclassified"): 
        continue
    names=ALIAS.get(g,[g])
    if g not in seen16s:
        ab,cols=genus_abund(names)
        seen16s[g]=(ab,cols)
    ab,cols=seen16s[g]
    found = len(cols)>0 and sum(ab.values())>0
    ms=""
    magS=mag_season(m)
    if found:
        U,z,d,meanS,meanW,ns,nw=mwu(ab)
        # concordance: MAG winter-ish vs 16S winter dir
        mag_dir = "winter" if magS in ("winter","winter-lean") else ("summer" if magS in ("summer","summer-lean") else "year-round")
        if mag_dir in ("winter","summer"):
            ntest+=1
            conc = "YES" if mag_dir==d else "no"
            if conc=="YES": nconc+=1
        else:
            conc="(yr-round)"
        print(f"{m['MAG']:<20}{g:<16}{m['nBGC']:>5}{magS:>12} | {d:>8}{meanS:>10.4f}{meanW:>10.4f}{z:>7.2f}{'yes':>8}  {conc}")
    else:
        print(f"{m['MAG']:<20}{g:<16}{m['nBGC']:>5}{magS:>12} | {'--':>8}{0:>10.4f}{0:>10.4f}{0:>7.2f}{'NO':>8}  not-in-16S")

print("-"*108)
print(f"Concordance among directional MAGs: {nconc}/{ntest} agree on seasonal direction with independent 16S")
print()
print("KEY GENERA CHECK (16S companion paper reported these winter-enriched):")
for g in ["Loktanella","Moraxella_A","Halomonas_B"]:
    names=ALIAS.get(g,[g]); ab,cols=genus_abund(names)
    if sum(ab.values())>0:
        U,z,d,meanS,meanW,ns,nw=mwu(ab)
        print(f"  {g:<14} 16S direction={d:<7} summer_mean={meanS:.4f} winter_mean={meanW:.4f}  (matched-as: {','.join(names)})")
    else:
        print(f"  {g:<14} NOT FOUND in 16S under names {names}")
