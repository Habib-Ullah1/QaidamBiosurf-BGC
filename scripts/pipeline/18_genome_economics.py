#!/usr/bin/env python3
import csv, glob, os
from math import sqrt

B="/data/habib/metagenome/biosurfactant"
OUT=f"{B}/integration/bgc_mag_mapping"

# ---- 1. gather CheckM2 metrics across all 4 samples ----
# Name in checkm2 = "concoct_N"; our MAG id = "{S}_concoct_N"
ck={}
for qr in glob.glob(f"{B}/*/10_checkm2/quality_report.tsv"):
    S=qr.split("/")[-3]   # C3/C6/H2/H6
    with open(qr) as f:
        r=csv.DictReader(f,delimiter='\t')
        for row in r:
            mag=f"{S}_{row['Name']}"
            try:
                ck[mag]=dict(
                    compl=float(row['Completeness']),
                    cont=float(row['Contamination']),
                    cod_dens=float(row['Coding_Density']),
                    avg_gene=float(row['Average_Gene_Length']),
                    gsize=float(row['Genome_Size']),
                    gc=float(row['GC_Content']),
                    ncds=float(row['Total_Coding_Sequences']),
                )
                # completeness-corrected estimated genome size (Mb)
                ck[mag]['est_gsize_mb']=ck[mag]['gsize']/(ck[mag]['compl']/100)/1e6 if ck[mag]['compl']>0 else None
            except: pass

# ---- 2. BGC richness + season from master table ----
master={}
with open(f"{OUT}/mag_master_table.tsv") as f:
    for m in csv.DictReader(f,delimiter='\t'):
        wc=max(float(m['H2cov']),float(m['H6cov'])); sc=max(float(m['C3cov']),float(m['C6cov']))
        if wc>0.7 and sc>0.7: sea="year-round"
        elif wc>0.7 and sc<0.3: sea="winter"
        elif wc>0.7: sea="winter-lean"
        elif sc>0.7 and wc<0.3: sea="summer"
        elif sc>0.7: sea="summer-lean"
        else: sea="sporadic"
        master[m['MAG']]=dict(nbgc=int(m['nBGC']), ncls=int(m['nClass']),
            genus=m['genus'], phylum=m['phylum'], season=sea,
            wab=(float(m['H2rel'])+float(m['H6rel']))/2,
            sab=(float(m['C3rel'])+float(m['C6rel']))/2)

# ---- 3. join, restrict to MAGs with both checkm + master ----
rows=[]
for mag in master:
    if mag in ck and ck[mag]['est_gsize_mb']:
        d=dict(MAG=mag, **master[mag], **ck[mag]); rows.append(d)

def spearman(x,y):
    n=len(x)
    def rank(v):
        s=sorted(range(n),key=lambda i:v[i]); r=[0]*n; i=0
        while i<n:
            j=i
            while j+1<n and v[s[j+1]]==v[s[i]]: j+=1
            for k in range(i,j+1): r[s[k]]=(i+j)/2.0+1
            i=j+1
        return r
    rx,ry=rank(x),rank(y); mx=sum(rx)/n; my=sum(ry)/n
    num=sum((rx[i]-mx)*(ry[i]-my) for i in range(n))
    den=sqrt(sum((rx[i]-mx)**2 for i in range(n))*sum((ry[i]-my)**2 for i in range(n)))
    rho=num/den if den else 0
    t=rho*sqrt((n-2)/(1-rho*rho)) if abs(rho)<1 else float('inf')
    return rho,t

print("="*72)
print(f"GENOME ECONOMICS  (n={len(rows)} MAGs with CheckM2 + abundance data)")
print("="*72)
print("Hypothesis: BGC-rich MAGs = larger, less-streamlined genomes")
print("(if true -> biosynthesis is a costly trait of non-streamlined specialists)\n")

nbgc=[r['nbgc'] for r in rows]
for label,key,pred in [
    ("estimated genome size (Mb)","est_gsize_mb","POSITIVE expected"),
    ("coding density","cod_dens","NEGATIVE expected (streamlined=dense)"),
    ("avg gene length","avg_gene","~"),
    ("GC content","gc","~"),
    ("total CDS count","ncds","POSITIVE expected"),
]:
    y=[r[key] for r in rows]
    rho,t=spearman(nbgc,y)
    sig="**" if abs(t)>2.02 else "  "
    print(f"  BGC richness vs {label:<26} rho={rho:+.3f} t={t:+.2f} {sig}  [{pred}]")

print("\n  (n=%d, df=%d, |t|>2.02 = p<0.05 = significant **)"%(len(rows),len(rows)-2))

# ---- 4. by cohort: are winter specialists larger-genomed? ----
print("\n--- genome size & streamlining by seasonal cohort ---")
import statistics as st
groups={}
for r in rows:
    g="winter" if r['season'] in ("winter","winter-lean") else ("summer" if r['season'] in ("summer","summer-lean") else "year-round")
    groups.setdefault(g,[]).append(r)
print(f"  {'cohort':<12}{'n':>3}{'med_BGC':>9}{'med_Gsize_Mb':>14}{'med_codDens':>12}")
for g in ("winter","year-round","summer"):
    if g in groups:
        G=groups[g]
        print(f"  {g:<12}{len(G):>3}{st.median(r['nbgc'] for r in G):>9.0f}"
              f"{st.median(r['est_gsize_mb'] for r in G):>14.2f}{st.median(r['cod_dens'] for r in G):>12.3f}")

# ---- 5. the dominant vs specialist contrast (the money table) ----
print("\n--- TOP 5 most abundant (winter) vs TOP 5 most BGC-rich ---")
print("  [most abundant winter:]")
for r in sorted(rows,key=lambda x:-x['wab'])[:5]:
    print(f"    {r['genus']:<16} ab={r['wab']:.2f}%  BGC={r['nbgc']:>2}  Gsize={r['est_gsize_mb']:.2f}Mb  codDens={r['cod_dens']:.3f}")
print("  [most BGC-rich:]")
for r in sorted(rows,key=lambda x:-x['nbgc'])[:5]:
    print(f"    {r['genus']:<16} ab={r['wab']:.2f}%  BGC={r['nbgc']:>2}  Gsize={r['est_gsize_mb']:.2f}Mb  codDens={r['cod_dens']:.3f}")

# write master economics table
with open(f"{OUT}/genome_economics.tsv","w") as o:
    o.write("MAG\tgenus\tphylum\tseason\tnBGC\tnClass\twinterAb\tsummerAb\test_gsize_mb\tcoding_density\tavg_gene_len\tGC\tnCDS\tcompleteness\n")
    for r in rows:
        o.write(f"{r['MAG']}\t{r['genus']}\t{r['phylum']}\t{r['season']}\t{r['nbgc']}\t{r['ncls']}\t"
                f"{r['wab']:.3f}\t{r['sab']:.3f}\t{r['est_gsize_mb']:.3f}\t{r['cod_dens']:.3f}\t"
                f"{r['avg_gene']:.1f}\t{r['gc']:.3f}\t{r['ncds']:.0f}\t{r['compl']:.1f}\n")
print(f"\nwritten: {OUT}/genome_economics.tsv")
