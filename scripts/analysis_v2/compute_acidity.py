# compute_acidity.py  ->  writes proteome_acidity.tsv
import glob, collections, os
base="/data/habib/metagenome/biosurfactant"
c2m={}
for l in open(f"{base}/integration/bgc_mag_mapping/contig_to_bin.tsv"):
    p=l.rstrip("\n").split("\t")
    if len(p)>=3: c2m[(p[0],p[1])]=f"{p[0]}_{p[2]}"
acc=collections.defaultdict(lambda:[0,0,0])
for S in ["C3","C6","H2","H6"]:
    fa=glob.glob(f"{base}/{S}/01_prodigal/*.faa")[0]
    mag=None; keep=False
    for line in open(fa):
        if line[0]==">":
            orf=line[1:].split()[0]; contig=orf.rsplit("_",1)[0]
            mag=c2m.get((S,contig)); keep=mag is not None
        elif keep:
            s=line.strip(); a=acc[mag]
            a[0]+=len(s); a[1]+=s.count("D")+s.count("E"); a[2]+=s.count("K")+s.count("R")
out=f"{base}/analysis_v2_community_abundance/proteome_acidity.tsv"
with open(out,"w") as o:
    o.write("MAG\tnres\tpctDE\tpctKR\tacidic_minus_basic\n")
    for mag,(n,de,kr) in sorted(acc.items()):
        if n: o.write(f"{mag}\t{n}\t{100*de/n:.2f}\t{100*kr/n:.2f}\t{100*(de-kr)/n:.3f}\n")
print("wrote", out)
