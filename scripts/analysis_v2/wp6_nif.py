#!/usr/bin/env python3
# WP6 — is nitrogen fixation a lone Chloroflexota diazotroph or community-wide?
# Sum nif TPM per sample; map nif contigs to bins; report binned vs unbinned share.
import os, glob, collections
BASE="/data/habib/metagenome/biosurfactant"; S=["C3","C6","H2","H6"]
SEASON={"C3":"summer","C6":"summer","H2":"winter","H6":"winter"}
NIF={"nifH","nifD","nifK","nifE","nifN","nifB"}
c2b={}
for ln in open(f"{BASE}/integration/bgc_mag_mapping/contig_to_bin.tsv"):
    p=ln.rstrip("\n").split("\t")
    if len(p)>=3: c2b[(p[0],p[1])]=f"{p[0]}_{p[2]}"
mag_tax={}
ge=f"{BASE}/integration/bgc_mag_mapping/genome_economics.tsv"
if os.path.exists(ge):
    h=open(ge); next(h)
    for ln in h:
        p=ln.rstrip("\n").split("\t")
        if len(p)>=3: mag_tax[p[0]]=(p[1],p[2])
tpm_sum=collections.defaultdict(float); tpm_binned=collections.defaultdict(float)
persamp=collections.defaultdict(lambda:collections.defaultdict(float))
mag_nif=collections.defaultdict(collections.Counter)
for s in S:
    egg=f"{BASE}/{s}/05_eggnog/{s}_eggnog.emapper.annotations"
    nif_orf={}  # orf -> gene
    for ln in open(egg):
        if ln.startswith("#"): continue
        p=ln.rstrip("\n").split("\t")
        if len(p)<=11: continue
        if p[8] in NIF: nif_orf[p[0]]=p[8]
    tpm=f"{BASE}/integration/functional/tpm/{s}.tpm.tsv"
    T={}
    for ln in open(tpm):
        if ln.startswith("#"): continue
        p=ln.rstrip("\n").split("\t")
        if len(p)<2: continue
        try: T[p[0]]=float(p[-1])
        except ValueError: continue
    for orf,gene in nif_orf.items():
        t=T.get(orf,0.0); contig=orf.rsplit("_",1)[0]
        tpm_sum[s]+=t; persamp[s][gene]+=t
        b=c2b.get((s,contig))
        if b: tpm_binned[s]+=t; mag_nif[b][gene]+=1
print("== community nif TPM per sample ==")
for s in S: print(f"  {s} ({SEASON[s]}): total={tpm_sum[s]:.1f}  binned={tpm_binned[s]:.1f} "
                  f"({100*tpm_binned[s]/(tpm_sum[s] or 1):.1f}% binned)")
su=(tpm_sum['C3']+tpm_sum['C6'])/2; wi=(tpm_sum['H2']+tpm_sum['H6'])/2
print(f"  winter:summer nif TPM ratio = {wi/(su or 1e-9):.2f}  (manuscript claims 0.20, summer-enriched)")
print("\n== MAGs carrying nif genes (is there a complete-nifHDK diazotroph?) ==")
for b,c in sorted(mag_nif.items()):
    tax=mag_tax.get(b,("?","?"))
    complete="COMPLETE-HDK" if {"nifH","nifD","nifK"}<=set(c) else ""
    print(f"  {b:22s} {tax[1]:16s} {tax[0]:14s} {dict(c)} {complete}")
print("\nNOTE: raw physchem not on cluster (only correlation summary). Verify Table 1")
print("C6 Na=2.0 vs EC=2840 manually against the ion-chromatography sheet.")
