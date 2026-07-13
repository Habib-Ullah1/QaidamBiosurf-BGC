#!/usr/bin/env python3
# WP1 part1 — compositional-artifact control: build KO x sample count matrices
# (whole community + bacterial-only) and a per-function per-sample table.
import os, glob, re, sys, collections

BASE    = "/data/habib/metagenome/biosurfactant"
SAMPLES = ["C3","C6","H2","H6"]
SEASON  = {"C3":"summer","C6":"summer","H2":"winter","H6":"winter"}
OUTDIR  = f"{BASE}/integration/wp1_compositional"; os.makedirs(OUTDIR, exist_ok=True)
COUNTS  = "{S}/../integration/functional/tpm/{S}.counts.tsv"   # confirmed path
EGG     = "{S}/05_eggnog/{S}_eggnog.emapper.annotations"
Q,PREF,OGS,KO = 0,8,4,11          # 0-based eggNOG columns; PFAMs = last col

# target functions: match on KO set OR preferred-name regex OR PFAM substring
TARGETS = {
 "phospholipaseA":   {"ko":{"K01058","K01059"},"name":r"^(pldA|plaA|plb|pldB)$","pfam":None},
 "lysophospholipase":{"ko":{"K01048","K14676","K06999"},"name":r"^(lypA|tesA)$","pfam":None},
 "FA_desaturase":    {"ko":set(),"name":r"^(desA|desB|desC)$","pfam":"FA_desaturase"},
 "cardiolipin_cfa":  {"ko":{"K06131","K06132","K00574"},"name":r"^(cls|clsA|clsB|cfa)$","pfam":None},
 "wza":              {"ko":{"K01991"},"name":r"^wza$","pfam":None},
 "wzc":              {"ko":{"K16692","K16081"},"name":r"^wzc$","pfam":None},
 "kps_group2":       {"ko":set(),"name":r"^(kpsE|kpsT|kpsM|kpsD|kpsC|kpsS)$","pfam":None},
 "wzx_wzy_control":  {"ko":set(),"name":r"^(wzx|wzy)$","pfam":None},
 "ectoine":          {"ko":{"K06718","K06720","K00836","K10674"},"name":r"^(ectA|ectB|ectC|ectD)$","pfam":None},
 "glycine_betaine":  {"ko":{"K00108","K00130","K06904"},"name":r"^(betA|betB|gbsA)$","pfam":None},
 "trehalose":        {"ko":{"K00697","K01087"},"name":r"^(otsA|otsB)$","pfam":None},
 "K_uptake":         {"ko":{"K03498","K03499","K01546","K01547","K01548","K03305"},"name":r"^(trkA|trkH|kdpA|kdpB|kdpC|ktrA|ktrB)$","pfam":None},
 "rhodopsin":        {"ko":{"K04641"},"name":r"^(bop|hop|xop)$","pfam":None},
}
def matches(kos,pref,pfam):
    out=[]
    for g,r in TARGETS.items():
        ok = bool(r["ko"] & set(kos))
        if not ok and r["name"] and pref not in ("","-") and re.match(r["name"],pref,re.I): ok=True
        if not ok and r["pfam"] and r["pfam"].lower() in (pfam or "").lower(): ok=True
        if ok: out.append(g)
    return out
def contig_of(orf): return orf.rsplit("_",1)[0]

ko_all=collections.defaultdict(lambda:collections.defaultdict(int))
ko_bact=collections.defaultdict(lambda:collections.defaultdict(int))
tgt_all=collections.defaultdict(lambda:collections.defaultdict(int))
tgt_bact=collections.defaultdict(lambda:collections.defaultdict(int))
lib_all=collections.defaultdict(int); lib_bact=collections.defaultdict(int)
domsum=collections.defaultdict(collections.Counter)

for S in SAMPLES:
    cf=os.path.normpath(os.path.join(BASE,COUNTS.format(S=S)))
    ef=os.path.join(BASE,EGG.format(S=S))
    if not os.path.exists(cf): sys.exit(f"[FATAL] no counts {cf}")
    # counts: orf -> int (featureCounts: skip #; take col0 & last col)
    counts={}
    for ln in open(cf):
        if ln.startswith("#"): continue
        p=ln.rstrip("\n").split("\t")
        if len(p)<2 or p[0]=="Geneid": continue
        try: c=int(round(float(p[-1])))
        except ValueError: continue
        if c>0: counts[p[0]]=c
    # pass 1: contig domain vote from eggNOG_OGs
    cdom=collections.defaultdict(lambda:[0,0,0])   # [arch,bact,unk]
    for ln in open(ef):
        if ln.startswith("#"): continue
        p=ln.rstrip("\n").split("\t")
        if len(p)<=KO: continue
        ogs=p[OGS]; a=("|Archaea" in ogs); b=("|Bacteria" in ogs)
        d=0 if (a and not b) else (1 if b else 2)
        cdom[contig_of(p[Q])][d]+=1
    dom={c:("A" if v[0]>v[1] and v[0]>=v[2] else "B") for c,v in cdom.items()}  # unknown->B (conservative)
    del cdom
    # pass 2: accumulate
    for ln in open(ef):
        if ln.startswith("#"): continue
        p=ln.rstrip("\n").split("\t")
        if len(p)<=KO: continue
        orf=p[Q]; c=counts.get(orf,0)
        if c<=0: continue
        kos=[k.replace("ko:","") for k in p[KO].split(",") if k not in ("","-")]
        pref=p[PREF]; pfam=p[-1]
        isb = dom.get(contig_of(orf),"B")!="A"
        lib_all[S]+=c; domsum[S]["A" if not isb else "B"]+=c
        if isb: lib_bact[S]+=c
        if kos:
            ko_all[kos[0]][S]+=c
            if isb: ko_bact[kos[0]][S]+=c
        for g in matches(kos,pref,pfam):
            tgt_all[g][S]+=c
            if isb: tgt_bact[g][S]+=c
    del counts,dom
    print(f"  {S}: lib_all={lib_all[S]:,} lib_bact={lib_bact[S]:,}")

def wmat(m,path):
    with open(path,"w") as o:
        o.write("KO\t"+"\t".join(SAMPLES)+"\n")
        for ko in sorted(m): o.write(ko+"\t"+"\t".join(str(m[ko].get(S,0)) for S in SAMPLES)+"\n")
wmat(ko_all,  f"{OUTDIR}/ko_counts_wholecommunity.tsv")
wmat(ko_bact, f"{OUTDIR}/ko_counts_bacterialonly.tsv")
with open(f"{OUTDIR}/domain_composition.tsv","w") as o:
    o.write("sample\tseason\tarch_reads\tbact_reads\tpct_archaea\n")
    for S in SAMPLES:
        a=domsum[S]["A"]; b=domsum[S]["B"]; t=a+b or 1
        o.write(f"{S}\t{SEASON[S]}\t{a}\t{b}\t{100*a/t:.2f}\n")
def rpm(x,l): return 1e6*x/l if l else 0.0
def rat(v): return ((v[2]+v[3])/2)/(((v[0]+v[1])/2) or 1e-9)
with open(f"{OUTDIR}/target_functions_persample.tsv","w") as o:
    o.write("function\t"+"\t".join(f"{S}_whole" for S in SAMPLES)+"\t"+
            "\t".join(f"{S}_bact" for S in SAMPLES)+"\tW:S_whole\tW:S_bact\n")
    for g in TARGETS:
        w=[rpm(tgt_all[g].get(S,0),lib_all[S]) for S in SAMPLES]
        bb=[rpm(tgt_bact[g].get(S,0),lib_bact[S]) for S in SAMPLES]
        o.write(g+"\t"+"\t".join(f"{x:.2f}" for x in w)+"\t"+
                "\t".join(f"{x:.2f}" for x in bb)+f"\t{rat(w):.2f}\t{rat(bb):.2f}\n")
print("[DONE]", OUTDIR)
print("KEY: in target_functions_persample.tsv compare W:S_whole vs W:S_bact.")
print("A function that stays winter-enriched (W:S_bact still >1) is real; one that")
print("drops toward 1.0 was a compositional artifact of the archaeal collapse.")
