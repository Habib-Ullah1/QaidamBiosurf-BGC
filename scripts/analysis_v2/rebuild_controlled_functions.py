#!/usr/bin/env python3
import os, re, collections
BASE="/data/habib/metagenome/biosurfactant"; AN=f"{BASE}/analysis_v2_community_abundance"
OUT=f"{BASE}/integration/wp1_compositional"; os.makedirs(OUT, exist_ok=True)
S=["C3","C6","H2","H6"]; SUM=["C3","C6"]; WIN=["H2","H6"]
Q,PREF,OGS,KO=0,8,4,11
CARRIER={("H6","concoct_8"):"Moraxella_A",("H6","concoct_40"):"Loktanella",
         ("H2","concoct_19"):"Roseovarius",("H2","concoct_7"):"Halomonas_B"}
F={
 "phospholipaseA":  ({"K01058"}, r"^(pldA|plaA)$","PLA1",False,"Membrane","robust"),
 "lysophospholipase":({"K01048"},r"^(pldB|plbA|lypA)$","Hydrolase_4",False,"Membrane","supporting"),
 "FA_desaturase":   (set(), r"^(des|desA|desB|desC|desK)$","FA_desaturase",True,"Membrane","robust"),
 "cardiolipin_cfa": ({"K06131","K06132","K00574"}, r"^(cls|clsA|clsB|cfa)$",None,False,"Membrane","dropped"),
 "EPS_wza":         ({"K01991"}, r"^wza$",None,False,"EPS","robust"),
 "EPS_wzc":         ({"K16692"}, r"^wzc$",None,False,"EPS","robust"),
 "EPS_wzb":         ({"K01104"}, r"^wzb$",None,False,"EPS","robust"),
 "EPS_kps_group2":  ({"K01992","K09689","K16554"}, r"^(kpsE|kpsT|kpsD|kpsM|exoP)$",None,False,"EPS","supporting"),
 "EPS_wzxwzy_ctrl": (set(), r"^(wzx|wzy|wzm|wzt)$",None,False,"Control","control"),
 "ectoine":         ({"K06718","K06720","K00836","K10674"}, None, None, False,"Solute","robust"),
 "glycine_betaine": (set(), r"^(betA|betB|betT|betS|betI|gbsA|gbsB|opuD|opuAA|opuAB|opuAC|opuBA|opuCA|proV|proW|proX)$",None,False,"Solute","supporting"),
 "trehalose":       (set(), r"^(otsA|otsB|treY|treZ|treS|treP|treT|treR)$",None,False,"Solute","robust"),
 "K_uptake":        (set(), r"^(trkA|trkG|trkH|kdpA|kdpB|kdpC|kdpD|kup|ktrA|ktrB|ktrC|ktrD|ktrE)$",None,False,"SaltIn","robust"),
 "rhodopsin":       (set(), r".","Bac_rhodopsin",True,"SaltIn","robust"),
 "AAP_bacteriochl": (set(), r"^(pufL|pufM|pufA|pufB|pufC|puhA|bch[CDEFGHILMNOPXYZ]|acsF|crtB|crtI)$",None,False,"Energy","robust"),
 "N_fixation":      (set(), r"^(nifH|nifD|nifK|nifE|nifN|nifB)$",None,False,"Nitrogen","robust"),
 "denitrification": (set(), r"^(narG|narH|narI|napA|napB|nirK|nirS|norB|norC|nosZ)$",None,False,"Nitrogen","robust"),
 "sulfur_ox_sox":   (set(), r"^(soxA|soxB|soxC|soxX|soxY|soxZ)$",None,False,"Sulfur","robust"),
}
CARRIER_PANEL={"phospholipaseA":"Phospholipase A","lysophospholipase":"Lysophospholipase",
 "FA_desaturase":"Desaturase","EPS_wza":"EPS export","EPS_wzc":"EPS export",
 "EPS_wzb":"EPS export","EPS_kps_group2":"EPS export"}
def hit(fn,kos,pref,pfam):
    kset,rx,pf,req,_,_=F[fn]
    if req: return bool(pf and pf.lower() in (pfam or "").lower())
    if kset & set(kos): return True
    if rx and pref not in ("","-") and re.match(rx,pref,re.I): return True
    if pf and pf.lower() in (pfam or "").lower(): return True
    return False
c_of={}
for ln in open(f"{BASE}/integration/bgc_mag_mapping/contig_to_bin.tsv"):
    p=ln.rstrip("\n").split("\t")
    if len(p)>=3 and (p[0],p[2]) in CARRIER: c_of[(p[0],p[1])]=CARRIER[(p[0],p[2])]
tpm=collections.defaultdict(lambda:collections.defaultdict(float))
tpmb=collections.defaultdict(lambda:collections.defaultdict(float))
carr=collections.defaultdict(lambda:collections.defaultdict(float)); arch_frac={}
for s in S:
    egg=f"{BASE}/{s}/05_eggnog/{s}_eggnog.emapper.annotations"
    cdom=collections.defaultdict(lambda:[0,0]); ann={}
    for ln in open(egg):
        if ln.startswith("#"): continue
        p=ln.rstrip("\n").split("\t")
        if len(p)<=KO: continue
        orf=p[Q]; ct=orf.rsplit("_",1)[0]; ogs=p[OGS]
        cdom[ct][0 if ("|Archaea" in ogs and "|Bacteria" not in ogs) else 1]+=1
        ann[orf]=([k.replace("ko:","") for k in p[KO].split(",") if k not in ("","-")],p[PREF],p[-1])
    dom={c:("A" if v[0]>v[1] else "B") for c,v in cdom.items()}
    T={}
    for ln in open(f"{BASE}/integration/functional/tpm/{s}.tpm.tsv"):
        if ln.startswith("#"): continue
        p=ln.rstrip("\n").split("\t")
        if len(p)<2: continue
        try: T[p[0]]=float(p[-1])
        except ValueError: continue
    arch=sum(t for o,t in T.items() if dom.get(o.rsplit("_",1)[0],"B")=="A")
    bact=sum(t for o,t in T.items() if dom.get(o.rsplit("_",1)[0],"B")!="A") or 1
    arch_frac[s]=100*arch/(arch+bact)
    for orf,t in T.items():
        if t<=0: continue
        a=ann.get(orf)
        if not a: continue
        kos,pref,pfam=a; ct=orf.rsplit("_",1)[0]; isb=dom.get(ct,"B")!="A"; gen=c_of.get((s,ct))
        for fn in F:
            if hit(fn,kos,pref,pfam):
                tpm[fn][s]+=t
                if isb: tpmb[fn][s]+=t*1e6/bact
                if gen and fn in CARRIER_PANEL: carr[gen][CARRIER_PANEL[fn]]+=t
    print(f"  {s}: archaea {arch_frac[s]:.1f}%")
def ws(d):
    a=sum(d.get(x,0) for x in SUM); b=sum(d.get(x,0) for x in WIN); return (b/a if a else float("nan"))
print("\n%-20s %8s %8s %8s"%("function","WS_whole","WS_bact","tier"))
with open(f"{OUT}/function_table_controlled.tsv","w") as o:
    o.write("function\tclass\ttier\t"+"\t".join(S)+"\tWS_whole\tWS_bact\n")
    for fn in F:
        w=ws(tpm[fn]); b=ws(tpmb[fn])
        o.write(f"{fn}\t{F[fn][4]}\t{F[fn][5]}\t"+"\t".join(f"{tpm[fn].get(x,0):.2f}" for x in S)+f"\t{w:.3f}\t{b:.3f}\n")
        print("%-20s %8.3f %8.3f %8s"%(fn,w,b,F[fn][5]))
with open(f"{OUT}/archaeal_fraction.tsv","w") as o:
    o.write("sample\tseason\tpct_archaea\n")
    for s in S: o.write(f"{s}\t{'summer' if s in SUM else 'winter'}\t{arch_frac[s]:.2f}\n")
with open(f"{OUT}/carrier_table.tsv","w") as o:
    o.write("genus\tfunction\twinter_TPM\n")
    for g in carr:
        for fnc,v in carr[g].items(): o.write(f"{g}\t{fnc}\t{v:.2f}\n")
print("\n[DONE] clean single run ->", OUT)
