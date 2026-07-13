#!/usr/bin/env bash
set -euo pipefail
B=/data/habib/metagenome/biosurfactant
OUT=$B/integration/reground; mkdir -p "$OUT"
PERORF=$OUT/biosurfactant_specific_orfs.tsv
printf 'sample\torf\tcontig\tKO\tgene\tclass\n' > "$PERORF"

KOMAP=$(mktemp)
cat > "$KOMAP" <<'MAP'
K18100;rhlA;Rhamnolipid;A
K18101;rhlB_rhamnosyltransf;Rhamnolipid;A
K18102;rhlC;Rhamnolipid;A
K15654;srfAA;Surfactin;A
K15655;srfAB;Surfactin;A
K15656;srfAC;Surfactin;A
K15657;srfAD_TE;Surfactin;A
K15664;ppsA;Fengycin;A
K15665;ppsB;Fengycin;A
K15666;ppsC;Fengycin;A
K15667;ppsD;Fengycin;A
K15668;ppsE;Fengycin;A
K15659;ituA;Iturin;A
K15660;ituB;Iturin;A
K15661;ituC;Iturin;A
K01784;galE;UDP-glc-4-epim[mislabeled_rhlA];X
K00059;fabG;3-oxoacyl-ACP-red[generic];X
K11927;rhlE;RNA-helicase[false_friend];X
K03732;rhlB_helicase;RNA-helicase-RhlB[false_friend];X
MAP

EGG=$(find "$B" -name "*.emapper.annotations" | sort)
echo "eggNOG files:"; echo "$EGG" | sed 's/^/  /'; echo

awk -F'\t' -v komap="$KOMAP" -v perorf="$PERORF" '
BEGIN{ while((getline l < komap)>0){ n=split(l,a,";"); ko=a[1];
  gene[ko]=a[2]; cls[ko]=a[3]; tier[ko]=a[4]; want[ko]=1; order[++nk]=ko } }
FNR==1{ s=FILENAME; sub(/.*\//,"",s); split(s,sp,"_"); samp=sp[1]; kc=0;pc=0;dc=0 }
/^#query/{ for(i=1;i<=NF;i++){ h=$i; sub(/^#/,"",h);
  if(h=="KEGG_ko")kc=i; if(h=="Preferred_name")pc=i; if(h=="Description")dc=i } next }
/^#/{next}
{ if(kc==0) next; kk=$kc; if(kk=="-"||kk=="") next; gsub(/ko:/,"",kk);
  m=split(kk,kos,","); for(j=1;j<=m;j++){ ko=kos[j]; if(ko in want){
    cnt[ko SUBSEP samp]++; tot[ko]++;
    pn=(pc?$pc:"-"); if(pn=="")pn="-"; key=ko SUBSEP pn; pnc[key]++;
    if(!(key in pseen)){ pseen[key]=1; names[ko]=names[ko](names[ko]==""?"":" ")pn }
    if(!(ko in dseen)){ descex[ko]=(dc?$dc:"-"); dseen[ko]=1 }
    if(tier[ko]=="A"){ c=$1; sub(/_[0-9]+$/,"",c);
      print samp"\t"$1"\t"c"\t"ko"\t"gene[ko]"\t"cls[ko] >> perorf } } } }
END{
  print "=========================================================================";
  print "SECTION 1  VERIFICATION: what does eggNOG actually call each KO?";
  print "=========================================================================";
  printf "%-9s %-24s %-5s %-6s %s\n","KO","expected","tier","ORFs","eggNOG_Preferred_name[count]";
  for(i=1;i<=nk;i++){ ko=order[i]; ln="";
    if(ko in names){ x=split(names[ko],ar," "); for(t=1;t<=x;t++) ln=ln ar[t]"["pnc[ko SUBSEP ar[t]]"] " } else ln="(no hits)";
    printf "%-9s %-24s %-5s %-6d %s\n",ko,gene[ko],tier[ko],(tot[ko]?tot[ko]:0),ln }
  print "";
  print "  one example Description per KO with hits:";
  for(i=1;i<=nk;i++){ ko=order[i]; if(ko in descex) printf "    %-9s %s\n",ko,descex[ko] }
  print "";
  print "=========================================================================";
  print "SECTION 2  INFLATION WATCH (non-specific / false-friend KOs, tier X)";
  print "=========================================================================";
  for(i=1;i<=nk;i++){ ko=order[i]; if(tier[ko]=="X") printf "  %-9s %-42s ORFs=%d\n",ko,gene[ko],(tot[ko]?tot[ko]:0) }
  print "";
  print "=========================================================================";
  print "SECTION 3  VERIFIED-SPECIFIC markers, per sample (tier A)";
  print "=========================================================================";
  printf "  %-9s %-24s %-12s %5s %5s %5s %5s %7s\n","KO","gene","class","C3","C6","H2","H6","TOTAL";
  for(i=1;i<=nk;i++){ ko=order[i]; if(tier[ko]!="A") continue;
    printf "  %-9s %-24s %-12s %5d %5d %5d %5d %7d\n",ko,gene[ko],cls[ko],
      cnt[ko SUBSEP "C3"],cnt[ko SUBSEP "C6"],cnt[ko SUBSEP "H2"],cnt[ko SUBSEP "H6"],(tot[ko]?tot[ko]:0) }
}' $EGG
rm -f "$KOMAP"
echo; echo "per-ORF tier-A table: $PERORF  (rows: $(($(wc -l < "$PERORF")-1)))"
