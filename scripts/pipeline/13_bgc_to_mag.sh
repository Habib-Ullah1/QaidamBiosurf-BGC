#!/usr/bin/env bash
set -euo pipefail
B=/data/habib/metagenome/biosurfactant
OUT=$B/integration/bgc_mag_mapping; mkdir -p "$OUT"

echo "STEP 1: Build contig-to-bin lookup from 44 MAG FASTAs..."
> "$OUT/contig_to_bin.tsv"
for fa in "$B"/*/11_gtdbtk/filtered_bins/*.fa; do
  sample=$(echo "$fa" | awk -F/ '{for(i=1;i<=NF;i++) if($i~/^(C3|C6|H2|H6)$/) print $i}')
  bin=$(basename "$fa" .fa)
  grep "^>" "$fa" | awk -v s="$sample" -v b="$bin" '{c=$1; sub(/^>/,"",c); print s"\t"c"\t"b}' >> "$OUT/contig_to_bin.tsv"
done
echo "  contigs in bins: $(wc -l < "$OUT/contig_to_bin.tsv")"

echo "STEP 2: Build GTDB-Tk taxonomy lookup..."
> "$OUT/gtdbtk_combined.tsv"
for tsv in "$B"/*/11_gtdbtk/output/classify/gtdbtk.*.summary.tsv "$B"/*/11_gtdbtk/output/gtdbtk.*.summary.tsv; do
  [ -f "$tsv" ] || continue
  tail -n+2 "$tsv" >> "$OUT/gtdbtk_combined.tsv"
done
awk -F'\t' '!seen[$1]++' "$OUT/gtdbtk_combined.tsv" > "$OUT/gtdbtk_dedup.tsv"
echo "  unique MAG taxonomies: $(wc -l < "$OUT/gtdbtk_dedup.tsv")"

echo "STEP 3: Extract BGC regions with product class from antiSMASH..."
> "$OUT/bgc_regions.tsv"
for S in C3 C6 H2 H6; do
  find "$B/$S/02_antismash" -name "*.region*.gbk" 2>/dev/null | while read gbk; do
    region=$(basename "$gbk" .gbk)
    contig=$(echo "$region" | sed 's/\.region[0-9]*//')
    product=$(grep '/product=' "$gbk" 2>/dev/null | head -1 | sed 's/.*product="\([^"]*\)".*/\1/' || echo "unknown")
    echo -e "${S}\t${contig}\t${region}\t${product}"
  done
done >> "$OUT/bgc_regions.tsv"
echo "  total BGC regions: $(wc -l < "$OUT/bgc_regions.tsv")"

echo "STEP 4: Join BGCs to bins and taxonomy..."
awk -F'\t' '
NR==FNR && FILENAME ~ /contig_to_bin/ { bin[$1 SUBSEP $2] = $3; next }
NR==FNR && FILENAME ~ /gtdbtk_dedup/  { tax[$1] = $2; next }
{
  key = $1 SUBSEP $2
  b = (key in bin) ? bin[key] : "unbinned"
  t = (b in tax) ? tax[b] : "unclassified"
  print $0"\t"b"\t"t
}
' "$OUT/contig_to_bin.tsv" "$OUT/gtdbtk_dedup.tsv" "$OUT/bgc_regions.tsv" > "$OUT/bgc_mag_taxonomy.tsv"

echo "STEP 5: Summary tables..."
echo
echo "============================================================"
echo "A. BGCs per MAG (top 20 MAGs by BGC count)"
echo "============================================================"
awk -F'\t' '$5!="unbinned"' "$OUT/bgc_mag_taxonomy.tsv" \
  | awk -F'\t' '{print $1"_"$5"\t"$6}' \
  | sort | uniq -c | sort -rn | head -20 \
  | awk '{printf "  %4d BGCs  %-25s %s\n",$1,$2,$3}'

echo
echo "============================================================"
echo "B. BGC class distribution: binned vs unbinned"
echo "============================================================"
awk -F'\t' '{
  if($5=="unbinned") ub[$4]++; else bn[$4]++;
  tot[$4]++
} END {
  printf "  %-25s %6s %6s %6s %6s\n","product","binned","unbind","total","%binned";
  PROCINFO["sorted_in"]="@val_num_desc";
  for(p in tot) printf "  %-25s %6d %6d %6d %5.1f%%\n",p,(bn[p]?bn[p]:0),(ub[p]?ub[p]:0),tot[p],(bn[p]?bn[p]:0)*100/tot[p]
}' "$OUT/bgc_mag_taxonomy.tsv" | sort -t'%' -k1 -rn | head -20

echo
echo "============================================================"
echo "C. BGC class per sample (all 1468)"
echo "============================================================"
awk -F'\t' '{s[$1][$4]++; cls[$4]++} END {
  printf "  %-20s %5s %5s %5s %5s %7s\n","product","C3","C6","H2","H6","TOTAL";
  PROCINFO["sorted_in"]="@val_num_desc";
  for(c in cls) printf "  %-20s %5d %5d %5d %5d %7d\n",c,(s["C3"][c]?s["C3"][c]:0),(s["C6"][c]?s["C6"][c]:0),(s["H2"][c]?s["H2"][c]:0),(s["H6"][c]?s["H6"][c]:0),cls[c]
}' "$OUT/bgc_mag_taxonomy.tsv"

echo
echo "============================================================"
echo "D. Phylum-level BGC attribution (binned only)"
echo "============================================================"
awk -F'\t' '$5!="unbinned" { split($6,a,";"); ph=a[2]; sub(/^p__/,"",ph); if(ph=="") ph="unclassified"; cnt[ph]++ } END {
  PROCINFO["sorted_in"]="@val_num_desc";
  for(p in cnt) printf "  %-25s %d\n",p,cnt[p]
}' "$OUT/bgc_mag_taxonomy.tsv"

echo
echo "============================================================"
echo "E. Per-MAG BGC class matrix (binned BGCs)"
echo "============================================================"
awk -F'\t' '$5!="unbinned" {
  mag=$1"_"$5; cls[$4]++; m[mag][$4]++; mtot[mag]++
} END {
  n=asorti(cls,clist);
  printf "  %-25s","MAG";
  for(i=1;i<=n;i++) printf " %6s",substr(clist[i],1,6);
  printf " %6s\n","TOTAL";
  PROCINFO["sorted_in"]="@val_num_desc";
  for(mag in mtot){
    printf "  %-25s",mag;
    for(i=1;i<=n;i++) printf " %6d",(m[mag][clist[i]]?m[mag][clist[i]]:0);
    printf " %6d\n",mtot[mag]
  }
}' "$OUT/bgc_mag_taxonomy.tsv"

echo
echo "Files written to: $OUT/"
echo "  bgc_mag_taxonomy.tsv  (full join: sample, contig, region, product, bin, taxonomy)"
