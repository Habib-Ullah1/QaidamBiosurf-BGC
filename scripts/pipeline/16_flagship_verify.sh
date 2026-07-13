#!/usr/bin/env bash
set -uo pipefail
B=/data/habib/metagenome/biosurfactant
OUT=$B/integration/flagship_verify; mkdir -p "$OUT"
PFAM=$B/databases/pfam/Pfam-A.hmm

# --- locate hmmscan (try base, then antismash/bigscape envs) ---
if ! command -v hmmscan >/dev/null 2>&1; then
  source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null
  for e in antismash antismash_env bigscape_env; do
    conda activate "$e" 2>/dev/null && command -v hmmscan >/dev/null 2>&1 && break
  done
fi
echo "hmmscan: $(command -v hmmscan || echo NOT-FOUND)"

HALO_FAA=$(ls $B/H2/prokka_concoct7/*.faa 2>/dev/null | head -1)
MORX_FAA=$(ls $B/H6/prokka_concoct8/*.faa 2>/dev/null | head -1)
echo "Halomonas faa: $HALO_FAA"
echo "Moraxella faa: $MORX_FAA"

FLAG=$OUT/flagship_proteins.faa; : > "$FLAG"

extract_tag () { awk -v tag="$2" -v lab="$3" '
  /^>/{p=0; if(index($0,tag)){p=1; print ">"lab" "substr($0,2)} next} p{print}' "$1"; }
extract_prod () { awk -v lab="$3" -v re="$2" '
  /^>/{p=0; if($0 ~ re){p=1; print ">"lab" "substr($0,2)} next} p{print}' "$1"; }

extract_tag  "$HALO_FAA" "PMNBPKPI_03120" "Halomonas_srfAB" >> "$FLAG"
grep -q "Halomonas_srfAB" "$FLAG" || extract_prod "$HALO_FAA" "[Ss]urfactin synthase" "Halomonas_srfAB" >> "$FLAG"
extract_tag  "$MORX_FAA" "BMDNJEMN_00474" "Moraxella_WSD" >> "$FLAG"
grep -q "Moraxella_WSD" "$FLAG" || extract_prod "$MORX_FAA" "O-acyltransferase WSD" "Moraxella_WSD" >> "$FLAG"

echo; echo "=== Extracted flagship proteins ==="; grep "^>" "$FLAG"
echo "--- lengths ---"
awk '/^>/{if(n)print name": "n" aa"; name=substr($1,2); n=0; next}{n+=length($0)}END{if(n)print name": "n" aa"}' "$FLAG"

echo; echo "=== Running hmmscan vs Pfam-A (few proteins, ~1-2 min) ==="
hmmscan --cpu 4 --domtblout "$OUT/flagship_pfam.domtbl" -E 1e-5 "$PFAM" "$FLAG" > "$OUT/flagship_pfam.log" 2>&1
echo "exit: $?  domtbl: $OUT/flagship_pfam.domtbl"

echo; echo "=== DOMAIN ARCHITECTURE (ordered N->C by envelope position) ==="
awk '!/^#/{printf "%-22s %-22s acc=%-10s i-Eval=%s  env=%s-%s\n",$4,$1,$2,$13,$20,$21}' "$OUT/flagship_pfam.domtbl" \
  | sort -k1,1 -k6,6
echo
echo "=== KEY DOMAIN CHECK ==="
echo "-- Halomonas_srfAB should carry NRPS module (Condensation/AMP-binding/PP-binding):"
awk '!/^#/ && $4=="Halomonas_srfAB"{print "   "$1" ("$2")  i-E="$13}' "$OUT/flagship_pfam.domtbl" \
  | grep -iE "Condensation|AMP-binding|PP-binding|Thioester|NRPS|ACP" || echo "   (no NRPS domains found!)"
echo "-- Moraxella_WSD should carry wax-ester synthase / DGAT (PF03007):"
awk '!/^#/ && $4=="Moraxella_WSD"{print "   "$1" ("$2")  i-E="$13}' "$OUT/flagship_pfam.domtbl" \
  | grep -iE "PF03007|wax|DGAT|acyltransf|Diacylglyc|WES" || echo "   (no wax-ester domain found!)"

echo; echo "=== NAPAA BGC (Moraxella_A, k141_497210) — antiSMASH's own domain calls ==="
GBK=$(find "$B/H6/02_antismash" -name "k141_497210.region*.gbk" 2>/dev/null | head -1)
echo "gbk: $GBK"
if [ -n "$GBK" ]; then
  echo "--- region product ---"; grep -m1 "/product=" "$GBK"
  echo "--- aSDomains / NRPS-PKS domains present ---"
  grep -oE "aSDomain=\"[^\"]+\"" "$GBK" | sort | uniq -c
  echo "--- functional domain mentions (AMP-binding/PP-binding/Condensation/acyltransferase/GNAT) ---"
  grep -oiE "AMP-binding|PP-binding|Condensation|Acetyltransf|GNAT|TIGR02353|N-acyl" "$GBK" | sort | uniq -c
fi
