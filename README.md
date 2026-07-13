# QaidamBiosurf-BGC

Genome-resolved metagenomic analysis of microbial adaptation in the Qaidam Basin
(hyperarid, hypersaline sediments), companion to a 16S rRNA amplicon study.
Part of a project investigating the role of surface-active/amphiphilic lipids in
microbial adaptation to extreme environments.

## Study design
- 4 sediment metagenomes: **C3, C6 = summer**; **H2, H6 = winter**
- Companion 16S dataset: 13 sediment samples (summer/winter)
- Assembly: MEGAHIT (`k141_` contigs); genome binning: metaWRAP/CONCOCT (44 MAGs)

## Repository layout
- `scripts/` — original BGC-focused pipeline (assembly, binning, antiSMASH, BiG-SCAPE,
  GTDB-Tk, CheckM2, CoverM, NRPS domain analysis).
- `analysis_v2_community_abundance/` — whole-community re-analysis and verified results.
  Contains the corrected abundance, functional profiling, taxonomy, CAZy/EPS, and the
  biosurfactant verdict. See its `README.md` for a file-by-file description.

## Key findings
1. **Seasonal osmoadaptation switch.** Summer is dominated by salt-in haloarchaea and
   *Salinibacter* (acidic proteomes, rhodopsin photoheterotrophy); winter by salt-out
   bacteria accumulating ectoine/glycine-betaine/trehalose. Demonstrated by an
   annotation-free proteome-acidity test and whole-community functional abundance.
2. **Surface-active amphiphile response (project thesis link).** Winter-enriched and
   multi-modal: membrane-lipid remodeling (phospholipase/lysophospholipase ~3–4×,
   fatty-acid desaturase ~4.5× [PFAM-verified]) matching lysophospholipids detected by
   metabolomics, plus secreted EPS/capsule export (Wza/Wzc/Kps, 4–15× winter). Carried
   by the winter salt-out cohort (*Moraxella_A*, *Loktanella*, *Roseovarius*,
   *Halomonas_B*) and, at community scale, by unrecovered lineages.
3. **Biosurfactant (lipopeptide) hypothesis — tested and refined.** Dedicated lipopeptide
   biosynthesis is essentially absent (1 of 1,468 BGCs with a lipopeptide starter-C
   domain), confirmed four independent ways (KO gene names, antiSMASH region domains,
   full C-domain inventory, and C-domain phylogeny against biosurfactant references).
   The BGC pool is dominated by terpene/RiPP-like classes; novelty is at the sequence
   level, not in compound class. The community's amphiphile adaptation is membrane
   remodeling + EPS, not secreted biosurfactant production.

## Methodological notes (important for interpretation)
- **Abundance is coverage/TPM-based over the whole assembly (binned + unbinned).**
  Earlier MAG-only relative abundance (CoverM genome mode) normalized within the
  ~6% binned fraction and is superseded; functional abundance uses per-ORF TPM
  (featureCounts on full-assembly BAMs). TPM is compositional — "winter-enriched"
  means a higher community *share*, not absolute cell counts.
- **Taxonomy uses sample-aware GTDB-Tk r214** (per-sample `11_gtdbtk/`); the deduplicated
  summary drops sample prefixes and must not be used for contig/gene joins.
- **Annotation caveats:** eggNOG KO/CAZy calls are conservative (lower bounds; consistent
  across samples, so seasonal ratios hold). Individual gene-name matches were verified by
  KO+PFAM before being attributed to named organisms.

## Data availability
Raw reads, BAMs, and per-ORF matrices are large and kept on compute storage, not in this
repository; summarized/analyzable tables are provided here.
