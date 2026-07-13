suppressPackageStartupMessages({library(tidyverse)})
DIR <- "/data/habib/metagenome/biosurfactant/analysis_v2_community_abundance"
S16 <- "/data/habib/metagenome/biosurfactant/16S_data/level-6.csv"
clean_genus <- function(x) str_replace(x, "_[A-Z]$", "")

# ---- metagenome: MAG abundance per genus, by season (winterAbund/summerAbund = cols 14/15) ----
mag <- read_tsv(file.path(DIR,"meeting/3_outputs/mag_master_table.tsv"), show_col_types=FALSE) %>%
  filter(!is.na(genus), genus!="?", genus!="") %>%
  mutate(genus = clean_genus(genus)) %>%
  group_by(genus, phylum) %>%
  summarise(mg_summer = sum(summerAbund, na.rm=TRUE),
            mg_winter = sum(winterAbund, na.rm=TRUE), .groups="drop")

# ---- 16S: genus relative abundance per season ----
raw  <- readr::read_csv(S16, show_col_types=FALSE)
syn  <- c("Billgrantia"="Halomonas", "Franzmannia"="Halomonas")   # GTDB reclassification bridge
s16 <- raw %>% rename(sample=index) %>%
  mutate(across(-sample, as.numeric)) %>%
  pivot_longer(-sample, names_to="taxon", values_to="count") %>%
  filter(!is.na(count)) %>%
  mutate(genus  = str_match(taxon, "g__([^;]*)")[,2],
         season = ifelse(substr(sample,1,1)=="C","summer","winter")) %>%
  filter(!is.na(genus), genus!="") %>%
  mutate(genus = ifelse(genus %in% names(syn), syn[genus], genus)) %>%
  group_by(sample) %>% mutate(rel = 100*count/sum(count)) %>% ungroup() %>%
  group_by(season, genus) %>% summarise(v = mean(rel), .groups="drop") %>%
  pivot_wider(names_from=season, values_from=v, names_prefix="s16_", values_fill=0)

# ---- join on shared genera ----
tab <- inner_join(mag, s16, by="genus") %>%
  mutate(
    mg_pref  = ifelse(mg_winter  > mg_summer,  "winter", "summer"),
    s16_pref = ifelse(s16_winter > s16_summer, "winter", "summer"),
    concordant = ifelse(mg_pref == s16_pref, "Yes", "No"),
    display_genus = ifelse(genus=="Halomonas", "Halomonas (16S: Billgrantia/Franzmannia)", genus)
  ) %>%
  arrange(desc(pmax(mg_winter, mg_summer)))

message("shared genera: ", nrow(tab))
message("concordant: ", sum(tab$concordant=="Yes"), " / ", nrow(tab))

# ---- Spearman rho across all genus x season points ----
pts <- bind_rows(
  tab %>% transmute(mg=mg_summer, s16=s16_summer),
  tab %>% transmute(mg=mg_winter, s16=s16_winter))
rho <- cor(pts$s16, pts$mg, method="spearman")
message(sprintf("Spearman rho = %.3f", rho))

# ---- write raw TSV (all numbers) ----
tab %>% transmute(
  Genus = display_genus, Phylum = phylum,
  MAG_summer = round(mg_summer,3), MAG_winter = round(mg_winter,3),
  `16S_summer_pct` = round(s16_summer,3), `16S_winter_pct` = round(s16_winter,3),
  Metagenome_prefers = mg_pref, `16S_prefers` = s16_pref, Concordant = concordant
) %>% write_tsv(file.path(DIR,"Table2_mag_16S_concordance.tsv"))

# ---- write a formatted, paste-ready version ----
fmt <- tab %>% transmute(
  Genus = display_genus, Phylum = phylum,
  `Metagenome (S / W)` = sprintf("%.2f / %.2f", mg_summer, mg_winter),
  `16S %% (S / W)`     = sprintf("%.2f / %.2f", s16_summer, s16_winter),
  `Seasonal preference (MG | 16S)` = sprintf("%s | %s", mg_pref, s16_pref),
  Concordant = concordant)
writeLines(c(
  sprintf("Table 2. Concordance between metagenome MAG abundance and 16S rRNA relative abundance across recovered genera (n = %d). Spearman rho = %.2f. Concordant = %d/%d genera agree on seasonal preference.",
          nrow(tab), rho, sum(tab$concordant=="Yes"), nrow(tab)),
  "",
  paste(colnames(fmt), collapse=" | "),
  paste(rep("---", ncol(fmt)), collapse=" | "),
  apply(fmt, 1, function(r) paste(r, collapse=" | "))
), file.path(DIR,"Table2_formatted.md"))

cat("\n===== TABLE 2 PREVIEW =====\n")
print(as.data.frame(fmt), row.names=FALSE)
cat(sprintf("\nSpearman rho = %.3f | concordant %d/%d\n", rho, sum(tab$concordant=="Yes"), nrow(tab)))
