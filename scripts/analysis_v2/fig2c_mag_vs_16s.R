source("fig_common.R")   # run from analysis_v2_community_abundance/
DIR  <- "/data/habib/metagenome/biosurfactant/analysis_v2_community_abundance"
S16  <- "/data/habib/metagenome/biosurfactant/16S_data/level-6.csv"

# --- metagenome MAG abundance per genus, by season ---
mag <- read_tsv(file.path(DIR,"meeting/3_outputs/mag_master_table.tsv"), show_col_types=FALSE)
message("mag_master_table cols: ", paste(colnames(mag), collapse=", "))
# expects genus + per-season abundance columns; adjust names below to match the printout
mg <- mag %>%
  mutate(genus = clean_genus(genus)) %>%
  group_by(genus) %>%
  summarise(mg_summer = sum(summerAbund, na.rm=TRUE),
            mg_winter = sum(winterAbund, na.rm=TRUE), .groups="drop")

# --- 16S genus relative abundance per season (from level-6.csv) ---
raw <- readr::read_csv(S16, show_col_types=FALSE)
long <- raw %>% rename(sample=index) %>%
  mutate(across(-sample, as.numeric)) %>%
  pivot_longer(-sample, names_to="taxon", values_to="count") %>%
  filter(!is.na(count)) %>%
  mutate(genus = str_match(taxon, "g__([^;]*)")[,2],
         season = ifelse(substr(sample,1,1)=="C","summer","winter")) %>%
  filter(!is.na(genus), genus!="") %>%
  group_by(sample) %>% mutate(rel = 100*count/sum(count)) %>% ungroup() %>%
  group_by(season, genus) %>% summarise(s16 = mean(rel), .groups="drop") %>%
  pivot_wider(names_from=season, values_from=s16, names_prefix="s16_", values_fill=0)

# --- GTDB genus synonym bridge (16S name -> MAG name); extend from the diagnostic ---
syn <- c("Billgrantia"="Halomonas", "Franzmannia"="Halomonas")
long <- long %>% mutate(genus = ifelse(genus %in% names(syn), syn[genus], genus))
# collapse any genera merged by the synonym map
long <- long %>% group_by(genus) %>%
  summarise(s16_summer=sum(s16_summer), s16_winter=sum(s16_winter), .groups="drop")

# --- join and reshape to one point per genus x season ---
comp <- inner_join(mg, long, by="genus")
message("shared genera: ", nrow(comp))
pts <- bind_rows(
  comp %>% transmute(genus, season="summer", mg=mg_summer, s16=s16_summer),
  comp %>% transmute(genus, season="winter", mg=mg_winter, s16=s16_winter)) %>%
  filter(mg>0 | s16>0)

rho <- cor(pts$s16, pts$mg, method="spearman")
message("Spearman rho (16S vs metagenome): ", round(rho,3))

p <- ggplot(pts, aes(s16, mg, color=season)) +
  geom_abline(slope=1, intercept=0, linetype=3, color="grey70") +
  geom_point(size=2.6, alpha=.9) +
  ggrepel::geom_text_repel(aes(label=genus), size=2.6, max.overlaps=15, show.legend=FALSE) +
  scale_color_manual(values=SEASON, name=NULL) +
  scale_x_continuous(trans="log1p") + scale_y_continuous(trans="log1p") +
  labs(x="16S relative abundance (%)", y="Metagenome MAG abundance (%)",
       subtitle=sprintf("Recovered genera: 16S vs metagenome (Spearman \u03c1 = %.2f)", rho)) +
  theme_pub() + theme(legend.position="top")

ggsave(file.path(DIR,"Fig2C_mag_vs_16S.pdf"), p, width=6, height=5.5)
ggsave(file.path(DIR,"Fig2C_mag_vs_16S.png"), p, width=6, height=5.5, dpi=300)
message("wrote Fig2C")
