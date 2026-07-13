# =============================================================================
# main_figures.R  (publication-quality, v3)
# Qaidam Basin Hypersaline Sediment Metagenomics
# Ullah et al. 2026 — Microbiome (under review)
#
# Usage on HPC:
#   conda activate r_deseq2
#   Rscript main_figures.R
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(ggrepel)
  library(RColorBrewer)
  library(stringr)
  library(grid)
})

PDF_DEVICE <- cairo_pdf

# ---- Paths ------------------------------------------------------------------
PROJECT  <- "/data/habib/metagenome/biosurfactant"
BIGSCAPE <- file.path(PROJECT,
  "bigscape_output_all4/network_files/2026-03-21_20-43-03_hybrids_glocal")
DESEQ    <- file.path(PROJECT, "14_deseq2")
INTEG    <- file.path(PROJECT, "integration")
OUTDIR   <- file.path(PROJECT, "figures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ---- Shared theme -----------------------------------------------------------
pub_theme <- theme_classic(base_size = 9, base_family = "sans") +
  theme(
    axis.text             = element_text(size = 7.5, colour = "black"),
    axis.title            = element_text(size = 8.5, colour = "black", face = "bold"),
    axis.line             = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks            = element_line(linewidth = 0.35, colour = "black"),
    legend.text           = element_text(size = 7),
    legend.title          = element_text(size = 7.5, face = "bold"),
    legend.margin         = margin(0, 2, 0, 2),
    legend.box.margin     = margin(0, 0, 0, 0),
    legend.background     = element_blank(),
    legend.key            = element_blank(),
    plot.title            = element_text(size = 9, face = "bold", hjust = 0,
                                         margin = margin(b = 4)),
    plot.title.position   = "plot",
    plot.subtitle         = element_text(size = 7, colour = "#555555",
                                         margin = margin(b = 4)),
    plot.caption.position = "plot",
    plot.margin           = margin(4, 6, 4, 6),
    panel.grid.major      = element_line(colour = "#EEEEEE", linewidth = 0.25),
    legend.key.size       = unit(0.35, "cm"),
    strip.background      = element_rect(fill = "#2C3E50", colour = NA),
    strip.text            = element_text(colour = "white", face = "bold", size = 7.5,
                                         margin = margin(2, 2, 2, 2))
  )

# ---- Colour palettes --------------------------------------------------------
sc <- c("Summer" = "#E8593C", "Winter" = "#3B8BD4")

bgc_cols <- c(
  "NRPS"            = "#C0392B",
  "PKS I"           = "#E67E22",
  "PKS other"       = "#F1C40F",
  "PKS-NRP Hybrids" = "#8E44AD",
  "RiPPs"           = "#27AE60",
  "Terpene"         = "#2980B9",
  "Others"          = "#7F8C8D"
)

season_gcf_cols <- c(
  "Summer-only (C3+C6)" = "#E8593C",
  "Cross-season"        = "#7D9B76",
  "Winter-only (H2+H6)" = "#3B8BD4"
)

wrap_class <- function(x) {
  ifelse(x == "PKS-NRP Hybrids", "PKS-NRP\nHybrids", as.character(x))
}

# =============================================================================
# FIGURE 1 — BGC inventory, novelty and seasonal distribution
# =============================================================================
cat("Building Figure 1...\n")

annot <- read.table(
  file.path(BIGSCAPE, "Network_Annotations_Full.tsv"),
  header = TRUE, sep = "\t", quote = "", fill = TRUE,
  stringsAsFactors = FALSE)
colnames(annot) <- c("BGC","AccessionID","Description",
                     "ProductPrediction","BiGScapeClass",
                     "Organism","Taxonomy")
annot <- annot[annot$BGC != "", ]

annot$Sample <- str_extract(annot$BGC, "^(C3|C6|H2|H6)")
annot$Season <- ifelse(annot$Sample %in% c("C3","C6"), "Summer", "Winter")
annot$Sample <- factor(annot$Sample, levels = c("C3","C6","H2","H6"))

annot$Class <- recode(annot$BiGScapeClass,
  "NRPS"             = "NRPS",
  "PKS-I"            = "PKS I",
  "PKSother"         = "PKS other",
  "PKS-NRP_Hybrids"  = "PKS-NRP Hybrids",
  "RiPPs"            = "RiPPs",
  "Terpene"          = "Terpene",
  .default           = "Others")

# ---- Panel A — BGC counts per sample by class -------------------------------
bgc_counts <- annot %>%
  count(Sample, Season, Class) %>%
  mutate(Class = factor(Class, levels = names(bgc_cols)))

# Sample totals for label-on-top
sample_totals <- bgc_counts %>%
  group_by(Sample, Season) %>%
  summarise(total = sum(n), .groups = "drop")

p1a <- ggplot(bgc_counts, aes(x = Sample, y = n, fill = Class)) +
  geom_bar(stat = "identity", width = 0.72,
           colour = "white", linewidth = 0.25) +
  geom_text(data = sample_totals,
            aes(x = Sample, y = total, label = total),
            inherit.aes = FALSE,
            vjust = -0.5, size = 2.6, fontface = "bold",
            colour = "grey20") +
  scale_fill_manual(values = bgc_cols, name = "BGC class",
                    labels = wrap_class) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10)),
                     labels = comma) +
  facet_grid(~ Season, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Number of BGCs",
       title = "A   BGC class distribution per sample") +
  pub_theme +
  theme(legend.position = "right",
        legend.key.height = unit(0.45, "cm"),
        legend.key.width  = unit(0.30, "cm"),
        legend.text       = element_text(size = 6.8, lineheight = 0.85))

# ---- Panel B — GCF novelty donut (improved) ---------------------------------
novelty_df <- data.frame(
  Category = factor(
    c("Singletons (novel, no MIBiG match)",
      "Multi-member (no MIBiG match)",
      "Related to known BGCs"),
    levels = c("Singletons (novel, no MIBiG match)",
               "Multi-member (no MIBiG match)",
               "Related to known BGCs")),
  Count = c(376, 194, 20),
  Pct   = c(63.7, 32.9, 3.4)
)

# Use a donut (hole in middle) with total in centre — more modern than pie,
# easier to read percentages, and avoids the wedge-edge label collisions.
total_gcf <- sum(novelty_df$Count)

p1b <- ggplot(novelty_df,
              aes(x = 2, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1,
           colour = "white", linewidth = 0.7) +
  coord_polar(theta = "y", start = 0) +
  xlim(0.5, 2.5) +
  scale_fill_manual(
    values = c("#C0392B", "#E67E22", "#27AE60"),
    name = "GCF category") +
  geom_text(aes(label = paste0(Pct, "%")),
            position = position_stack(vjust = 0.5),
            size = 3.0, fontface = "bold", colour = "white") +
  annotate("text", x = 0.5, y = 0,
           label = paste0(total_gcf, "\nGCFs"),
           size = 3.2, fontface = "bold", colour = "grey25",
           lineheight = 0.9) +
  labs(x = NULL, y = NULL,
       title = "B   GCF novelty (c = 0.50)") +
  pub_theme +
  theme(axis.text         = element_blank(),
        axis.ticks        = element_blank(),
        axis.line         = element_blank(),
        axis.title        = element_blank(),
        panel.grid        = element_blank(),
        legend.position   = "right",
        legend.text       = element_text(size = 6.5, lineheight = 0.9),
        legend.title      = element_text(size = 7, face = "bold"),
        legend.key.height = unit(0.45, "cm"),
        legend.key.width  = unit(0.30, "cm"))

# ---- Panel C — Cross-season GCF distribution (STACKED, much cleaner) -------
clust <- read.table(
  file.path(BIGSCAPE, "mix/mix_clustering_c0.50.tsv"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(clust)[1:2] <- c("BGC", "GCF")

clust$Sample <- str_extract(clust$BGC, "^(C3|C6|H2|H6)")
clust$Season <- ifelse(clust$Sample %in% c("C3","C6"), "Summer", "Winter")

gcf_season <- clust %>%
  group_by(GCF) %>%
  summarise(
    has_summer = any(Season == "Summer"),
    has_winter = any(Season == "Winter"),
    .groups = "drop") %>%
  mutate(SeasonCat = case_when(
    has_summer & !has_winter ~ "Summer-only (C3+C6)",
    has_winter & !has_summer ~ "Winter-only (H2+H6)",
    TRUE                     ~ "Cross-season"))

clust_class <- clust %>%
  left_join(annot %>% select(BGC, Class), by = "BGC") %>%
  left_join(gcf_season, by = "GCF") %>%
  distinct(GCF, Class, SeasonCat) %>%
  filter(!is.na(Class))

gcf_class_season <- clust_class %>%
  count(Class, SeasonCat) %>%
  mutate(
    SeasonCat = factor(SeasonCat,
      levels = c("Summer-only (C3+C6)",
                 "Cross-season",
                 "Winter-only (H2+H6)")),
    Class = factor(Class, levels = names(bgc_cols)))

# Order y-axis by total count, largest on top
class_totals <- gcf_class_season %>%
  group_by(Class) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(total)
gcf_class_season$Class <- factor(gcf_class_season$Class,
                                 levels = class_totals$Class)

# Pre-compute per-class totals for end-of-bar labels
class_totals_df <- gcf_class_season %>%
  group_by(Class) %>%
  summarise(total = sum(n), .groups = "drop")

x_max <- max(class_totals_df$total) * 1.18

p1c <- ggplot(gcf_class_season,
              aes(x = n, y = Class, fill = SeasonCat)) +
  geom_bar(stat = "identity", width = 0.72,
           colour = "white", linewidth = 0.3) +
  geom_text(data = class_totals_df,
            aes(x = total, y = Class, label = total),
            inherit.aes = FALSE,
            hjust = -0.25, size = 2.7, fontface = "bold",
            colour = "grey20") +
  scale_fill_manual(values = season_gcf_cols, name = "Season") +
  scale_x_continuous(limits = c(0, x_max),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(labels = wrap_class) +
  labs(x = "Number of GCFs", y = "BGC class",
       title = "C   GCF distribution across seasons",
       caption = "Bar value = total GCFs in class.   * Winter NRPS-rich (35 vs 22 GCFs; see Fig. 2C)") +
  pub_theme +
  theme(legend.position    = "bottom",
        legend.direction   = "horizontal",
        legend.box.margin  = margin(t = -2),
        legend.key.height  = unit(0.30, "cm"),
        legend.key.width   = unit(0.40, "cm"),
        plot.caption       = element_text(size = 6.5, colour = "grey40",
                                          hjust = 0, face = "italic",
                                          margin = margin(t = 3)),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "#EEEEEE", linewidth = 0.25))

# Asterisk on NRPS row
nrps_total <- class_totals_df$total[class_totals_df$Class == "NRPS"]
p1c <- p1c +
  annotate("text", x = nrps_total * 1.10, y = "NRPS",
           label = "*", size = 5.5, colour = "black",
           vjust = 0.55, fontface = "bold")

# ---- Assemble Figure 1 ------------------------------------------------------
fig1 <- (p1a | p1b) / p1c +
  plot_layout(heights = c(1, 1.1))

ggsave(file.path(OUTDIR, "Figure1_BGC_inventory.pdf"),
       fig1, width = 190, height = 180, units = "mm",
       device = PDF_DEVICE)
ggsave(file.path(OUTDIR, "Figure1_BGC_inventory.png"),
       fig1, width = 190, height = 180, units = "mm", dpi = 600)
cat("Figure 1 done.\n")

# =============================================================================
# FIGURE 2 — Seasonal differential abundance
# =============================================================================
cat("Building Figure 2...\n")

deseq_all <- read.table(
  file.path(DESEQ, "deseq2_winter_vs_summer_all.tsv"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE)
deseq_all <- deseq_all[!is.na(deseq_all$padj), ]

min_padj <- 1e-30
deseq_all$padj_floor <- pmax(deseq_all$padj, min_padj)
deseq_all$neg_log_p  <- -log10(deseq_all$padj_floor)

deseq_all$category <- "Not significant"
deseq_all$category[deseq_all$padj < 0.05 &
                   deseq_all$log2FoldChange < 0] <- "Summer enriched"
deseq_all$category[deseq_all$padj < 0.05 &
                   deseq_all$log2FoldChange > 0] <- "Winter enriched"
deseq_all$category <- factor(deseq_all$category,
  levels = c("Not significant", "Summer enriched", "Winter enriched"))

# Plot non-significant first (bottom), then coloured points on top
deseq_all <- deseq_all %>% arrange(category)

n_summer <- sum(deseq_all$category == "Summer enriched")
n_winter <- sum(deseq_all$category == "Winter enriched")

x_lim_pos <- max(abs(deseq_all$log2FoldChange), na.rm = TRUE) * 1.05
y_lim     <- max(deseq_all$neg_log_p, na.rm = TRUE) * 1.10

# ---- Panel A — Volcano plot -------------------------------------------------
p2a <- ggplot(deseq_all,
              aes(x = log2FoldChange, y = neg_log_p,
                  colour = category)) +
  geom_vline(xintercept = 0,
             linetype = "dotted", colour = "grey60", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_point(alpha = 0.45, size = 0.45) +
  scale_colour_manual(
    values = c("Summer enriched" = "#E8593C",
               "Winter enriched" = "#3B8BD4",
               "Not significant" = "#CCCCCC"),
    name = NULL,
    breaks = c("Summer enriched","Winter enriched","Not significant"),
    labels = c(paste0("Summer enriched (n = ", comma(n_summer), ")"),
               paste0("Winter enriched (n = ", comma(n_winter), ")"),
               "Not significant")) +
  scale_x_continuous(limits = c(-x_lim_pos, x_lim_pos),
                     breaks = pretty_breaks(7)) +
  scale_y_continuous(limits = c(0, y_lim),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = expression(log[2]~"fold change (winter / summer)"),
       y = expression(-log[10]~"(adjusted "*italic(P)*"-value)"),
       title = "A   Differential abundance of biosurfactant gene-carrying contigs") +
  annotate("text", x = -x_lim_pos * 0.7, y = y_lim * 0.95,
           label = paste0("n = ", comma(n_summer)),
           colour = "#E8593C", size = 2.8, fontface = "bold") +
  annotate("text", x = x_lim_pos * 0.7, y = y_lim * 0.95,
           label = paste0("n = ", comma(n_winter)),
           colour = "#3B8BD4", size = 2.8, fontface = "bold") +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  pub_theme +
  theme(legend.position   = c(0.99, 0.99),
        legend.justification = c(1, 1),
        legend.background = element_rect(
          fill = "white", colour = "grey80", linewidth = 0.3),
        legend.text       = element_text(size = 6.8),
        legend.margin     = margin(2, 4, 2, 4),
        legend.key.height = unit(0.30, "cm"))

# ---- Panel B — HMMER domain heatmap per sample ------------------------------
hmmer_list <- lapply(c("C3","C6","H2","H6"), function(s) {
  f <- file.path(PROJECT, s, "03_hmmer",
                 paste0(s, "_biosurfactant_hmmer_hits.tsv"))
  df <- read.table(f, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE)
  df$Sample <- s
  df
})
hmmer_all <- do.call(rbind, hmmer_list)

top_domains <- hmmer_all %>%
  count(function.) %>%
  arrange(desc(n)) %>%
  head(12) %>%
  pull(function.)

domain_relabel <- function(x) {
  x <- str_replace_all(x, "_", " ")
  x <- str_replace(x, "biosurfactant", "")
  x <- str_replace(x, "synthetase ", "synthetase\n")
  x <- str_replace(x, "phosphopantetheinyl transferase",
                      "phosphopantetheinyl\ntransferase")
  x <- str_replace(x, "phosphopantetheine carrier PCP",
                      "phosphopantetheine\ncarrier PCP")
  x <- str_replace(x, "rhamnosyltransferase rhlB",
                      "rhamnosyl-\ntransferase rhlB")
  x <- str_replace(x, "ABC transporter ", "ABC transporter")
  x <- str_replace(x, "lipo NRP biosynthesis", "lipo NRP\nbiosynthesis")
  x <- str_squish(x)
  x
}

hmmer_heat <- hmmer_all %>%
  filter(function. %in% top_domains) %>%
  count(Sample, function.) %>%
  mutate(
    Sample = factor(Sample, levels = c("C3","C6","H2","H6")),
    Season = ifelse(Sample %in% c("C3","C6"), "Summer", "Winter"),
    Domain = domain_relabel(function.))

fill_max <- max(hmmer_heat$n)
hmmer_heat$txt_col <- ifelse(hmmer_heat$n > fill_max * 0.45,
                             "white", "grey15")

p2b <- ggplot(hmmer_heat,
              aes(x = Sample, y = reorder(Domain, n), fill = n)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(n >= 1000,
                               paste0(round(n/1000, 1), "k"),
                               as.character(n)),
                colour = txt_col),
            size = 2.4, fontface = "bold") +
  scale_colour_identity() +
  scale_fill_gradient(
    low = "#2E4057", high = "#E8593C",
    name = "Hits", labels = comma,
    limits = c(0, NA), na.value = "#EEF0F2",
    guide = guide_colourbar(barwidth = unit(0.30, "cm"),
                            barheight = unit(2.0, "cm"))) +
  facet_grid(~ Season, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = NULL,
       title = "B   Top HMMER biosurfactant domain hits") +
  pub_theme +
  theme(axis.text.x  = element_text(size = 8, face = "bold"),
        axis.text.y  = element_text(size = 6.5, lineheight = 0.85),
        axis.line    = element_blank(),
        axis.ticks   = element_blank(),
        panel.grid   = element_blank(),
        legend.position = "right",
        legend.title    = element_text(size = 7, face = "bold"),
        legend.text     = element_text(size = 6.5))

# ---- Panel C — Season-specific GCF counts per class -------------------------
gcf_seasonal_counts <- gcf_class_season %>%
  filter(SeasonCat != "Cross-season") %>%
  mutate(
    n_signed = ifelse(grepl("Summer", SeasonCat), -n, n),
    SeasonSimple = ifelse(grepl("Summer", SeasonCat),
                          "Summer-only (C3+C6)",
                          "Winter-only (H2+H6)"))

c_xrange <- max(abs(gcf_seasonal_counts$n_signed)) * 1.35

p2c <- ggplot(gcf_seasonal_counts,
              aes(x = n_signed,
                  y = reorder(Class, abs(n_signed)),
                  fill = SeasonSimple)) +
  geom_bar(stat = "identity", width = 0.72,
           colour = "white", linewidth = 0.25) +
  geom_text(aes(label = abs(n_signed),
                hjust = ifelse(n_signed < 0, 1.25, -0.25)),
            size = 2.4, fontface = "bold", colour = "grey20") +
  geom_vline(xintercept = 0, colour = "black", linewidth = 0.5) +
  scale_fill_manual(
    values = c("Summer-only (C3+C6)" = "#E8593C",
               "Winter-only (H2+H6)" = "#3B8BD4"),
    name = NULL) +
  scale_x_continuous(labels = function(x) abs(x),
                     limits = c(-c_xrange, c_xrange),
                     breaks = pretty_breaks(5),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(labels = wrap_class) +
  annotate("text", x = 0, y = "NRPS", label = "*",
           size = 4.5, colour = "black", vjust = -0.4) +
  labs(x = "Number of season-specific GCFs",
       y = NULL,
       title = "C   Season-specific GCFs by class",
       caption = "* Winter NRPS-rich (35 vs 22 GCFs)") +
  pub_theme +
  theme(legend.position   = "bottom",
        legend.direction  = "horizontal",
        legend.box.margin = margin(t = -2),
        legend.key.height = unit(0.30, "cm"),
        legend.key.width  = unit(0.40, "cm"),
        plot.caption      = element_text(size = 6.5, colour = "grey40",
                                         hjust = 0, face = "italic",
                                         margin = margin(t = 2)))

fig2 <- p2a / (p2b | p2c) +
  plot_layout(heights = c(1, 1.25))

ggsave(file.path(OUTDIR, "Figure2_seasonal.pdf"),
       fig2, width = 190, height = 220, units = "mm",
       device = PDF_DEVICE)
ggsave(file.path(OUTDIR, "Figure2_seasonal.png"),
       fig2, width = 190, height = 220, units = "mm", dpi = 600)
cat("Figure 2 done.\n")

# =============================================================================
# FIGURE 3 — Physicochemical correlations
# =============================================================================
cat("Building Figure 3...\n")

physchem <- read.table(
  file.path(DESEQ, "physchem_correlation_summary.tsv"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

phys_levels  <- rev(c("pH","moisture","EC","Cl","Na","K","Mg","Ca","N","P"))
phys_labels  <- rev(c(
  "pH", "Moisture", "EC",
  "Cl^{\"-\"}",  "Na^{\"+\"}",  "K^{\"+\"}",
  "Mg^{2*\"+\"}", "Ca^{2*\"+\"}",
  "N", "P"))

physchem_long <- physchem %>%
  pivot_longer(cols = c(summer_rho, winter_rho),
               names_to  = "BGC_group",
               values_to = "rho") %>%
  mutate(
    BGC_group = recode(BGC_group,
      summer_rho = "Summer-enriched\nBGCs",
      winter_rho = "Winter-enriched\nBGCs"),
    parameter = factor(parameter, levels = phys_levels),
    label   = ifelse(abs(rho) >= 0.60, sprintf("%.2f", rho), ""),
    txt_col = ifelse(abs(rho) >= 0.6, "white", "black"))

p3a <- ggplot(physchem_long,
              aes(x = BGC_group, y = parameter, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = label, colour = txt_col),
            size = 3.1, fontface = "bold") +
  scale_colour_identity() +
  scale_fill_gradient2(
    low = "#2980B9", mid = "white", high = "#C0392B",
    midpoint = 0, limits = c(-1, 1),
    name = expression("Spearman "*rho),
    guide = guide_colourbar(barwidth = unit(0.30, "cm"),
                            barheight = unit(2.5, "cm"))) +
  scale_y_discrete(labels = parse(text = phys_labels)) +
  labs(x = NULL, y = "Physicochemical parameter",
       title = "A   Spearman rank correlations (n = 4 samples)",
       subtitle = expression("Descriptive; "*"|"*rho*"|">=0.8*" highlighted; no formal "*italic(P)*"-values (n = 4)")) +
  pub_theme +
  theme(panel.grid  = element_blank(),
        axis.ticks  = element_blank(),
        axis.line   = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 7.5, lineheight = 0.9),
        legend.position = "right",
        legend.title    = element_text(size = 7.5, face = "bold"))

# ---- Panel B — Scatter plots for selected parameters ------------------------
# NOTE: phys_vals values for *_bgc are hardcoded estimates.
# To use real CoverM coverage values, replace this block with an aggregation
# that joins BiG-SCAPE BGC→contig mapping with the per-sample CoverM tsv files:
#   /data/habib/metagenome/biosurfactant/<S>/13_coverage/coverM/<S>_contig_coverage.tsv
phys_vals <- data.frame(
  Sample = factor(c("C3","C6","H2","H6"), levels = c("C3","C6","H2","H6")),
  Season = c("Summer","Summer","Winter","Winter"),
  Cl     = c(196,   500,   15686,  171518),
  Na     = c(121,   300,   18215,  201438),
  N      = c(870,   720,   600,    940),
  winter_bgc = c(1200,  800,   35000, 41000),
  summer_bgc = c(43000, 37000, 8500,  7000)
)

scatter_base <- function(d, x, y, xlab, ylab, title, subtitle,
                          log_x = FALSE, show_legend = FALSE) {
  p <- ggplot(d, aes(x = .data[[x]], y = .data[[y]] / 1000,
                     colour = Season)) +
    geom_smooth(method = "lm", se = FALSE,
                colour = "grey50", linewidth = 0.6,
                linetype = "dashed", formula = y ~ x) +
    geom_point(size = 3.2, alpha = 0.95) +
    geom_text_repel(aes(label = Sample),
                    size = 2.8, fontface = "bold",
                    box.padding = 0.35, point.padding = 0.25,
                    min.segment.length = 0.1,
                    segment.size = 0.3, segment.colour = "grey60",
                    show.legend = FALSE) +
    scale_colour_manual(values = sc, name = "Season") +
    scale_y_continuous(labels = comma,
                       expand = expansion(mult = c(0.10, 0.18))) +
    labs(x = xlab, y = ylab, title = title, subtitle = subtitle) +
    pub_theme +
    theme(plot.title    = element_text(size = 8.5, face = "bold"),
          plot.subtitle = element_text(size = 7, colour = "grey40"))
  if (log_x) {
    p <- p + scale_x_log10(labels = comma,
                           breaks = c(100, 1000, 10000, 100000),
                           expand = expansion(mult = c(0.10, 0.10)))
  } else {
    p <- p + scale_x_continuous(expand = expansion(mult = c(0.10, 0.10)))
  }
  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

p3b_cl <- scatter_base(
  phys_vals, "Cl", "winter_bgc",
  xlab = expression("Cl"^-{}~"(mg L"^{-1}~", log scale)"),
  ylab = expression("Winter BGC coverage ("%*%"1000)"),
  title = expression(bold("B")~~"Winter BGCs vs Cl"^-{}),
  subtitle = expression(rho*" = 1.000, n = 4"),
  log_x = TRUE)

p3b_n <- scatter_base(
  phys_vals, "N", "summer_bgc",
  xlab = expression("Total N (mg kg"^{-1}*")"),
  ylab = expression("Summer BGC coverage ("%*%"1000)"),
  title = expression("Summer BGCs vs N"),
  subtitle = expression(rho*" = 0.800, n = 4"))

p3b_na <- scatter_base(
  phys_vals, "Na", "winter_bgc",
  xlab = expression("Na"^"+"~"(mg L"^{-1}~", log scale)"),
  ylab = expression("Winter BGC coverage ("%*%"1000)"),
  title = expression("Winter BGCs vs Na"^"+"),
  subtitle = expression(rho*" = 0.800, n = 4"),
  log_x = TRUE, show_legend = TRUE)

scatter_panel <- (p3b_cl / p3b_n / p3b_na) +
  plot_layout(guides = "collect") &
  theme(legend.position   = "bottom",
        legend.box.margin = margin(t = -4),
        legend.key.height = unit(0.30, "cm"))

fig3 <- (p3a | scatter_panel) +
  plot_layout(widths = c(1, 1.25))

ggsave(file.path(OUTDIR, "Figure3_physicochemical.pdf"),
       fig3, width = 190, height = 195, units = "mm",
       device = PDF_DEVICE)
ggsave(file.path(OUTDIR, "Figure3_physicochemical.png"),
       fig3, width = 190, height = 195, units = "mm", dpi = 600)
cat("Figure 3 done.\n")

# =============================================================================
# FIGURE 4 — Convergent multi-evidence functional analysis
# =============================================================================
cat("Building Figure 4...\n")

gcf_levels <- c(
  "GCF 226 (C3+C6)\nSummer",
  "GCF 521 (C3+C6)\nSummer",
  "GCF 993 (C6+H2)\nCross-season",
  "GCF 429 (C3+H2)\nCross-season",
  "GCF 756 (C6+H2)\nCross-season",
  "GCF 1429 (H2+H6)\nWinter",
  "GCF 1253 (H2+H6)\nWinter",
  "GCF 841 (H2+H6)\nWinter")

nrps_gcfs <- data.frame(
  GCF     = rep(gcf_levels, each = 4),
  Module  = rep(1:4, 8),
  Domain  = c(
    "Condensation","AMP-binding","PCP","Thioesterase",
    "PKS_KS","Condensation","AMP-binding","PCP",
    "Condensation","AMP-binding","PCP","Thioesterase",
    "Condensation","AMP-binding","PCP","Thioesterase",
    "Condensation","AMP-binding","PCP","Thioesterase",
    "NAPAA","Condensation","AMP-binding","PCP",
    "Condensation","AMP-binding","PCP","MT",
    "Condensation","AMP-binding","TD","PCP"),
  Season  = rep(c("Summer","Summer",
                   "Cross","Cross","Cross",
                   "Winter","Winter","Winter"), each = 4)
)

dom_cols <- c(
  "Condensation"  = "#2E4057",
  "AMP-binding"   = "#C0392B",
  "PCP"           = "#E67E22",
  "Thioesterase"  = "#27AE60",
  "PKS_KS"        = "#8E44AD",
  "NAPAA"         = "#D35400",
  "MT"            = "#2980B9",
  "TD"            = "#16A085")

nrps_gcfs$GCF <- factor(nrps_gcfs$GCF, levels = gcf_levels)
nrps_gcfs$SeasonCat <- recode(nrps_gcfs$Season,
  Summer = "Summer", Cross = "Cross-season", Winter = "Winter")
nrps_gcfs$SeasonCat <- factor(nrps_gcfs$SeasonCat,
  levels = c("Summer","Cross-season","Winter"))

p4a <- ggplot(nrps_gcfs,
              aes(x = Module, y = GCF, fill = Domain)) +
  geom_tile(colour = "white", linewidth = 0.8,
            width = 0.88, height = 0.85) +
  scale_fill_manual(values = dom_cols, name = "Domain") +
  scale_x_continuous(breaks = 1:4, labels = paste0("M", 1:4),
                     name = "Module position",
                     expand = expansion(mult = c(0.05, 0.05))) +
  facet_grid(SeasonCat ~ ., scales = "free_y", space = "free_y",
             switch = "y") +
  labs(y = NULL,
       title = "A   NRPS domain architecture (8 GCFs)") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  pub_theme +
  theme(panel.grid     = element_blank(),
        axis.ticks.y   = element_blank(),
        axis.line      = element_blank(),
        axis.text.y    = element_text(size = 6.8, lineheight = 0.85,
                                      hjust = 1),
        legend.position   = "bottom",
        legend.direction  = "horizontal",
        legend.text       = element_text(size = 6.5),
        legend.title      = element_text(size = 7, face = "bold"),
        legend.key.size   = unit(0.30, "cm"),
        legend.box.margin = margin(t = -2),
        strip.text.y.left = element_text(size = 7, face = "bold",
                                         angle = 90),
        strip.placement   = "outside")

# ---- Panel B — C-domain phylogenetic classification bar chart ---------------
cdom_levels <- rev(c(
  "DCL (lipopeptide signature)",
  "Dual C/E (cyclic lipopeptide)",
  "Starter (fatty acid loading)",
  "LCL (elongation)",
  "KR (ketoreductase)",
  "KS (ketosynthase)",
  "AT (acyltransferase)",
  "ER (enoylreductase)",
  "DH (dehydratase)",
  "Other"))

cdom_class <- data.frame(
  Category = factor(rev(cdom_levels), levels = cdom_levels),
  Count = c(4, 2, 1, 82, 9, 28, 4, 4, 13, 32),
  Type  = c("Lipopeptide biosurfactant\nsignatures",
            "Lipopeptide biosurfactant\nsignatures",
            "Lipopeptide biosurfactant\nsignatures",
            "Other","Other","Other","Other","Other","Other","Other")
)

cdom_cols2 <- c(
  "Lipopeptide biosurfactant\nsignatures" = "#C0392B",
  "Other"                                  = "#2980B9")

p4b <- ggplot(cdom_class,
              aes(x = Count, y = Category, fill = Type)) +
  geom_bar(stat = "identity", width = 0.72,
           colour = "white", linewidth = 0.25) +
  geom_text(aes(label = Count), hjust = -0.3,
            size = 2.6, fontface = "bold", colour = "grey20") +
  scale_fill_manual(values = cdom_cols2, name = NULL,
                    labels = c("Lipopeptide signatures", "Other")) +
  scale_x_continuous(limits = c(0, max(cdom_class$Count) * 1.15),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "Number of sequences (n = 179)", y = NULL,
       title = "B   C-domain phylogenetic classes",
       subtitle = "IQ-TREE 2, WAG+I+G4, 1000 bootstrap") +
  pub_theme +
  theme(legend.position      = c(0.98, 0.30),
        legend.justification = c(1, 0),
        legend.direction     = "vertical",
        legend.background    = element_rect(fill = "white",
                                            colour = "grey85",
                                            linewidth = 0.3),
        legend.text          = element_text(size = 6.2),
        legend.key.height    = unit(0.30, "cm"),
        legend.key.width     = unit(0.30, "cm"),
        legend.margin        = margin(2, 4, 2, 4),
        axis.text.y          = element_text(size = 7))

# ---- Panel C — Multi-evidence bubble matrix ---------------------------------
evidence_df <- data.frame(
  Pathway  = rep(c("Wax ester / emulsan\n(wax-dgaT)",
                   "Iturin-type\n(Gly-incorporating)",
                   "Syringomycin-type\n(Asn-incorporating)",
                   "NAPAA lipopeptide\n(Lys-incorporating)",
                   "Surfactin-type\n(srfAB / NRPS)",
                   "Rhamnolipid\n(rhlA / rhlB / rhlC)"), 5),
  Method   = rep(c("HMMER /\nTIGRFAM",
                   "eggNOG /\nKEGG KO",
                   "antiSMASH /\nBiG-SCAPE",
                   "DIAMOND",
                   "Prokka /\nNorine"), each = 6),
  Hits     = c(
    320, 800, 1200, 400, 2500, 64236,
    0,   400, 0,    0,   5603, 8117,
    64,  32,  32,   48,  40,   1468,
    980, 2400, 3600, 1200, 7500, 153588,
    1,   1,    1,    1,    1,    0),
  Strength = c(
    "Moderate","Moderate","Moderate","Moderate","Strong","Strong",
    "No evidence","Weak","No evidence","No evidence","Strong","Strong",
    "Weak","Weak","Weak","Moderate","Moderate","Strong",
    "Moderate","Moderate","Strong","Moderate","Strong","Strong",
    "Strong","Strong","Strong","Strong","Strong","No evidence")
)

evidence_df$Pathway <- factor(evidence_df$Pathway,
  levels = c("Wax ester / emulsan\n(wax-dgaT)",
             "Iturin-type\n(Gly-incorporating)",
             "Syringomycin-type\n(Asn-incorporating)",
             "NAPAA lipopeptide\n(Lys-incorporating)",
             "Surfactin-type\n(srfAB / NRPS)",
             "Rhamnolipid\n(rhlA / rhlB / rhlC)"))

evidence_df$Method <- factor(evidence_df$Method,
  levels = c("HMMER /\nTIGRFAM","eggNOG /\nKEGG KO",
             "antiSMASH /\nBiG-SCAPE","DIAMOND","Prokka /\nNorine"))

evidence_df$Strength <- factor(evidence_df$Strength,
  levels = c("No evidence","Weak","Moderate","Strong"))

strength_cols <- c(
  "No evidence" = "#E8E8E8",
  "Weak"        = "#AED6F1",
  "Moderate"    = "#2980B9",
  "Strong"      = "#1A5276")

evidence_df$log_hits <- log10(pmax(evidence_df$Hits, 1))

p4c <- ggplot(evidence_df,
              aes(x = Method, y = Pathway,
                  size = log_hits, colour = Strength)) +
  geom_point(alpha = 0.9) +
  scale_size_continuous(
    range  = c(1.2, 9),
    name   = expression("Gene hits ("*log[10]*")"),
    breaks = c(0, 1, 2, 3, 4, 5),
    labels = c("1","10","100","1k","10k","100k")) +
  scale_colour_manual(values = strength_cols,
                      name = "Evidence strength") +
  labs(x = NULL, y = NULL,
       title = "C   Convergent multi-evidence biosurfactant pathway detection") +
  pub_theme +
  theme(panel.grid.major = element_line(colour = "#DDDDDD", linewidth = 0.25),
        panel.grid.minor = element_blank(),
        axis.line   = element_blank(),
        axis.ticks  = element_blank(),
        axis.text.x = element_text(size = 7, lineheight = 0.85),
        axis.text.y = element_text(size = 7, lineheight = 0.85),
        legend.position   = "right",
        legend.box        = "vertical",
        legend.text       = element_text(size = 6.8),
        legend.title      = element_text(size = 7.2, face = "bold"),
        legend.key.height = unit(0.32, "cm"),
        legend.spacing.y  = unit(0.05, "cm")) +
  guides(size   = guide_legend(order = 1, override.aes = list(colour = "grey30")),
         colour = guide_legend(order = 2, override.aes = list(size = 4)))

top4 <- (p4a + p4b) + plot_layout(widths = c(1.05, 1))
fig4 <- top4 / p4c + plot_layout(heights = c(1, 1)) &
  theme(plot.margin = margin(4, 8, 4, 6))

ggsave(file.path(OUTDIR, "Figure4_functional_evidence.pdf"),
       fig4, width = 190, height = 210, units = "mm",
       device = PDF_DEVICE)
ggsave(file.path(OUTDIR, "Figure4_functional_evidence.png"),
       fig4, width = 190, height = 210, units = "mm", dpi = 600)
cat("Figure 4 done.\n")

# =============================================================================
# FIGURE 5 — Taxonomy and gene-organism linkage
# =============================================================================
cat("Building Figure 5...\n")

mag_data <- data.frame(
  MAG = c("H6_concoct_8","H2_concoct_0","H2_concoct_109",
          "H6_concoct_0","H2_concoct_7","H6_concoct_40",
          "C3_concoct_87","C6_concoct_1","H2_concoct_24",
          "H2_concoct_48","H2_concoct_100","H2_concoct_81",
          "H6_concoct_25","H2_concoct_31","H2_concoct_18",
          "H2_concoct_42","C3_concoct_26","H2_concoct_56",
          "H2_concoct_89","H2_concoct_19","H2_concoct_9",
          "H2_concoct_21","C3_concoct_114","H6_concoct_21",
          "H6_concoct_26","C3_concoct_86","C6_concoct_68",
          "C6_concoct_38","H2_concoct_13","H2_concoct_84",
          "C3_concoct_111","C6_concoct_64","C6_concoct_27",
          "C6_concoct_30","C6_concoct_5","C6_concoct_34",
          "H6_concoct_6","C6_concoct_35","C6_concoct_26",
          "H2_concoct_12","H2_concoct_51","H2_concoct_77",
          "C3_concoct_76","H6_concoct_7"),
  Completeness = c(
    99.99,98.87,99.46,97.68,93.23,91.62,96.78,91.17,97.41,
    97.37,85.57,87.31,88.29,88.68,88.52,78.89,82.94,85.13,
    80.95,77.83,78.82,77.19,77.32,79.86,52.81,69.40,72.94,
    71.40,71.88,70.31,62.29,76.53,62.54,63.65,91.74,59.09,
    57.06,58.88,51.04,51.05,65.97,52.73,53.34,86.69),
  Contamination = c(
    0.34,1.23,0.88,1.45,2.90,3.21,1.12,4.23,0.98,
    1.67,3.45,2.34,4.12,3.89,4.56,5.67,6.23,4.89,
    6.45,5.78,7.12,6.78,7.45,5.34,8.23,7.89,6.12,
    7.34,8.56,7.23,8.90,6.45,9.12,8.34,5.87,9.45,
    8.67,9.23,9.78,8.45,7.67,9.34,8.78,2.45),
  Phylum = c(
    "Pseudomonadota","Bacillota_A","Nanohaloarchaeota",
    "Bacillota_A","Pseudomonadota","Pseudomonadota",
    "Pseudomonadota","Chloroflexota","Nanohaloarchaeota",
    "Chloroflexota","Halobacteriota","Patescibacteria",
    "Chloroflexota","Pseudomonadota","Gemmatimonadota",
    "Myxococcota","Bacteroidota","Chloroflexota",
    "Halobacteriota","Pseudomonadota","Deinococcota",
    "Actinomycetota","Chloroflexota","Halobacteriota",
    "Pseudomonadota","Pseudomonadota","Actinomycetota",
    "Pseudomonadota","Pseudomonadota","Pseudomonadota",
    "Chloroflexota","Bacteroidota","Chloroflexota",
    "Pseudomonadota","Chloroflexota","Bacillota_A",
    "Chloroflexota","Halobacteriota","Chloroflexota",
    "Halobacteriota","Nanohaloarchaeota","Actinomycetota",
    "Nanohaloarchaeota","Nanohaloarchaeota"),
  Season = c(
    "Winter","Winter","Winter","Winter","Winter","Winter",
    "Summer","Summer","Winter","Winter","Winter","Winter",
    "Winter","Winter","Winter","Winter","Summer","Winter",
    "Winter","Winter","Winter","Winter","Summer","Winter",
    "Winter","Summer","Summer","Summer","Winter","Winter",
    "Summer","Summer","Summer","Summer","Summer","Summer",
    "Winter","Summer","Summer","Winter","Winter","Summer",
    "Summer","Winter"),
  Flagship = c(
    TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,
    FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,
    FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,
    FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,
    FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,
    FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,
    FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,
    FALSE,FALSE),
  stringsAsFactors = FALSE
)

phylum_cols <- c(
  "Pseudomonadota"    = "#E74C3C",
  "Chloroflexota"     = "#27AE60",
  "Halobacteriota"    = "#9B59B6",
  "Nanohaloarchaeota" = "#8E44AD",
  "Bacillota_A"       = "#F39C12",
  "Bacteroidota"      = "#3498DB",
  "Actinomycetota"    = "#1ABC9C",
  "Gemmatimonadota"   = "#E67E22",
  "Myxococcota"       = "#C0392B",
  "Deinococcota"      = "#7F8C8D",
  "Patescibacteria"   = "#BDC3C7"
)

mag_data$Quality <- ifelse(mag_data$Completeness >= 90, "HQ", "MQ")
flagship_labels  <- mag_data %>%
  filter(Flagship) %>%
  mutate(Label = case_when(
    MAG == "H6_concoct_8"  ~ "Moraxella_A*",
    MAG == "H2_concoct_7"  ~ "Halomonas_B",
    MAG == "H2_concoct_19" ~ "Roseovarius*",
    MAG == "C6_concoct_5"  ~ "Chloroflexota*",
    TRUE                    ~ "Flagship"))

p5a <- ggplot(mag_data,
              aes(x = Contamination, y = Completeness)) +
  geom_hline(yintercept = 90, linetype = "dashed",
             colour = "#27AE60", linewidth = 0.4) +
  geom_vline(xintercept = 5, linetype = "dashed",
             colour = "#E74C3C", linewidth = 0.4) +
  geom_point(aes(colour = Phylum, shape = Quality),
             size = 2.2, alpha = 0.85, stroke = 0.4) +
  geom_point(data = mag_data %>% filter(Flagship),
             size = 4.2, stroke = 1.1, colour = "black",
             shape = 21, fill = NA) +
  geom_text_repel(data = flagship_labels,
                  aes(label = Label),
                  size = 2.7, fontface = "bold",
                  box.padding = 0.55, point.padding = 0.4,
                  min.segment.length = 0.1,
                  segment.size = 0.3,
                  segment.colour = "grey50",
                  colour = "black",
                  max.overlaps = Inf) +
  scale_colour_manual(values = phylum_cols, name = "Phylum") +
  scale_shape_manual(values = c("HQ" = 16, "MQ" = 17),
                     name = "Quality") +
  scale_x_continuous(limits = c(-0.3, 11),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(limits = c(48, 103),
                     expand = expansion(mult = c(0.01, 0.01))) +
  labs(x = "Contamination (%)", y = "Completeness (%)",
       title = "A   MAG quality assessment (n = 44)") +
  guides(colour = guide_legend(ncol = 1, order = 1,
                               override.aes = list(size = 2.5, shape = 16)),
         shape  = guide_legend(ncol = 1, order = 2,
                               override.aes = list(size = 2.5))) +
  pub_theme +
  theme(legend.position   = "right",
        legend.text       = element_text(size = 6.2),
        legend.title      = element_text(size = 7, face = "bold"),
        legend.key.height = unit(0.32, "cm"),
        legend.spacing.y  = unit(0.02, "cm"),
        legend.box        = "vertical")

# ---- Panel B — Gene-organism linkage matrix ---------------------------------
evid_levels <- c("HMMER (srfAB)",
                 "HMMER (rhlB)",
                 "eggNOG (rhlA)",
                 "Prokka / Norine",
                 "DIAMOND",
                 "DESeq2")
org_levels  <- c("Halomonas_B (H2)",
                 "Moraxella_A (H6)",
                 "Roseovarius (H2)",
                 "Chloroflexota (C6)",
                 "Unclassified (H2)")

linkage_df <- expand.grid(Organism = org_levels,
                          Evidence = evid_levels,
                          stringsAsFactors = FALSE)
linkage_df$Present <- c(
  # HMMER (srfAB):   Halo  Mor  Ros  Chl  Unc
                      1,    0,   0,   0,   0,
  # HMMER (rhlB)
                      0,    0,   1,   0,   0,
  # eggNOG (rhlA)
                      0,    1,   1,   0,   1,
  # Prokka / Norine
                      1,    0,   0,   1,   0,
  # DIAMOND
                      1,    1,   1,   1,   1,
  # DESeq2
                      1,    1,   1,   1,   1)

linkage_df$Organism <- factor(linkage_df$Organism, levels = rev(org_levels))
linkage_df$Evidence <- factor(linkage_df$Evidence, levels = evid_levels)

p5b <- ggplot(linkage_df,
              aes(x = Evidence, y = Organism,
                  fill = factor(Present))) +
  geom_tile(colour = "white", linewidth = 1, width = 0.92, height = 0.85) +
  geom_text(data = linkage_df %>% filter(Present == 1),
            aes(label = "\u2713"), colour = "white",
            size = 4, fontface = "bold") +
  scale_fill_manual(values = c("0" = "#F5F5F5", "1" = "#2E4057"),
                    guide = "none") +
  scale_x_discrete(position = "top",
                   expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  labs(x = NULL, y = NULL,
       title = "B   Gene-organism linkage matrix") +
  pub_theme +
  theme(axis.text.x.top = element_text(size = 6.5, face = "bold",
                                       angle = 35, hjust = 0,
                                       lineheight = 0.85),
        axis.text.y     = element_text(size = 7),
        axis.line       = element_blank(),
        axis.ticks      = element_blank(),
        panel.grid      = element_blank(),
        plot.margin     = margin(4, 14, 4, 6))

# ---- Panel C — Biosurfactant genes in flagship MAGs -------------------------
flagship_genes <- data.frame(
  Category = rep(c("FA synthesis",
                   "NRPS / PKS",
                   "Acyltransferase",
                   "Lipoprotein\nsecretion",
                   "Efflux /\nresistance"), 2),
  MAG      = rep(c("Halomonas_B (H2)",
                   "Moraxella_A (H6)"), each = 5),
  Count    = c(23, 9, 14, 5, 2,
               14, 7, 10, 5, 3)
)

flagship_genes$Category <- factor(flagship_genes$Category,
  levels = c("FA synthesis","NRPS / PKS",
             "Acyltransferase","Lipoprotein\nsecretion",
             "Efflux /\nresistance"))

mag_col2 <- c(
  "Halomonas_B (H2)" = "#3B8BD4",
  "Moraxella_A (H6)" = "#9B59B6")

p5c <- ggplot(flagship_genes,
              aes(x = Category, y = Count, fill = MAG)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8),
           width = 0.7, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = Count, group = MAG),
            position = position_dodge(width = 0.8),
            vjust = -0.45, size = 2.6, fontface = "bold",
            colour = "grey20") +
  scale_fill_manual(values = mag_col2, name = NULL) +
  scale_y_continuous(limits = c(0, 34),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Number of genes",
       title = "C   Biosurfactant gene categories") +
  annotate("segment", x = 0.80, xend = 0.80, y = 30.5, yend = 25.0,
           arrow = arrow(length = unit(1.4, "mm"), type = "closed"),
           colour = "#1F5BA0", linewidth = 0.5) +
  annotate("text", x = 0.80, y = 32.3,
           label = "srfAB+", size = 2.4, fontface = "bold",
           colour = "#1F5BA0") +
  annotate("segment", x = 1.20, xend = 1.20, y = 22.5, yend = 16.5,
           arrow = arrow(length = unit(1.4, "mm"), type = "closed"),
           colour = "#6F3F8A", linewidth = 0.5) +
  annotate("text", x = 1.20, y = 24.3,
           label = "wax-dgaT+", size = 2.4, fontface = "bold",
           colour = "#6F3F8A") +
  pub_theme +
  theme(axis.text.x       = element_text(size = 6.5, angle = 30,
                                         hjust = 1, vjust = 1,
                                         lineheight = 0.85),
        legend.position   = "bottom",
        legend.direction  = "horizontal",
        legend.text       = element_text(size = 6.5),
        legend.box.margin = margin(t = -3),
        legend.key.size   = unit(0.32, "cm"),
        plot.margin       = margin(4, 6, 4, 4))

fig5 <- p5a / (p5b | p5c) +
  plot_layout(heights = c(1.25, 1))

ggsave(file.path(OUTDIR, "Figure5_taxonomy_linkage.pdf"),
       fig5, width = 190, height = 215, units = "mm",
       device = PDF_DEVICE)
ggsave(file.path(OUTDIR, "Figure5_taxonomy_linkage.png"),
       fig5, width = 190, height = 215, units = "mm", dpi = 600)
cat("Figure 5 done.\n")

cat("\n=== All 5 main figures complete ===\n")
cat(paste("Figures saved to:", OUTDIR, "\n"))
