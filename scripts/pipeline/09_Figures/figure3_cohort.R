#!/usr/bin/env Rscript
# Figure 3 — Seasonal cohort structure and distribution of biosynthetic potential
suppressWarnings(suppressMessages({
  need <- c("ggplot2","ggbeeswarm","ggrepel","patchwork","dplyr","scales","grid")
  miss <- need[!need %in% rownames(installed.packages())]
  if (length(miss)) install.packages(miss, repos="https://cloud.r-project.org")
  library(ggplot2); library(ggbeeswarm); library(ggrepel)
  library(patchwork); library(dplyr); library(scales); library(grid)
}))

B <- "/data/habib/metagenome/biosurfactant"
dat <- read.delim(file.path(B,"integration/bgc_mag_mapping/genome_economics.tsv"),
                  stringsAsFactors=FALSE)

dat$cohort <- ifelse(dat$season %in% c("winter","winter-lean"),"Winter",
              ifelse(dat$season %in% c("summer","summer-lean"),"Summer","Year-round"))
dat$cohort <- factor(dat$cohort, levels=c("Winter","Year-round","Summer"))
coh_pal <- c(Winter="#2166AC", `Year-round`="#9A9A9A", Summer="#D55E00")

# floor tiny abundances so the log axis behaves (true values <0.002 shown as <0.002)
dat$wAb <- pmax(dat$winterAb, 0.002); dat$sAb <- pmax(dat$summerAb, 0.002)
lab3 <- dat %>%
  filter(genus %in% c("Ralstonia","Halomonas_B","Persicimonas","Gillisia","Moraxella_A","Loktanella")) %>%
  group_by(genus) %>% slice_max(nBGC, n=1) %>% ungroup()

th <- theme_classic(base_size=8, base_family="sans") +
  theme(axis.text=element_text(colour="#444444", size=7.5),
        axis.title=element_text(colour="#111111", size=8.5),
        axis.line=element_line(colour="#444444", linewidth=0.4),
        axis.ticks=element_line(colour="#444444", linewidth=0.4),
        plot.tag=element_text(face="bold", size=12),
        legend.position="none", plot.margin=margin(5,7,3,5))

set.seed(7)
# ---- Panel a: seasonal abundance switch ----
pa <- ggplot(dat, aes(sAb, wAb)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="#888888", linewidth=0.45) +
  annotate("text", x=0.0035, y=2.6, label="winter-enriched", hjust=0, size=2.5,
           colour="#2166AC", fontface="italic") +
  annotate("text", x=0.5, y=0.0045, label="summer-enriched", hjust=1, size=2.5,
           colour="#D55E00", fontface="italic") +
  geom_point(aes(fill=cohort, size=nBGC), shape=21, stroke=0.4, colour="white", alpha=0.9) +
  geom_text_repel(data=lab3, aes(sAb, wAb, label=genus), fontface="italic", size=2.8,
                  segment.colour="#888888", segment.size=0.3, min.segment.length=0.2,
                  box.padding=0.55, point.padding=0.4, max.overlaps=Inf, seed=3,
                  nudge_y=ifelse(lab3$genus=="Halomonas_B",0.7,ifelse(lab3$genus=="Gillisia",0.25,0)),
                  nudge_x=ifelse(lab3$genus=="Halomonas_B",0.15,0)) +
  scale_x_log10(breaks=c(0.002,0.01,0.1,1), labels=c("<0.002","0.01","0.1","1"), limits=c(0.0018,1.2)) +
  scale_y_log10(breaks=c(0.002,0.01,0.1,1,4), labels=c("<0.002","0.01","0.1","1","4"), limits=c(0.0018,4)) +
  scale_fill_manual(values=coh_pal, name="Seasonal cohort") +
  scale_size_continuous(range=c(1.3,5), breaks=c(2,5,8), name="BGCs") +
  labs(x="Summer relative abundance (%)", y="Winter relative abundance (%)", tag="a") +
  guides(fill=guide_legend(order=1, override.aes=list(size=3, colour="white")), size=guide_legend(order=2, override.aes=list(shape=21, fill="#555555", colour="white", stroke=0.4))) +
  th + theme(legend.position="right", legend.title=element_text(size=7,face="bold"),
             legend.text=element_text(size=7), legend.key.size=unit(0.30,"cm"),
             legend.spacing.y=unit(0.05,"cm"), legend.margin=margin(2,2,2,2))

# ---- Panel b: BGC richness across cohorts ----
meds <- dat %>% group_by(cohort) %>% summarise(n=n(), .groups="drop")
pb <- ggplot(dat, aes(cohort, nBGC, fill=cohort)) +
  geom_boxplot(width=0.55, outlier.shape=NA, alpha=0.18, colour="#888888", linewidth=0.4) +
  geom_quasirandom(aes(size=est_gsize_mb), shape=21, stroke=0.4, colour="white", width=0.18, alpha=0.9) +
  geom_text(data=meds, aes(cohort, -0.7, label=paste0("n=",n)), size=2.5,
            colour="#666666", inherit.aes=FALSE) +
  scale_fill_manual(values=coh_pal, guide="none") +
  scale_size_continuous(range=c(1.3,4.2), guide="none") +
  scale_y_continuous(breaks=seq(0,9,3), limits=c(-1.1,9.6)) +
  labs(x=NULL, y="BGCs per genome", tag="b") + th

fig <- (pa | pb) + plot_layout(widths=c(1.5,1)) +
  plot_annotation(caption="Cohort defined by MAG covered-fraction across seasons; point size = BGC count (a) or genome size (b)") &
  theme(plot.caption=element_text(size=6.3, colour="#666666", hjust=1))

dir.create(file.path(B,"figures"), showWarnings=FALSE)
ggsave(file.path(B,"figures/Figure3_cohort.pdf"), fig, width=180, height=92, units="mm", device=cairo_pdf)
ggsave(file.path(B,"figures/Figure3_cohort.png"), fig, width=180, height=92, units="mm", dpi=400)
cat("written: figures/Figure3_cohort.pdf (+ .png)\n")
