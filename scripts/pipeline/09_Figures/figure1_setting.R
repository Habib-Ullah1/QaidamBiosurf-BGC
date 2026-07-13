#!/usr/bin/env Rscript
# Figure 1 — Study system: 16S seasonal community turnover + MAG recovery
suppressWarnings(suppressMessages({
  need <- c("ggplot2","ggrepel","patchwork","dplyr","grid")
  miss <- need[!need %in% rownames(installed.packages())]
  if (length(miss)) install.packages(miss, repos="https://cloud.r-project.org")
  library(ggplot2); library(ggrepel); library(patchwork); library(dplyr); library(grid)
}))
B <- "/data/habib/metagenome/biosurfactant"
pc  <- read.delim(file.path(B,"integration/bgc_mag_mapping/pcoa_season.tsv"), stringsAsFactors=FALSE)
rec <- read.delim(file.path(B,"integration/bgc_mag_mapping/mag_recovery.tsv"), stringsAsFactors=FALSE)
sea_pal <- c(Summer="#D55E00", Winter="#2166AC")
R2 <- 0.193; pval <- 0.001; ax1 <- 23.0; ax2 <- 19.1   # PERMANOVA + PCoA axis variance (from QIIME2 output)

th <- theme_classic(base_size=8, base_family="sans") +
  theme(axis.text=element_text(colour="#444444", size=7.3), axis.title=element_text(colour="#111111", size=8.3),
        axis.line=element_line(colour="#444444", linewidth=0.4), axis.ticks=element_line(colour="#444444", linewidth=0.4),
        plot.tag=element_text(face="bold", size=12), plot.title=element_text(size=8.3,face="bold"),
        legend.position="none", plot.margin=margin(10,8,4,6))

# ---- Panel a: PCoA seasonal turnover ----
hulls <- pc %>% group_by(season) %>% slice(chull(PCo1, PCo2))
pa <- ggplot(pc, aes(PCo1, PCo2)) +
  geom_polygon(data=hulls, aes(fill=season, group=season), alpha=0.12, colour=NA) +
  geom_vline(xintercept=0, colour="#E8E8E8", linewidth=0.4) + geom_hline(yintercept=0, colour="#E8E8E8", linewidth=0.4) +
  geom_point(aes(fill=season, shape=type), size=2.6, stroke=0.5, colour="white") +
  scale_shape_manual(values=c(Sediment=21, Water=24, Rock=22), name="Sample type") +
  scale_fill_manual(values=sea_pal, name="Season") +
  annotate("text", x=-0.30, y=0.62, label=paste0("PERMANOVA\nR\u00b2 = ", R2, ", P = ", pval),
           hjust=0, size=2.7, colour="#222222", fontface="bold", lineheight=0.95) +
  annotate("text", x=0.34, y=0.30, label="Winter", colour="#2166AC", fontface="bold", size=3) +
  annotate("text", x=-0.32, y=-0.27, label="Summer", colour="#D55E00", fontface="bold", size=3) +
  labs(title="16S community turnover", x=paste0("PCo1 (",ax1,"%)"), y=paste0("PCo2 (",ax2,"%)"), tag="a") +
  guides(fill=guide_legend(order=1, override.aes=list(shape=21, size=3)),
         shape=guide_legend(order=2, override.aes=list(fill="#888888"))) +
  th + theme(legend.position="right", legend.title=element_text(size=7,face="bold"),
             legend.text=element_text(size=7), legend.key.size=unit(0.32,"cm"))

# ---- Panel b: MAG recovery (in-panel title + manual legend to avoid patchwork clipping) ----
rec$sample <- factor(rec$sample, levels=c("C3","C6","H2","H6"))
recl <- rbind(
  data.frame(sample=rec$sample, season=rec$season, cat="Total bins", n=rec$bins),
  data.frame(sample=rec$sample, season=rec$season, cat=">=70% complete", n=rec$hq70))
recl$cat <- factor(recl$cat, levels=c("Total bins",">=70% complete"))
ymax <- max(rec$bins)*1.18
pb <- ggplot(recl, aes(sample, n, fill=season, alpha=cat)) +
  geom_col(position="identity", width=0.66, colour="white", linewidth=0.4) +
  geom_text(data=rec, aes(sample, bins, label=bins), vjust=-0.5, size=2.6, colour="#555555", inherit.aes=FALSE) +
  geom_text(data=rec, aes(sample, hq70, label=hq70), vjust=-0.5, size=2.6, colour="white", fontface="bold", inherit.aes=FALSE) +
  scale_fill_manual(values=sea_pal, guide="none") +
  scale_alpha_manual(values=c("Total bins"=0.38, ">=70% complete"=1), guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0,0.02)), limits=c(0,ymax)) +
  annotate("text", x=0.5, y=ymax*0.99, label="MAG recovery", hjust=0, fontface="bold", size=3.0, colour="#111111") +
  annotate("rect", xmin=1.7, xmax=1.9, ymin=ymax*0.72, ymax=ymax*0.78, fill="#888888", alpha=0.38) +
  annotate("text", x=1.97, y=ymax*0.75, label="total bins", hjust=0, size=2.3, colour="#444444") +
  annotate("rect", xmin=1.7, xmax=1.9, ymin=ymax*0.61, ymax=ymax*0.67, fill="#888888", alpha=1) +
  annotate("text", x=1.97, y=ymax*0.64, label="\u226570% complete", hjust=0, size=2.3, colour="#444444") +
  labs(x=NULL, y="Metagenome bins", tag="b") + th + theme(plot.title=element_blank())

fig <- (pa | pb) + plot_layout(widths=c(1.5,1)) +
  plot_annotation(caption="17 samples (11 summer, 6 winter); Bray-Curtis PCoA. Bar shading: light = total bins, solid = \u226570% complete.",
                  theme=theme(plot.caption=element_text(size=6.2, colour="#666666", hjust=0, family="sans")))

dir.create(file.path(B,"figures"), showWarnings=FALSE)
ggsave(file.path(B,"figures/Figure1_setting.pdf"), fig, width=180, height=88, units="mm", device=cairo_pdf)
ggsave(file.path(B,"figures/Figure1_setting.png"), fig, width=180, height=88, units="mm", dpi=400)
cat("written: figures/Figure1_setting.pdf (+ .png)\n")
