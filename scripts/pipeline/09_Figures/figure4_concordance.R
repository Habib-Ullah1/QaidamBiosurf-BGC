#!/usr/bin/env Rscript
# Figure 4 — Cross-method concordance: MAG seasonal cohort (metagenome) vs
#            independent 16S seasonal direction (17-sample amplicon survey).
suppressWarnings(suppressMessages({
  need <- c("ggplot2","scales","grid")
  miss <- need[!need %in% rownames(installed.packages())]
  if (length(miss)) install.packages(miss, repos="https://cloud.r-project.org")
  library(ggplot2); library(scales); library(grid)
}))

B <- "/data/habib/metagenome/biosurfactant"
d <- read.delim(file.path(B,"integration/bgc_mag_mapping/concordance_fig.tsv"), stringsAsFactors=FALSE)
d$log2FC <- log2((d$winter16S+1)/(d$summer16S+1))
d$cohort <- factor(d$cohort, levels=c("Winter","Year-round","Summer"))
coh_pal <- c(Winter="#2166AC", `Year-round`="#9A9A9A", Summer="#D55E00")
d <- d[order(d$log2FC),]; d$rank <- seq_len(nrow(d))
flag <- c("Halomonas_B","Moraxella_A"); key <- c("Halomonas_B","Moraxella_A","Loktanella")
d$outline <- ifelse(d$genus %in% flag,"#111111","white")
d$segcol  <- coh_pal[as.character(d$cohort)]
dir <- d[d$cohort %in% c("Winter","Summer"),]
agree <- sum((dir$cohort=="Winter" & dir$log2FC>0) | (dir$cohort=="Summer" & dir$log2FC<0))
n <- nrow(d)

th <- theme_classic(base_size=8, base_family="sans") +
  theme(axis.text.x=element_text(colour="#444444", size=7.5),
        axis.text.y=element_text(colour="#333333", size=7.3, face="italic"),
        axis.title=element_text(colour="#111111", size=8.5),
        axis.line.y=element_blank(), axis.ticks.y=element_blank(),
        axis.line.x=element_line(colour="#444444", linewidth=0.4),
        axis.ticks.x=element_line(colour="#444444", linewidth=0.4),
        legend.position="right", legend.title=element_text(size=7,face="bold"),
        legend.text=element_text(size=7), legend.key.size=unit(0.32,"cm"),
        plot.margin=margin(5,7,3,5),
        plot.caption=element_text(size=6, colour="#777777", hjust=0))

p <- ggplot(d, aes(log2FC, rank)) +
  annotate("rect", xmin=0,  xmax=11.5, ymin=-Inf, ymax=Inf, fill="#2166AC", alpha=0.045) +
  annotate("rect", xmin=-7, xmax=0,    ymin=-Inf, ymax=Inf, fill="#D55E00", alpha=0.045) +
  geom_vline(xintercept=0, colour="#666666", linewidth=0.5) +
  geom_segment(aes(x=0, xend=log2FC, y=rank, yend=rank), colour=d$segcol, linewidth=0.55) +
  geom_point(aes(fill=cohort, size=nBGC), shape=21, stroke=0.6, colour=d$outline) +
  scale_fill_manual(values=coh_pal, name="MAG cohort\n(metagenome)") +
  scale_size_continuous(range=c(1.6,5), breaks=c(2,5,8), name="BGCs") +
  scale_y_continuous(breaks=d$rank, labels=d$genus, expand=expansion(add=c(0.6,1.4))) +
  scale_x_continuous(breaks=seq(-6,10,2), limits=c(-7,11.5)) +
  annotate("text", x=10.5, y=11.5, label="winter-enriched (16S)", hjust=1, size=2.5, colour="#2166AC", fontface="italic") +
  annotate("text", x=-6.8, y=4.5,  label="summer-enriched (16S)", hjust=0, size=2.5, colour="#D55E00", fontface="italic") +
  annotate("text", x=-6.6, y=n+0.9, label=paste0(agree,"/9 directional MAGs agree with 16S"),
           hjust=0, vjust=1, size=2.7, colour="#222222", fontface="bold") +
  guides(fill=guide_legend(order=1, override.aes=list(size=3)), size=guide_legend(order=2)) +
  labs(x=expression("16S seasonal direction   "*log[2]*"(winter / summer)"), y=NULL,
       caption="Black outline = flagship MAG.\n16S: independent 17-sample survey; MAG cohort: 4-sample metagenome.") + th

dir.create(file.path(B,"figures"), showWarnings=FALSE)
ggsave(file.path(B,"figures/Figure4_concordance.pdf"), p, width=150, height=100, units="mm", device=cairo_pdf)
ggsave(file.path(B,"figures/Figure4_concordance.png"), p, width=150, height=100, units="mm", dpi=400)
cat("written: figures/Figure4_concordance.pdf  (concordance", agree, "/9)\n")
