suppressWarnings(suppressMessages({
  need <- c("ggplot2","ggbeeswarm","ggrepel","patchwork","dplyr","scales","grid")
  miss <- need[!need %in% rownames(installed.packages())]
  if (length(miss)) install.packages(miss, repos="https://cloud.r-project.org")
  library(ggplot2); library(ggbeeswarm); library(ggrepel)
  library(patchwork); library(dplyr); library(scales); library(grid)
}))
B <- "/data/habib/metagenome/biosurfactant"
dat <- read.delim(file.path(B,"integration/bgc_mag_mapping/genome_economics.tsv"), stringsAsFactors=FALSE)
major <- c("Pseudomonadota","Bacteroidota","Chloroflexota","Halobacteriota","Actinomycetota","Myxococcota","Deinococcota")
pal <- c(Pseudomonadota="#0072B2", Bacteroidota="#E69F00", Chloroflexota="#009E73", Halobacteriota="#CC79A7",
         Actinomycetota="#D55E00", Myxococcota="#56B4E9", Deinococcota="#000000", `other / unclassified`="#BFBFBF")
dat$phylum_grp <- factor(ifelse(dat$phylum %in% major, dat$phylum, "other / unclassified"), levels=names(pal))
flag_main <- c("H2_concoct_7","H6_concoct_8")
dat$flag <- factor(ifelse(dat$MAG %in% flag_main,"flagship","other"), levels=c("other","flagship"))
sA <- suppressWarnings(cor.test(dat$nBGC, dat$winterAb, method="spearman", exact=FALSE))
sB <- suppressWarnings(cor.test(dat$nBGC, dat$est_gsize_mb, method="spearman", exact=FALSE))
labA<-sprintf("rho==%.2f",sA$estimate); pAt<-ifelse(sA$p.value<0.001,"P<0.001",sprintf("P==%.2f",sA$p.value))
labB<-sprintf("rho==%.2f",sB$estimate); pBt<-ifelse(sB$p.value<0.001,"P<0.001",sprintf("P==%.2f",sB$p.value))
gm  <- 10^mean(log10(dat$winterAb))
lab_taxa <- dat %>% filter(MAG %in% c("C3_concoct_26","H2_concoct_31","H2_concoct_7","H6_concoct_8")) %>%
  mutate(lab=c(C3_concoct_26="Gillisia",H2_concoct_31="Ralstonia",H2_concoct_7="Halomonas_B",H6_concoct_8="Moraxella_A")[MAG],
         nx=c(C3_concoct_26=4.4,H2_concoct_31=9.2,H2_concoct_7=6.5,H6_concoct_8=1.3)[MAG],
         ny=c(C3_concoct_26=3.7,H2_concoct_31=0.045,H2_concoct_7=0.9,H6_concoct_8=1.8)[MAG])
th <- theme_classic(base_size=8, base_family="sans") +
  theme(axis.text=element_text(colour="#444444", size=7.5), axis.title=element_text(colour="#111111", size=8.5),
        axis.line=element_line(colour="#444444", linewidth=0.4), axis.ticks=element_line(colour="#444444", linewidth=0.4),
        plot.tag=element_text(face="bold", size=12), legend.position="none", plot.margin=margin(5,7,3,5))
set.seed(42)

# PANEL A (was b): BGC richness vs genome size  -- cited FIRST in text
pa <- ggplot(dat, aes(est_gsize_mb, nBGC)) +
  geom_smooth(method="lm", se=FALSE, colour="#333333", linewidth=0.55, alpha=0.6) +
  geom_point(aes(fill=phylum_grp, colour=flag), position=position_jitter(height=0.18, width=0, seed=42), shape=21, size=2.3, stroke=0.5, alpha=0.85) +
  annotate("label", x=6.6, y=9.7, hjust=1, vjust=1, parse=TRUE, size=2.9, colour="#222222", fontface="bold", label.size=0, fill=NA, label=labB) +
  annotate("text", x=6.6, y=8.5, hjust=1, vjust=1, size=2.6, colour="#222222", label=pBt) +
  scale_x_continuous(breaks=1:6, limits=c(0.8,6.7)) + scale_y_continuous(breaks=seq(0,9,3), limits=c(-0.6,10)) +
  scale_fill_manual(values=pal, drop=FALSE, name=NULL) +
  scale_colour_manual(values=c(other="white", flagship="#111111"), guide="none") +
  labs(x="Estimated genome size (Mb)", y="BGCs per genome", tag="a") +
  guides(fill=guide_legend(order=1, nrow=2, byrow=TRUE, override.aes=list(size=2.7, colour="white"))) +
  th

# PANEL B (was a): BGC richness vs winter abundance  -- cited SECOND in text
pb <- ggplot(dat, aes(nBGC, winterAb)) +
  geom_quasirandom(aes(fill=phylum_grp, colour=flag, size=est_gsize_mb), shape=21, width=0.30, stroke=0.5, alpha=0.85) +
  geom_text_repel(data=lab_taxa, aes(nx, ny, label=lab), fontface="italic", size=2.9, segment.colour="#777777",
                  segment.size=0.3, min.segment.length=0.25, box.padding=0.3, point.padding=0.4, max.overlaps=Inf, seed=7) +
  annotate("label", x=9.5, y=4.6, hjust=1, vjust=1, parse=TRUE, size=2.9, colour="#C0392B", fontface="bold", label.size=0, fill=NA, label=labA) +
  annotate("text", x=9.5, y=3.05, hjust=1, vjust=1, parse=TRUE, size=2.6, colour="#C0392B", label=pAt) +
  annotate("text", x=9.5, y=2.32, hjust=1, vjust=1, size=2.4, colour="#888888", label="Spearman, n=44") +
  scale_y_log10(breaks=c(0.01,0.1,1,4), labels=c("0.01","0.1","1","4"), limits=c(0.006,4.9), expand=expansion(mult=c(0.02,0.03))) +
  scale_x_continuous(breaks=0:9, limits=c(-0.4,9.7)) +
  scale_fill_manual(values=pal, drop=FALSE, guide="none") +
  scale_colour_manual(values=c(other="white", flagship="#111111"), guide="none") +
  scale_size_continuous(range=c(1.4,4.6), breaks=c(2,4,6), name="Genome size (Mb)") +
  labs(x="Biosynthetic gene clusters per genome", y="Winter relative abundance (%)", tag="b") +
  guides(size=guide_legend(order=2, nrow=1)) +
  th + annotation_logticks(sides="l", colour="#AAAAAA", long=unit(0.10,"cm"), mid=unit(0.06,"cm"), short=unit(0.04,"cm"))

# layout: a (genome size) LEFT, b (abundance) RIGHT
fig <- (pa | pb) + plot_layout(widths=c(1,1.3), guides="collect") +
  plot_annotation(caption="Black outline = high-quality flagship MAG") &
  theme(legend.position="bottom", legend.box="horizontal", legend.box.just="center",
        legend.title=element_text(size=6.8, face="bold"), legend.text=element_text(size=6.8),
        legend.key.size=unit(0.30,"cm"), legend.margin=margin(2,2,2,2), legend.spacing.x=unit(0.05,"cm"),
        legend.box.spacing=unit(0.1,"cm"), plot.caption=element_text(size=6.3, colour="#666666", hjust=1))
dir.create(file.path(B,"figures"), showWarnings=FALSE)
ggsave(file.path(B,"figures/Figure2_decoupling.pdf"), fig, width=180, height=105, units="mm", device=cairo_pdf)
ggsave(file.path(B,"figures/Figure2_decoupling.png"), fig, width=180, height=105, units="mm", dpi=400)
cat("rendered v4 (panels swapped: a=genome size, b=abundance)\n")
