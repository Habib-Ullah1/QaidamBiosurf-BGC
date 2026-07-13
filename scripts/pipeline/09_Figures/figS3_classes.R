#!/usr/bin/env Rscript
# Supplementary Figure — BGC class composition across the community inventory
suppressWarnings(suppressMessages({
  need <- c("ggplot2")
  miss <- need[!need %in% rownames(installed.packages())]
  if (length(miss)) install.packages(miss, repos="https://cloud.r-project.org")
  library(ggplot2)
}))
B <- "/data/habib/metagenome/biosurfactant"

# read class counts directly from bgc_regions.tsv (column 4 = product class)
reg <- read.delim(file.path(B,"integration/bgc_mag_mapping/bgc_regions.tsv"),
                  header=FALSE, stringsAsFactors=FALSE)
tab <- as.data.frame(table(reg$V4), stringsAsFactors=FALSE)
colnames(tab) <- c("class","n")

cat_of <- function(cls){
  if (grepl("NRPS|NAPAA|NRP-metallophore", cls)) return("NRPS")
  if (grepl("PKS|PUFA|hglE", cls)) return("PKS")
  if (grepl("RiPP|RRE|lanthipeptide|lassopeptide|thiopeptide|LAP|ranthipeptide|proteusin|thioamitides|linaridin|triceptide|sactipeptide|microviridin|nucleoside", cls)) return("RiPP")
  if (grepl("terpene", cls)) return("Terpene")
  return("Other")
}
tab$category <- sapply(tab$class, cat_of)
tab <- tab[order(tab$n),]
tab$class <- factor(tab$class, levels=tab$class)

cat_pal <- c(NRPS="#0072B2", PKS="#D55E00", RiPP="#009E73", Terpene="#CC79A7", Other="#999999")

p <- ggplot(tab, aes(n, class, fill=category)) +
  geom_col(width=0.72, colour="white", linewidth=0.25) +
  geom_text(aes(label=n), hjust=-0.25, size=2.3, colour="#333333") +
  scale_fill_manual(values=cat_pal, name="Biosynthetic\ncategory") +
  scale_x_continuous(expand=expansion(mult=c(0,0.10)), breaks=c(0,100,200,300,400,500,600)) +
  labs(x="Biosynthetic gene clusters (n)", y=NULL,
       title=paste0("BGC class composition (community inventory, n = ", sum(tab$n), ")")) +
  theme_classic(base_size=8, base_family="sans") +
  theme(axis.text.y=element_text(size=6.5, colour="#333333"),
        axis.text.x=element_text(size=7, colour="#444444"),
        axis.title.x=element_text(size=8),
        plot.title=element_text(size=8, face="bold", hjust=0.5),
        axis.line=element_line(colour="#444444", linewidth=0.4),
        axis.ticks=element_line(colour="#444444", linewidth=0.4),
        legend.position=c(0.82,0.35), legend.title=element_text(size=7,face="bold"),
        legend.text=element_text(size=7), legend.key.size=unit(0.32,"cm"),
        legend.background=element_rect(fill="white", colour="#DDDDDD", linewidth=0.3),
        plot.margin=margin(8,12,6,6))
ggsave(file.path(B,"figures/FigS1_classes.pdf"), p, width=130, height=150, units="mm", device=cairo_pdf)
ggsave(file.path(B,"figures/FigS1_classes.png"), p, width=130, height=150, units="mm", dpi=400)
cat("written: figures/FigS1_classes.pdf (+ .png)\n")
