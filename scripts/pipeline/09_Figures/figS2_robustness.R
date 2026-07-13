suppressWarnings(suppressMessages({ library(ggplot2); library(patchwork) }))
B <- "/data/habib/metagenome/biosurfactant"
d <- read.delim(file.path(B,"integration/bgc_mag_mapping/genome_economics.tsv"))
d$meanAb <- (d$winterAb + d$summerAb)/2
conc <- read.delim(file.path(B,"integration/bgc_mag_mapping/concordance_fig.tsv"))
th <- theme_classic(base_size=8, base_family="sans") +
  theme(axis.text=element_text(colour="#444444",size=7), axis.title=element_text(size=8),
        plot.tag=element_text(face="bold",size=11), plot.title=element_text(size=8,hjust=0.5,colour="#333333"),
        plot.margin=margin(6,8,4,6))
sc <- function(df,x,y,xlab,tag,title){
  df <- df[is.finite(df[[x]]) & is.finite(df[[y]]),]
  rho <- suppressWarnings(cor(df[[x]], df[[y]], method="spearman"))
  p   <- suppressWarnings(cor.test(df[[x]], df[[y]], method="spearman")$p.value)
  xmax<-max(df[[x]]); ymax<-max(df[[y]])
  ggplot(df, aes(.data[[x]], .data[[y]])) +
    geom_point(shape=21, fill="#0072B2", colour="white", size=2, stroke=0.4, alpha=0.85) +
    annotate("text", x=xmax, y=ymax,      hjust=1, vjust=1, label=sprintf("rho = %.2f", rho), size=2.7, colour="#222222") +
    annotate("text", x=xmax, y=ymax*0.88, hjust=1, vjust=1, label=sprintf("P = %.2f", p),     size=2.5, colour="#666666") +
    labs(x=xlab, y="BGC richness", tag=tag, title=title) + th
}
p1 <- sc(d,"winterAb","nBGC","Winter relative abundance (%)","a","vs winter abundance")
p2 <- sc(d,"summerAb","nBGC","Summer relative abundance (%)","b","vs summer abundance")
p3 <- sc(d,"meanAb","nBGC","Mean relative abundance (%)","c","vs mean abundance")
p4 <- sc(conc,"winter16S","nBGC","16S winter abundance","d","vs independent 16S abundance")
fig <- (p1|p2)/(p3|p4)
ggsave(file.path(B,"figures/FigS5_robustness.pdf"), fig, width=150, height=130, units="mm", device=cairo_pdf)
ggsave(file.path(B,"figures/FigS5_robustness.png"), fig, width=150, height=130, units="mm", dpi=400)
cat("written FigS5\n")
