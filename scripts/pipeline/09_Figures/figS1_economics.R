suppressWarnings(suppressMessages({ library(ggplot2); library(patchwork) }))
B <- "/data/habib/metagenome/biosurfactant"
d <- read.delim(file.path(B,"integration/bgc_mag_mapping/genome_economics.tsv"))
d$meanAb <- (d$winterAb+d$summerAb)/2
th <- theme_classic(base_size=8, base_family="sans") +
  theme(axis.text=element_text(colour="#444444",size=7), axis.title=element_text(size=8),
        plot.tag=element_text(face="bold",size=11), plot.title=element_text(size=8,hjust=0.5,colour="#333333"),
        plot.margin=margin(6,8,4,6))
sc <- function(df,x,y,xlab,ylab,tag,title,fit=FALSE){
  rho <- suppressWarnings(cor(df[[x]], df[[y]], method="spearman"))
  p   <- suppressWarnings(cor.test(df[[x]], df[[y]], method="spearman")$p.value)
  xmax<-max(df[[x]]); ymax<-max(df[[y]])
  g <- ggplot(df, aes(.data[[x]], .data[[y]]))
  if(fit) g <- g + geom_smooth(method="lm", se=TRUE, colour="#D55E00", fill="#F0D9C9", linewidth=0.5)
  g + geom_point(shape=21, fill="#0072B2", colour="white", size=2, stroke=0.4, alpha=0.85) +
    annotate("text", x=xmax, y=ymax,      hjust=1, vjust=1, label=sprintf("rho = %.2f", rho), size=2.7, colour="#222222") +
    annotate("text", x=xmax, y=ymax*0.88, hjust=1, vjust=1, label=sprintf("P = %.3f", p),     size=2.5, colour="#666666") +
    labs(x=xlab, y=ylab, tag=tag, title=title) + th
}
pa <- sc(d,"est_gsize_mb","nBGC","Genome size (Mb)","BGC richness","a","Richness scales with genome size", fit=TRUE)
pb <- sc(d,"nCDS","nBGC","Predicted CDS","BGC richness","b","Richness scales with coding capacity", fit=TRUE)
pc <- sc(d,"nBGC","coding_density","BGC richness","Coding density","c","No streamlining trade-off")
pd <- sc(d,"est_gsize_mb","meanAb","Genome size (Mb)","Mean relative abundance (%)","d","Genome size unrelated to abundance")
fig <- (pa|pb)/(pc|pd)
ggsave(file.path(B,"figures/FigS3_economics.pdf"), fig, width=150, height=130, units="mm", device=cairo_pdf)
ggsave(file.path(B,"figures/FigS3_economics.png"), fig, width=150, height=130, units="mm", dpi=400)
cat("written FigS3\n")
