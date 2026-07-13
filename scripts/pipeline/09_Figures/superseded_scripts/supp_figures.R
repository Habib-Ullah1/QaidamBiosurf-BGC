suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(ggrepel)
  library(scales)
  library(ape)
  library(pheatmap)
})

PROJECT <- "/data/habib/metagenome"
OUTDIR  <- file.path(PROJECT,"biosurfactant/figures")

pub_theme <- theme_classic(base_size=10,base_family="sans") +
  theme(
    axis.text=element_text(size=8,color="black"),
    axis.title=element_text(size=9,color="black",face="bold"),
    legend.text=element_text(size=8),
    legend.title=element_text(size=9,face="bold"),
    plot.title=element_text(size=9,face="bold",hjust=0),
    plot.subtitle=element_text(size=7,color="#555555"),
    panel.grid.major=element_line(color="#EEEEEE",linewidth=0.3),
    legend.key.size=unit(0.4,"cm"),
    strip.background=element_rect(fill="#2C3E50",color=NA),
    strip.text=element_text(color="white",face="bold",size=8)
  )

sc <- c("Summer"="#E8593C","Winter"="#3B8BD4")

# ============================================================
# SUPP FIGURE S1 — Assembly contig length distribution
# ============================================================
cat("Building Supplementary Figure S1...\n")

contig_data <- data.frame(
  Sample=factor(rep(c("C3","C6","H2","H6"),each=5),
                levels=c("C3","C6","H2","H6")),
  Season=factor(rep(c("Summer","Summer","Winter","Winter"),each=5),
                levels=c("Summer","Winter")),
  LengthBin=factor(rep(c("500-999","1000-1999","2000-4999",
                          "5000-9999",">=10000"),4),
                   levels=c("500-999","1000-1999","2000-4999",
                             "5000-9999",">=10000")),
  Count=c(
    926891,245598,69972,8984,2755,   # C3
    438138,95129,23861,2823,751,     # C6
    711556,184288,56040,9055,3084,   # H2
    167152,35993,9362,1521,474       # H6
  )
)

contig_pct <- contig_data %>%
  group_by(Sample) %>%
  mutate(Pct=Count/sum(Count)*100) %>%
  ungroup()

bin_colors <- c(
  "500-999"   ="#4292C6",
  "1000-1999" ="#2171B5",
  "2000-4999" ="#08519C",
  "5000-9999" ="#F16913",
  ">=10000"   ="#D94801"
)

pS1a <- ggplot(contig_pct,
               aes(x=Sample,y=Count/1000,fill=LengthBin)) +
  geom_bar(stat="identity",width=0.7,
           color="white",linewidth=0.2) +
  scale_fill_manual(values=bin_colors,name="Contig length (bp)") +
  scale_y_continuous(expand=expansion(mult=c(0,0.05)),
                     labels=comma) +
  facet_grid(~Season,scales="free_x",space="free_x") +
  labs(x=NULL,y="Number of contigs (x1000)",
       title="A. Contig count by length bin") +
  pub_theme +
  theme(legend.position="right")

pS1b <- ggplot(contig_pct,
               aes(x=Sample,y=Pct,fill=LengthBin)) +
  geom_bar(stat="identity",width=0.7,
           color="white",linewidth=0.2) +
  scale_fill_manual(values=bin_colors,name="Contig length (bp)") +
  scale_y_continuous(expand=expansion(mult=c(0,0.02)),
                     labels=function(x) paste0(x,"%")) +
  facet_grid(~Season,scales="free_x",space="free_x") +
  labs(x=NULL,y="Percentage of contigs (%)",
       title="B. Proportional contig length distribution") +
  pub_theme +
  theme(legend.position="right")

# N50 comparison
n50_data <- data.frame(
  Sample=factor(c("C3","C6","H2","H6"),
                levels=c("C3","C6","H2","H6")),
  Season=c("Summer","Summer","Winter","Winter"),
  N50=c(1010,891,1053,925),
  TotalMb=c(1259,515,1000,205)
)

pS1c <- ggplot(n50_data,
               aes(x=Sample,y=N50,fill=Season)) +
  geom_bar(stat="identity",width=0.6,
           color="white",linewidth=0.3) +
  geom_text(aes(label=paste0(N50," bp")),
            vjust=-0.4,size=3,fontface="bold") +
  scale_fill_manual(values=sc,name="Season") +
  scale_y_continuous(expand=expansion(mult=c(0,0.15))) +
  labs(x=NULL,y="N50 (bp)",
       title="C. Assembly N50 per sample") +
  pub_theme +
  theme(legend.position="none")

pS1d <- ggplot(n50_data,
               aes(x=Sample,y=TotalMb,fill=Season)) +
  geom_bar(stat="identity",width=0.6,
           color="white",linewidth=0.3) +
  geom_text(aes(label=paste0(TotalMb," Mb")),
            vjust=-0.4,size=3,fontface="bold") +
  scale_fill_manual(values=sc,name="Season") +
  scale_y_continuous(expand=expansion(mult=c(0,0.15))) +
  labs(x=NULL,y="Total assembly length (Mb)",
       title="D. Total assembly size per sample") +
  pub_theme +
  theme(legend.position="bottom",
        legend.direction="horizontal")

figS1 <- (pS1a | pS1b) / (pS1c | pS1d) +
  plot_layout(heights=c(1.5,1))

ggsave(file.path(OUTDIR,"FigureS1_assembly_stats.pdf"),
       figS1,width=170,height=160,units="mm",dpi=300)
ggsave(file.path(OUTDIR,"FigureS1_assembly_stats.png"),
       figS1,width=170,height=160,units="mm",dpi=300)
cat("Supplementary Figure S1 done.\n")

# ============================================================
# SUPP FIGURE S2 — BiG-SCAPE GCF network summary
# (BiG-SCAPE network image not directly extractable;
#  we create a publication-quality network statistics figure)
# ============================================================
cat("Building Supplementary Figure S2...\n")

# GCF size distribution
gcf_sizes <- data.frame(
  Size_Category=factor(c("Singleton\n(1 BGC)",
                          "Small\n(2-3 BGCs)",
                          "Medium\n(4-9 BGCs)",
                          "Large\n(>=10 BGCs)"),
                       levels=c("Singleton\n(1 BGC)",
                                 "Small\n(2-3 BGCs)",
                                 "Medium\n(4-9 BGCs)",
                                 "Large\n(>=10 BGCs)")),
  GCF_Count=c(376,178,30,6),
  Pct=c(63.7,30.2,5.1,1.0)
)

pS2a <- ggplot(gcf_sizes,
               aes(x=Size_Category,y=GCF_Count,
                   fill=Size_Category)) +
  geom_bar(stat="identity",width=0.7,
           color="white",linewidth=0.3) +
  geom_text(aes(label=paste0(GCF_Count,"\n(",Pct,"%)")),
            vjust=-0.3,size=3,fontface="bold") +
  scale_fill_manual(
    values=c("#C0392B","#E67E22","#2980B9","#27AE60"),
    guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0,0.18))) +
  labs(x="GCF size (number of member BGCs)",
       y="Number of GCFs",
       title="A. GCF size distribution (n=590 GCFs, c=0.50)") +
  pub_theme

# BGC class per GCF cross-season
gcf_class <- data.frame(
  Class=rep(c("NRPS","PKS other","RiPPs",
               "Terpene","Others","PKS I"),3),
  Season=factor(rep(c("Summer-only",
                       "Cross-season",
                       "Winter-only"),each=6),
                levels=c("Summer-only","Cross-season","Winter-only")),
  Count=c(22,22,75,58,44,5,
          3,11,17,51,12,2,
          35,31,52,47,71,3)
)

season_col <- c("Summer-only"="#E8593C",
                "Cross-season"="#7D9B76",
                "Winter-only"="#3B8BD4")

pS2b <- ggplot(gcf_class,
               aes(x=Class,y=Count,fill=Season)) +
  geom_bar(stat="identity",
           position=position_dodge(width=0.8),
           width=0.7,color="white",linewidth=0.2) +
  scale_fill_manual(values=season_col,name="Season") +
  scale_y_continuous(expand=expansion(mult=c(0,0.12))) +
  labs(x="BGC class",y="Number of GCFs",
       title="B. Seasonal distribution of GCFs by BGC class") +
  pub_theme +
  theme(legend.position="bottom",
        legend.direction="horizontal")

# Novelty by class
novelty_class <- data.frame(
  Class=factor(c("NRPS","PKS other","RiPPs",
                  "Terpene","Others","PKS I"),
               levels=c("NRPS","PKS other","RiPPs",
                         "Terpene","Others","PKS I")),
  Singleton=c(38,42,102,108,75,8),
  Multi_novel=c(19,20,40,43,34,2),
  MIBiG_related=c(3,2,2,5,1,0)
)

nov_long <- novelty_class %>%
  pivot_longer(cols=c(Singleton,Multi_novel,MIBiG_related),
               names_to="Category",values_to="Count") %>%
  mutate(Category=factor(
    recode(Category,
           Singleton="Singleton (novel)",
           Multi_novel="Multi-member (novel)",
           MIBiG_related="MIBiG-related"),
    levels=c("MIBiG-related","Multi-member (novel)",
              "Singleton (novel)")))

pS2c <- ggplot(nov_long,
               aes(x=Class,y=Count,fill=Category)) +
  geom_bar(stat="identity",width=0.7,
           color="white",linewidth=0.2) +
  scale_fill_manual(
    values=c("MIBiG-related"="#27AE60",
             "Multi-member (novel)"="#E67E22",
             "Singleton (novel)"="#C0392B"),
    name="GCF novelty") +
  scale_y_continuous(expand=expansion(mult=c(0,0.08))) +
  labs(x="BGC class",y="Number of GCFs",
       title="C. Novelty composition per BGC class") +
  pub_theme +
  theme(legend.position="bottom",
        legend.direction="horizontal")

# Cosine similarity distribution
# At cutoffs 0.30, 0.50, 0.70
cutoff_data <- data.frame(
  Cutoff=factor(c("c=0.30","c=0.50","c=0.70"),
                levels=c("c=0.30","c=0.50","c=0.70")),
  GCFs=c(520,590,623),
  Singletons=c(312,376,441),
  MultiMember=c(208,214,182)
)

cutoff_long <- cutoff_data %>%
  pivot_longer(cols=c(Singletons,MultiMember),
               names_to="Type",values_to="Count")

pS2d <- ggplot(cutoff_long,
               aes(x=Cutoff,y=Count,fill=Type)) +
  geom_bar(stat="identity",width=0.6,
           color="white",linewidth=0.3) +
  scale_fill_manual(
    values=c("Singletons"="#C0392B",
             "MultiMember"="#3B8BD4"),
    labels=c("Multi-member GCFs","Singleton GCFs"),
    name=NULL) +
  scale_y_continuous(expand=expansion(mult=c(0,0.1))) +
  labs(x="BiG-SCAPE cosine similarity cutoff",
       y="Number of GCFs",
       title="D. GCF count at three similarity cutoffs") +
  pub_theme +
  theme(legend.position="bottom",
        legend.direction="horizontal")

figS2 <- (pS2a | pS2b) / (pS2c | pS2d)

ggsave(file.path(OUTDIR,"FigureS2_bigscape_network.pdf"),
       figS2,width=170,height=160,units="mm",dpi=300)
ggsave(file.path(OUTDIR,"FigureS2_bigscape_network.png"),
       figS2,width=170,height=160,units="mm",dpi=300)
cat("Supplementary Figure S2 done.\n")

# ============================================================
# SUPP FIGURE S3 — DESeq2 MA plot
# ============================================================
cat("Building Supplementary Figure S3...\n")

deseq_all <- read.table(
  file.path(PROJECT,"biosurfactant/14_deseq2/deseq2_winter_vs_summer_all.tsv"),
  header=TRUE,sep="\t",stringsAsFactors=FALSE)
deseq_all <- deseq_all[!is.na(deseq_all$padj),]
deseq_all$category <- "Not significant"
deseq_all$category[deseq_all$padj<0.05 & deseq_all$log2FoldChange < -1] <- "Summer enriched"
deseq_all$category[deseq_all$padj<0.05 & deseq_all$log2FoldChange >  1] <- "Winter enriched"
deseq_all$category <- factor(deseq_all$category,
  levels=c("Summer enriched","Winter enriched","Not significant"))

n_s <- sum(deseq_all$category=="Summer enriched")
n_w <- sum(deseq_all$category=="Winter enriched")

# MA plot: log2FC vs log10(baseMean)
pS3a <- ggplot(deseq_all,
               aes(x=log10(baseMean+1),
                   y=log2FoldChange,
                   color=category)) +
  geom_point(alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(-1,0,1),
             linetype=c("dashed","solid","dashed"),
             color=c("#E8593C","black","#3B8BD4"),
             linewidth=c(0.5,0.8,0.5),alpha=0.7) +
  scale_color_manual(
    values=c("Summer enriched"="#E8593C",
             "Winter enriched"="#3B8BD4",
             "Not significant"="#CCCCCC"),
    name=NULL,
    labels=c(paste0("Summer enriched (n=",n_s,")"),
             paste0("Winter enriched (n=",n_w,")"),
             "Not significant")) +
  labs(x=expression(log[10]*"(mean normalised coverage + 1)"),
       y=expression(log[2]*" Fold Change (Winter/Summer)"),
       title="A. MA plot — biosurfactant contig differential abundance",
       subtitle="Horizontal dashed lines at log2FC = +/-1 (biological significance threshold)") +
  guides(color=guide_legend(override.aes=list(size=2.5,alpha=1))) +
  pub_theme +
  theme(legend.position=c(0.75,0.15),
        legend.background=element_rect(fill="white",
                                        color="grey80",
                                        linewidth=0.3),
        legend.text=element_text(size=7))

# P-value distribution
pS3b <- ggplot(deseq_all[!is.na(deseq_all$pvalue),],
               aes(x=pvalue)) +
  geom_histogram(bins=40,fill="#2980B9",
                 color="white",linewidth=0.2) +
  geom_vline(xintercept=0.05,linetype="dashed",
             color="red",linewidth=0.6) +
  scale_x_continuous(expand=expansion(mult=c(0,0.02))) +
  scale_y_continuous(expand=expansion(mult=c(0,0.05)),
                     labels=comma) +
  labs(x="Raw p-value",y="Number of contigs",
       title="B. Raw p-value distribution",
       subtitle="Uniform distribution with enrichment near 0 indicates valid test") +
  pub_theme

# log2FC distribution for significant contigs
sig_only <- deseq_all[deseq_all$padj<0.05,]
pS3c <- ggplot(sig_only,aes(x=log2FoldChange,fill=category)) +
  geom_histogram(bins=50,color="white",linewidth=0.1) +
  scale_fill_manual(
    values=c("Summer enriched"="#E8593C",
             "Winter enriched"="#3B8BD4"),
    name=NULL,guide="none") +
  geom_vline(xintercept=0,color="black",linewidth=0.6) +
  scale_x_continuous(breaks=seq(-30,30,10)) +
  scale_y_continuous(expand=expansion(mult=c(0,0.05))) +
  labs(x=expression(log[2]*" Fold Change"),
       y="Number of significant contigs",
       title="C. log2FC distribution of significant contigs",
       subtitle=paste0("n=",nrow(sig_only)," total significant (padj<0.05)")) +
  annotate("text",x=-20,y=Inf,
           label=paste0("Summer\nn=",n_s),
           color="#E8593C",vjust=1.2,
           size=3,fontface="bold") +
  annotate("text",x=15,y=Inf,
           label=paste0("Winter\nn=",n_w),
           color="#3B8BD4",vjust=1.2,
           size=3,fontface="bold") +
  pub_theme

# baseMean comparison between enrichment groups
pS3d <- ggplot(sig_only,
               aes(x=category,y=log10(baseMean+1),
                   fill=category)) +
  geom_boxplot(width=0.5,color="black",
               outlier.size=0.5,outlier.alpha=0.3) +
  scale_fill_manual(
    values=c("Summer enriched"="#E8593C",
             "Winter enriched"="#3B8BD4"),
    guide="none") +
  labs(x=NULL,
       y=expression(log[10]*"(mean coverage + 1)"),
       title="D. Mean coverage distribution by enrichment group",
       subtitle="Higher baseMean = consistently higher abundance across samples") +
  pub_theme +
  theme(axis.text.x=element_text(size=8))

figS3 <- (pS3a | pS3b) / (pS3c | pS3d)

ggsave(file.path(OUTDIR,"FigureS3_deseq2_supplement.pdf"),
       figS3,width=170,height=160,units="mm",dpi=300)
ggsave(file.path(OUTDIR,"FigureS3_deseq2_supplement.png"),
       figS3,width=170,height=160,units="mm",dpi=300)
cat("Supplementary Figure S3 done.\n")

# ============================================================
# SUPP FIGURE S4 — MAG phylogenetic tree from GTDB-Tk
# ============================================================
cat("Building Supplementary Figure S4...\n")

# Build from GTDB taxonomy strings rather than tree file
# (tree files are per-sample, not consolidated)
mag_tax <- data.frame(
  MAG=c("H6_concoct_8","H2_concoct_0","H2_concoct_109",
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
  Phylum=c("Pseudomonadota","Bacillota_A","Nanohaloarchaeota",
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
  Domain=c("Bacteria","Bacteria","Archaea",
           "Bacteria","Bacteria","Bacteria",
           "Bacteria","Bacteria","Archaea",
           "Bacteria","Archaea","Bacteria",
           "Bacteria","Bacteria","Bacteria",
           "Bacteria","Bacteria","Bacteria",
           "Archaea","Bacteria","Bacteria",
           "Bacteria","Bacteria","Archaea",
           "Bacteria","Bacteria","Bacteria",
           "Bacteria","Bacteria","Bacteria",
           "Bacteria","Bacteria","Bacteria",
           "Bacteria","Bacteria","Bacteria",
           "Bacteria","Archaea","Bacteria",
           "Archaea","Archaea","Bacteria",
           "Archaea","Archaea"),
  Completeness=c(99.99,98.87,99.46,97.68,93.23,91.62,
                 96.78,91.17,97.41,97.37,85.57,87.31,
                 88.29,88.68,88.52,78.89,82.94,85.13,
                 80.95,77.83,78.82,77.19,77.32,79.86,
                 52.81,69.40,72.94,71.40,71.88,70.31,
                 62.29,76.53,62.54,63.65,91.74,59.09,
                 57.06,58.88,51.04,51.05,65.97,52.73,
                 53.34,86.69),
  Season=c("Winter","Winter","Winter","Winter","Winter","Winter",
           "Summer","Summer","Winter","Winter","Winter","Winter",
           "Winter","Winter","Winter","Winter","Summer","Winter",
           "Winter","Winter","Winter","Winter","Summer","Winter",
           "Winter","Summer","Summer","Summer","Winter","Winter",
           "Summer","Summer","Summer","Summer","Summer","Summer",
           "Winter","Summer","Summer","Winter","Winter","Summer",
           "Summer","Winter"),
  Flagship=c(TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,
             FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,
             FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,
             FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,
             FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,
             FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,
             FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,
             FALSE,FALSE),
  stringsAsFactors=FALSE
)

phylum_colors2 <- c(
  "Pseudomonadota"    ="#E74C3C",
  "Chloroflexota"     ="#27AE60",
  "Halobacteriota"    ="#9B59B6",
  "Nanohaloarchaeota" ="#8E44AD",
  "Bacillota_A"       ="#F39C12",
  "Bacteroidota"      ="#3498DB",
  "Actinomycetota"    ="#1ABC9C",
  "Gemmatimonadota"   ="#E67E22",
  "Myxococcota"       ="#C0392B",
  "Deinococcota"      ="#7F8C8D",
  "Patescibacteria"   ="#BDC3C7"
)

# Phylum summary bar
phylum_summary <- mag_tax %>%
  count(Phylum,Domain) %>%
  arrange(Domain,desc(n))

phylum_summary$Phylum <- factor(phylum_summary$Phylum,
  levels=phylum_summary$Phylum)

pS4a <- ggplot(phylum_summary,
               aes(x=n,y=reorder(Phylum,n),
                   fill=Phylum)) +
  geom_bar(stat="identity",width=0.7,
           color="white",linewidth=0.3) +
  geom_text(aes(label=n),hjust=-0.3,size=3) +
  scale_fill_manual(values=phylum_colors2,guide="none") +
  scale_x_continuous(expand=expansion(mult=c(0,0.2))) +
  facet_grid(Domain~.,scales="free_y",space="free_y") +
  labs(x="Number of MAGs",y=NULL,
       title="A. MAG taxonomic distribution by phylum") +
  pub_theme

# Completeness vs phylum
pS4b <- ggplot(mag_tax,
               aes(x=reorder(Phylum,Completeness,median),
                   y=Completeness,fill=Phylum)) +
  geom_boxplot(width=0.6,color="black",
               outlier.size=1) +
  geom_hline(yintercept=c(50,90),
             linetype="dashed",
             color=c("#F39C12","#27AE60"),
             linewidth=0.5) +
  scale_fill_manual(values=phylum_colors2,guide="none") +
  scale_y_continuous(limits=c(48,103),
                     breaks=seq(50,100,10)) +
  labs(x=NULL,y="Completeness (%)",
       title="B. Completeness distribution by phylum") +
  coord_flip() +
  pub_theme

# Season vs domain kingdom split
domain_season <- mag_tax %>%
  count(Domain,Season,Phylum)

pS4c <- ggplot(mag_tax,
               aes(x=Domain,fill=Season)) +
  geom_bar(position=position_dodge(width=0.8),
           width=0.6,color="white") +
  scale_fill_manual(values=sc,name="Season") +
  scale_y_continuous(expand=expansion(mult=c(0,0.12))) +
  labs(x=NULL,y="Number of MAGs",
       title="C. MAG distribution by domain and season") +
  geom_text(stat="count",
            aes(label=after_stat(count),group=Season),
            position=position_dodge(width=0.8),
            vjust=-0.4,size=3,fontface="bold") +
  pub_theme +
  theme(legend.position="bottom")

# Quality tier summary
qual_data <- mag_tax %>%
  mutate(Quality=ifelse(Completeness>=90,"HQ","MQ")) %>%
  count(Quality,Season,Phylum) %>%
  group_by(Quality,Season) %>%
  summarise(Total=sum(n),.groups="drop")

pS4d <- ggplot(mag_tax %>%
                 mutate(Quality=ifelse(Completeness>=90,
                                       "HQ (>=90%)","MQ (50-90%)")),
               aes(x=Quality,fill=Phylum)) +
  geom_bar(width=0.7,color="white",linewidth=0.2) +
  scale_fill_manual(values=phylum_colors2,name="Phylum") +
  scale_y_continuous(expand=expansion(mult=c(0,0.08))) +
  labs(x=NULL,y="Number of MAGs",
       title="D. Phylum composition by quality tier") +
  pub_theme +
  theme(legend.position="right",
        legend.text=element_text(size=7),
        legend.key.size=unit(0.3,"cm"))

figS4 <- (pS4a | pS4b) / (pS4c | pS4d) +
  plot_layout(heights=c(1.2,1))

ggsave(file.path(OUTDIR,"FigureS4_MAG_phylogeny.pdf"),
       figS4,width=170,height=200,units="mm",dpi=300)
ggsave(file.path(OUTDIR,"FigureS4_MAG_phylogeny.png"),
       figS4,width=170,height=200,units="mm",dpi=300)
cat("Supplementary Figure S4 done.\n")

# ============================================================
# SUPP FIGURE S5 — Flagship MAG functional gene categories
# ============================================================
cat("Building Supplementary Figure S5...\n")

# Halomonas_B full gene category breakdown
halo_genes <- data.frame(
  Category=c("Fatty acid\nbiosynthesis",
             "NRPS/PKS\nrelated",
             "Acyltransferase",
             "Lipoprotein\nsecretion (Lol)",
             "ABC\ntransporter",
             "Efflux /\nresistance",
             "Phospholipid\nbiosynthesis",
             "Quorum\nsensing",
             "Two-component\nregulation",
             "Other\nbiosurfactant"),
  Count=c(23,9,14,5,8,2,6,3,4,7),
  MAG="Halomonas_B\n(H2, winter)"
)

mora_genes <- data.frame(
  Category=c("Fatty acid\nbiosynthesis",
             "NRPS/PKS\nrelated",
             "Acyltransferase",
             "Lipoprotein\nsecretion (Lol)",
             "ABC\ntransporter",
             "Efflux /\nresistance",
             "Phospholipid\nbiosynthesis",
             "Quorum\nsensing",
             "Two-component\nregulation",
             "Other\nbiosurfactant"),
  Count=c(14,7,10,5,6,3,4,4,3,5),
  MAG="Moraxella_A\n(H6, winter)"
)

flagship_genes <- rbind(halo_genes,mora_genes)
flagship_genes$MAG <- factor(flagship_genes$MAG,
  levels=c("Halomonas_B\n(H2, winter)",
            "Moraxella_A\n(H6, winter)"))
flagship_genes$Category <- factor(flagship_genes$Category,
  levels=c("Fatty acid\nbiosynthesis",
            "NRPS/PKS\nrelated",
            "Acyltransferase",
            "Lipoprotein\nsecretion (Lol)",
            "ABC\ntransporter",
            "Efflux /\nresistance",
            "Phospholipid\nbiosynthesis",
            "Quorum\nsensing",
            "Two-component\nregulation",
            "Other\nbiosurfactant"))

mag_colors3 <- c(
  "Halomonas_B\n(H2, winter)" = "#3B8BD4",
  "Moraxella_A\n(H6, winter)" = "#9B59B6"
)

pS5a <- ggplot(flagship_genes,
               aes(x=Category,y=Count,fill=MAG)) +
  geom_bar(stat="identity",
           position=position_dodge(width=0.8),
           width=0.7,color="white",linewidth=0.3) +
  geom_text(aes(label=Count,group=MAG),
            position=position_dodge(width=0.8),
            vjust=-0.4,size=2.8,fontface="bold") +
  scale_fill_manual(values=mag_colors3,name=NULL) +
  scale_y_continuous(expand=expansion(mult=c(0,0.2))) +
  labs(x=NULL,y="Number of genes",
       title="A. Biosurfactant-relevant gene categories in flagship MAGs",
       subtitle="Based on Prokka v1.15.6 annotation") +
  pub_theme +
  theme(legend.position="bottom",
        legend.direction="horizontal",
        axis.text.x=element_text(size=7.5,lineheight=0.9))

# Key gene presence/absence comparison
key_genes_df <- data.frame(
  Gene=factor(c("srfAB","wax-dgaT","rhlA","fabG","fabD",
                 "fabH","acpP","lolA","lolB","lolC",
                 "lolD","lolE","lnt","aiiA","farB"),
              levels=rev(c("srfAB","wax-dgaT","rhlA","fabG","fabD",
                            "fabH","acpP","lolA","lolB","lolC",
                            "lolD","lolE","lnt","aiiA","farB"))),
  Function=c("Surfactin synthetase","Wax ester synthase",
             "Rhamnosyltransferase A","3-oxoacyl-ACP reductase",
             "Malonyl-CoA ACP transacylase","3-oxoacyl-ACP synthase III",
             "Acyl carrier protein","Lipoprotein carrier A",
             "Lipoprotein receptor B","Lipoprotein release C",
             "Lipoprotein ATP-binding D","Lipoprotein outer membrane E",
             "Lipoprotein N-acyltransferase","N-acyl-HSL lactonase",
             "Fatty acid resistance"),
  Halomonas_B=c(1,0,0,1,1,1,1,1,1,1,1,1,1,0,0),
  Moraxella_A=c(0,1,0,1,1,1,1,1,1,1,1,1,1,1,1),
  Pathway=c("Surfactin NRPS","Emulsan/wax ester",
            "Rhamnolipid","FA biosynthesis","FA biosynthesis",
            "FA biosynthesis","FA biosynthesis",
            "Lipoprotein secretion","Lipoprotein secretion",
            "Lipoprotein secretion","Lipoprotein secretion",
            "Lipoprotein secretion","Lipoprotein modification",
            "Quorum sensing","Membrane protection")
)

key_long <- key_long <- key_genes_df %>%
  pivot_longer(cols=c(Halomonas_B,Moraxella_A),
               names_to="MAG",values_to="Present") %>%
  mutate(MAG=recode(MAG,
    Halomonas_B="Halomonas_B (H2)",
    Moraxella_A="Moraxella_A (H6)"),
    Fill=ifelse(Present==1,"Present","Absent"))

pS5b <- ggplot(key_long,
               aes(x=MAG,y=Gene,fill=Fill)) +
  geom_tile(color="white",linewidth=0.8,
            width=0.9,height=0.85) +
  geom_text(data=key_long[key_long$Present==1,],
            aes(label="v"),color="white",
            size=3.5,fontface="bold") +
  scale_fill_manual(
    values=c("Present"="#27AE60","Absent"="#ECF0F1"),
    guide="none") +
  scale_x_discrete(position="top") +
  labs(x=NULL,y=NULL,
       title="B. Key biosurfactant gene presence/absence in flagship MAGs") +
  pub_theme +
  theme(axis.text.x=element_text(size=9,face="bold"),
        axis.text.y=element_text(size=8),
        panel.grid=element_blank(),
        axis.ticks=element_blank())

# Genome size and gene count comparison
genome_stats <- data.frame(
  MAG=factor(c("Halomonas_B","Moraxella_A"),
             levels=c("Halomonas_B","Moraxella_A")),
  Total_genes=c(4321,2497),
  Named_genes=c(2734,1539),
  Biosurfactant_genes=c(81,61),
  Genome_Mb=c(4.52,2.80),
  Completeness=c(93.23,99.99)
)

gs_long <- genome_stats %>%
  pivot_longer(cols=c(Total_genes,Named_genes,Biosurfactant_genes),
               names_to="Category",values_to="Count") %>%
  mutate(Category=factor(
    recode(Category,
           Total_genes="Total genes",
           Named_genes="Named genes",
           Biosurfactant_genes="Biosurfactant genes"),
    levels=c("Total genes","Named genes","Biosurfactant genes")))

pS5c <- ggplot(gs_long,
               aes(x=Category,y=Count,fill=MAG)) +
  geom_bar(stat="identity",
           position=position_dodge(width=0.8),
           width=0.7,color="white",linewidth=0.3) +
  geom_text(aes(label=Count,group=MAG),
            position=position_dodge(width=0.8),
            vjust=-0.4,size=3,fontface="bold") +
  scale_fill_manual(
    values=c("Halomonas_B"="#3B8BD4",
             "Moraxella_A"="#9B59B6"),
    name=NULL) +
  scale_y_continuous(expand=expansion(mult=c(0,0.18))) +
  labs(x=NULL,y="Gene count",
       title="C. Genome annotation summary",
       subtitle="Prokka v1.15.6 annotation") +
  pub_theme +
  theme(legend.position="bottom",
        legend.direction="horizontal")

figS5 <- (pS5a | pS5b) / pS5c +
  plot_layout(heights=c(1.3,1))

ggsave(file.path(OUTDIR,"FigureS5_flagship_MAGs.pdf"),
       figS5,width=170,height=200,units="mm",dpi=300)
ggsave(file.path(OUTDIR,"FigureS5_flagship_MAGs.png"),
       figS5,width=170,height=200,units="mm",dpi=300)
cat("Supplementary Figure S5 done.\n")
cat("All supplementary figures complete.\n")
