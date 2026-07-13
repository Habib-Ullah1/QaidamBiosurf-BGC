suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(scales)
})

PROJECT <- "/data/habib/metagenome"
OUTDIR <- file.path(PROJECT, "biosurfactant/figures")

pub_theme <- theme_classic(base_size=10) +
  theme(
    axis.text=element_text(size=8, color="black"),
    axis.title=element_text(size=9, color="black", face="bold"),
    legend.text=element_text(size=7.5),
    plot.title=element_text(size=9, face="bold", hjust=0),
    plot.subtitle=element_text(size=6.5, color="#555555"),
    panel.grid.major=element_line(color="#EEEEEE", linewidth=0.3),
    legend.key.size=unit(0.35, "cm"),
    plot.margin=margin(5, 15, 5, 5)
  )

cat_levels <- c("FA biosynthesis", "NRPS/PKS related", "Acyltransferase",
                "Lipoprotein secretion", "ABC transporter", "Efflux/resistance",
                "Phospholipid biosyn.", "Quorum sensing",
                "Two-component reg.", "Other biosurfactant")

flagship_genes <- data.frame(
  Category=factor(rep(cat_levels, 2), levels=cat_levels),
  Count=c(23,9,14,5,8,2,6,3,4,7, 14,7,10,5,6,3,4,4,3,5),
  MAG=factor(rep(c("Halomonas_B (H2)","Moraxella_A (H6)"), each=10),
             levels=c("Halomonas_B (H2)","Moraxella_A (H6)"))
)

pS5a <- ggplot(flagship_genes, aes(x=Category, y=Count, fill=MAG)) +
  geom_bar(stat="identity", position=position_dodge(width=0.8),
           width=0.7, color="white", linewidth=0.3) +
  geom_text(aes(label=Count, group=MAG),
            position=position_dodge(width=0.8),
            vjust=-0.4, size=2.8, fontface="bold") +
  scale_fill_manual(values=c("Halomonas_B (H2)"="#3B8BD4",
                              "Moraxella_A (H6)"="#9B59B6"), name=NULL) +
  scale_y_continuous(expand=expansion(mult=c(0, 0.2))) +
  labs(x=NULL, y="Number of genes",
       title="A. Biosurfactant-relevant gene categories",
       subtitle="Prokka v1.15.6 annotation of flagship MAGs") +
  pub_theme +
  theme(legend.position="bottom", legend.direction="horizontal",
        axis.text.x=element_text(size=7, angle=45, hjust=1, vjust=1))

gene_levels <- c("srfAB","wax-dgaT","rhlA","fabG","fabD",
                 "fabH","acpP","lolA","lolB","lolC",
                 "lolD","lolE","lnt","aiiA","farB")

key_long <- data.frame(
  Gene=factor(rep(gene_levels, 2), levels=rev(gene_levels)),
  MAG=factor(rep(c("Halomonas_B","Moraxella_A"), each=15),
             levels=c("Halomonas_B","Moraxella_A")),
  Present=c(1,0,0,1,1,1,1,1,1,1,1,1,1,0,0,
            0,1,0,1,1,1,1,1,1,1,1,1,1,1,1),
  stringsAsFactors=FALSE
)
key_long$Fill <- ifelse(key_long$Present==1, "Present", "Absent")

pS5b <- ggplot(key_long, aes(x=MAG, y=Gene, fill=Fill)) +
  geom_tile(color="white", linewidth=0.8, width=0.9, height=0.85) +
  geom_text(data=key_long[key_long$Present==1,],
            aes(label="v"), color="white", size=3.5, fontface="bold") +
  scale_fill_manual(
    values=c("Present"="#27AE60","Absent"="#ECF0F1"), guide="none") +
  scale_x_discrete(position="top") +
  labs(x=NULL, y=NULL, title="B. Key biosurfactant genes in flagship MAGs") +
  pub_theme +
  theme(axis.text.x=element_text(size=9, face="bold"),
        axis.text.y=element_text(size=8),
        panel.grid=element_blank(), axis.ticks=element_blank(),
        plot.margin=margin(5, 20, 5, 5))

gs_long <- data.frame(
  MAG=factor(rep(c("Halomonas_B","Moraxella_A"), each=3),
             levels=c("Halomonas_B","Moraxella_A")),
  Category=factor(rep(c("Total genes","Named genes","Biosurfactant genes"), 2),
                  levels=c("Total genes","Named genes","Biosurfactant genes")),
  Count=c(4321,2734,81, 2497,1539,61)
)

pS5c <- ggplot(gs_long, aes(x=Category, y=Count, fill=MAG)) +
  geom_bar(stat="identity", position=position_dodge(width=0.8),
           width=0.7, color="white", linewidth=0.3) +
  geom_text(aes(label=Count, group=MAG),
            position=position_dodge(width=0.8),
            vjust=-0.4, size=3, fontface="bold") +
  scale_fill_manual(values=c("Halomonas_B"="#3B8BD4",
                              "Moraxella_A"="#9B59B6"), name=NULL) +
  scale_y_continuous(expand=expansion(mult=c(0, 0.18))) +
  labs(x=NULL, y="Gene count",
       title="C. Genome annotation summary",
       subtitle="Prokka v1.15.6 annotation") +
  pub_theme +
  theme(legend.position="bottom", legend.direction="horizontal")

figS5 <- (pS5a | pS5b) / pS5c + plot_layout(heights=c(1.3, 1))

ggsave(file.path(OUTDIR, "FigureS5_flagship_MAGs.pdf"),
       figS5, width=200, height=210, units="mm", dpi=300)
ggsave(file.path(OUTDIR, "FigureS5_flagship_MAGs.png"),
       figS5, width=200, height=210, units="mm", dpi=300)
cat("S5 done.\n")
