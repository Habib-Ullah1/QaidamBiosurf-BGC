suppressPackageStartupMessages({library(ggtree); library(ggplot2); library(treeio); library(dplyr); library(stringr)})
DIR <- "/data/habib/metagenome/biosurfactant/analysis_v2_community_abundance"
tr  <- read.tree(file.path(DIR,"figS2","figS2_Cdomain.treefile"))

grp <- function(x) case_when(
  str_detect(x, regex("Srf|Itu|Fen|Rhl|surfactin|iturin|fengycin|rhamnolipid", ignore_case=TRUE)) ~ "Biosurfactant reference",
  str_detect(x, "^(C3|C6|H2|H6)\\|") ~ "Query (this study)",
  TRUE ~ "MIBiG reference")
meta <- data.frame(label=tr$tip.label, group=grp(tr$tip.label))
message("tip composition:"); print(table(meta$group))

p <- ggtree(tr, size=0.3) %<+% meta +
  geom_tippoint(aes(color=group), size=1.5) +
  geom_tiplab(aes(subset=(group=="Biosurfactant reference"), label=label),
              size=2.6, color="#e34948", offset=0.03) +
  scale_color_manual(values=c("Biosurfactant reference"="#e34948",
                              "Query (this study)"="#2a78d6",
                              "MIBiG reference"="#b0aeb0"), name=NULL) +
  theme_tree2() + theme(legend.position="top")

ggsave(file.path(DIR,"FigS2_Cdomain_phylogeny.pdf"), p, width=8, height=11)
ggsave(file.path(DIR,"FigS2_Cdomain_phylogeny.png"), p, width=8, height=11, dpi=300)
message("wrote FigS2")
