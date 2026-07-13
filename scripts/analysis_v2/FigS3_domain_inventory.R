source("fig_common.R")
DIR <- "/data/habib/metagenome/biosurfactant/analysis_v2_community_abundance"

# condensation subtypes read from nrps_classified.tsv (starterC count per NRPS);
# FAAL/CAL totals from the genome-wide domain inventory (four-method analysis)
nc <- read_tsv(file.path(DIR,"nrps_classified.tsv"), col_names=FALSE, show_col_types=FALSE)
starter <- sum(nc$X3 > 0, na.rm=TRUE)      # NRPS with a lipopeptide starter-C domain
inv <- tibble::tribble(
  ~domain,                         ~count, ~group,
  "FAAL (fatty-acyl-AMP ligase)",   149,   "Fatty-acid loading",
  "CAL (co-enzyme A ligase)",        36,   "Fatty-acid loading",
  "Condensation - LCL",              75,   "Condensation (elongation)",
  "Condensation - DCL",              24,   "Condensation (elongation)",
  "Condensation - Dual",             11,   "Condensation (elongation)",
  "Condensation - Starter",     starter,   "Lipopeptide initiation")
p <- ggplot(inv, aes(count, reorder(domain,count), fill=group)) +
  geom_col(width=.66) +
  geom_text(aes(label=count), hjust=-0.2, size=3) +
  scale_fill_manual(values=c("Fatty-acid loading"="#2a78d6",
     "Condensation (elongation)"="#8a8f98","Lipopeptide initiation"="#e34948"), name=NULL) +
  scale_x_continuous(expand=expansion(mult=c(0,.12))) +
  labs(x="Domain count (all assemblies)", y=NULL,
       subtitle="Fatty-acid loading is abundant; lipopeptide-initiating starter domains are near-absent") +
  theme_pub() + theme(legend.position="top", plot.subtitle=element_text(size=9))
ggsave(file.path(DIR,"FigS3_domain_inventory.pdf"), p, width=8, height=4.5)
ggsave(file.path(DIR,"FigS3_domain_inventory.png"), p, width=8, height=4.5, dpi=300)
message("wrote FigS3; starter-C NRPS = ", starter)
