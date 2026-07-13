source("fig_common.R")
DIR <- "/data/habib/metagenome/biosurfactant/analysis_v2_community_abundance"
m <- read_tsv(file.path(DIR,"mag_master_table.tsv"), show_col_types=FALSE)
e <- read_tsv(file.path(DIR,"genome_economics.tsv"), show_col_types=FALSE)
eps <- 1e-3

m <- m %>% mutate(
  phylum = ifelse(is.na(phylum)|phylum=="","unclassified",phylum),
  genus_disp = clean_genus(ifelse(!is.na(genus)&genus!="", genus, MAG)),
  ratio = pmax(pmin(log2((winterAbund+eps)/(summerAbund+eps)),6),-6),
  rowlab = paste0(genus_disp," (",substr(MAG,1,2),")"))

p2a <- ggplot(m, aes(ratio, reorder(rowlab, ratio))) +
  geom_vline(xintercept=0, linetype=2, color="grey60") +
  geom_point(aes(size=nBGC, fill=phylum), shape=21, color="grey30", alpha=.9) +
  annotate("text", x=-3.5, y=nrow(m)+0.8, label="summer-favoring", size=3.1)+
  annotate("text", x= 3.5, y=nrow(m)+0.8, label="winter-favoring", size=3.1)+
  scale_fill_manual(values=PAL, name="Phylum") +
  scale_size_continuous(range=c(1.5,7), name="BGCs") +
  coord_cartesian(clip="off", ylim=c(1,nrow(m)+1.5)) +
  labs(x=expression(log[2]*"(winter / summer abundance)"), y=NULL) +
  theme_pub(11) + theme(axis.text.y=element_text(size=7.2, color="black"),
                        plot.margin=margin(18,10,10,10))

sz <- intersect(c("gsize_mb","est_gsize_mb","genome_size_mb","size_mb"), colnames(e))[1]
e2 <- e %>% mutate(gsize=.data[[sz]], phylum=ifelse(is.na(phylum)|phylum=="","unclassified",phylum))
p2b <- ggplot(e2, aes(gsize, nBGC, fill=phylum)) +
  geom_point(shape=21, size=3, color="grey30", alpha=.9) +
  scale_fill_manual(values=PAL, guide="none") +
  labs(x="Estimated genome size (Mb)", y="BGCs per genome") + theme_pub(11)

fig2 <- (p2a | p2b) + plot_layout(widths=c(1.5,1)) + plot_annotation(tag_levels="A")
ggsave(file.path(DIR,"Fig2_genome_landscape.pdf"), fig2, width=11.5, height=7.5)
ggsave(file.path(DIR,"Fig2_genome_landscape.png"), fig2, width=11.5, height=7.5, dpi=300)
message("wrote Fig2")
