source("fig_common.R")
DIR <- "/data/habib/metagenome/biosurfactant/analysis_v2_community_abundance"
sf <- read_tsv(file.path(DIR,"seasonal_function_summary.tsv"), show_col_types=FALSE)
message(">>> function rows: ", paste(sf[["function"]], collapse=", "))

want <- tribble(~pat,~lab,~cat,
  "rhodopsin","Rhodopsin *","Phototrophy",
  "AAP|bacteriochloro","AAP (bch)","Phototrophy",
  "N_fixation|nif","N fixation","Nitrogen",
  "denitrification","Denitrification","Nitrogen",
  "nitrification","Nitrification","Nitrogen",
  "sulfur_ox|sox","Sulfur oxidation","Sulfur")
d6 <- want %>% rowwise() %>% mutate(row=list(pick_fn(sf,pat))) %>% unnest(row) %>%
  mutate(l2=log2(winter_TPM/summer_TPM),
         dir=ifelse(l2>0,"winter","summer"))
stopifnot(nrow(d6)>0)
p6 <- ggplot(d6, aes(l2, reorder(lab,l2), color=dir)) +
  geom_vline(xintercept=0, color="grey60") +
  geom_segment(aes(x=0,xend=l2,yend=lab), linewidth=1) +
  geom_point(aes(size=pmax(summer_TPM,winter_TPM))) +
  facet_grid(cat~., scales="free_y", space="free_y", switch="y") +
  scale_color_manual(values=c(winter="#2a78d6", summer="#eda100"), name="enriched in")+
  scale_size_continuous(range=c(2,6), name="max TPM") +
  labs(x=expression(log[2]*"(winter / summer TPM)"), y=NULL,
       caption="* rhodopsin from single-KO proxy (K04641). Denitrification/nitrification near-baseline; see Methods.")+
  theme_pub() + theme(strip.placement="outside", strip.background=element_blank(),
     strip.text.y.left=element_text(angle=0, face="bold"),
     plot.caption=element_text(size=8,hjust=0), legend.position="right")
ggsave(file.path(DIR,"Fig6_element_cycling.pdf"), p6, width=8.5, height=6)
ggsave(file.path(DIR,"Fig6_element_cycling.png"), p6, width=8.5, height=6, dpi=300)
message("wrote Fig6")
