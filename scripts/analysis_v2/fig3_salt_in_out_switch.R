source("fig_common.R")
DIR <- "/data/habib/metagenome/biosurfactant/analysis_v2_community_abundance"
CTRL <- "/data/habib/metagenome/biosurfactant/integration/wp1_compositional"
m  <- read_tsv(file.path(DIR,"mag_master_table.tsv"), show_col_types=FALSE)
ac <- read_tsv(file.path(DIR,"proteome_acidity.tsv"), show_col_types=FALSE)

# ---------- Panel A: proteome acidity (UNCHANGED) ----------
eps <- 1e-3
d3 <- m %>% left_join(ac, by="MAG") %>% filter(!is.na(acidic_minus_basic)) %>%
  mutate(ratio=pmax(pmin(log2((winterAbund+eps)/(summerAbund+eps)),6),-6),
         genus_disp=clean_genus(genus),
         grp=case_when(phylum %in% c("Halobacteriota","Nanohaloarchaeota") ~ "Archaea (salt-in)",
                       genus=="Salinibacter" ~ "Salinibacter (salt-in bacterium)",
                       TRUE ~ "Bacteria (salt-out)"))
lab3 <- c("Haloarcula","Halobaculum","Haloferax","Salinibacter","Halomonas",
          "Loktanella","Moraxella","Persicimonas","Roseovarius","Nanosalina")
p3a <- ggplot(d3, aes(acidic_minus_basic, ratio)) +
  annotate("rect", xmin=5.5, xmax=Inf, ymin=-Inf, ymax=Inf, fill="#eda100", alpha=.06) +
  annotate("text", x=7.6, y=6.2, label="salt-in zone", color="#b8860b", size=3.2)+
  geom_vline(xintercept=5.5, linetype=2, color="#c9930f") +
  geom_hline(yintercept=0, linetype=2, color="grey60") +
  geom_point(aes(fill=grp, size=nBGC), shape=21, color="grey30", alpha=.9) +
  geom_text_repel(data=subset(d3, genus_disp %in% lab3), aes(label=genus_disp),
                  size=2.9, max.overlaps=Inf, min.segment.length=0, box.padding=.5,
                  segment.size=.25, segment.color="grey55") +
  scale_fill_manual(values=c("Archaea (salt-in)"="#eda100",
      "Salinibacter (salt-in bacterium)"="#1baf7a","Bacteria (salt-out)"="#2a78d6"),name=NULL)+
  scale_size_continuous(range=c(1.5,6), guide="none") +
  labs(x="Proteome acidity (% acidic \u2212 % basic residues)  \u2192 salt-in",
       y=expression(log[2]*"(winter / summer)")) + theme_pub() + theme(legend.position="top")

# ---------- Panel B: bacterial-compartment-controlled functions (v5) ----------
ft <- read_tsv(file.path(CTRL,"function_table_controlled.tsv"), show_col_types=FALSE) %>%
      rename(func=`function`)
L2 <- function(name, basis="bact"){
  r <- ft[ft$func==name,]; if(nrow(r)==0) return(NA_real_)
  log2(if(basis=="bact") r$WS_bact else r$WS_whole)
}
d3b <- tribble(
  ~lab, ~func, ~basis,
  "Ectoine (biosynthesis)","ectoine","bact",
  "Trehalose","trehalose","bact",
  "Glycine betaine","glycine_betaine","bact",
  "Rhodopsin *","rhodopsin","whole") %>%
  rowwise() %>% mutate(l2=L2(func,basis),
                       dir=ifelse(l2>0,"winter (salt-out)","summer (salt-in)")) %>% ungroup()
stopifnot(nrow(d3b)>0)
p3b <- ggplot(d3b, aes(l2, reorder(lab,l2), color=dir)) +
  geom_vline(xintercept=0, color="grey60") +
  geom_segment(aes(x=0,xend=l2,yend=lab), linewidth=1) + geom_point(size=4) +
  scale_color_manual(values=c("winter (salt-out)"="#2a78d6","summer (salt-in)"="#eda100"),name=NULL)+
  labs(x=expression(log[2]*"(winter / summer)"), y=NULL,
       caption="Solutes: bacterial compartment; rhodopsin: whole community; K-uptake omitted.\n* rhodopsin from single-KO proxy (K04641); see Methods") +
  theme_pub() + theme(legend.position="top", plot.caption=element_text(size=8,hjust=0))

# ---------- merge A | B ----------
fig3 <- (p3a | p3b) + plot_layout(widths=c(1.3,1)) + plot_annotation(tag_levels="A")
ggsave(file.path(DIR,"Fig3_saltin_saltout_switch.pdf"), fig3, width=11, height=5.5)
ggsave(file.path(DIR,"Fig3_saltin_saltout_switch.png"), fig3, width=11, height=5.5, dpi=300)
message("wrote merged Fig3 (A|B)")
