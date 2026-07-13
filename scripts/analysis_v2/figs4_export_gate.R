source("fig_common.R")
DIR <- "/data/habib/metagenome/biosurfactant/analysis_v2_community_abundance"

eps <- tibble::tribble(
  ~gene, ~summer, ~winter, ~system,
  "wza",   7.3,  56.0, "Capsule export (Wza\u2013Wzc)",
  "wzc",  21.5, 132.5, "Capsule export (Wza\u2013Wzc)",
  "wzb",   2.0,   7.8, "Capsule export (Wza\u2013Wzc)",
  "kpsE",  2.9,  42.8, "Group-2 capsule (Kps)",
  "kpsT",  4.0,  29.5, "Group-2 capsule (Kps)",
  "kpsD", 13.7,  55.0, "Group-2 capsule (Kps)",
  "kpsM",  2.4,  13.3, "Group-2 capsule (Kps)",
  "exoP",  1.2,   2.7, "Group-2 capsule (Kps)",
  "wzx",   6.0,   4.1, "Wzy-dependent O-antigen (control)",
  "wzy",   0.6,   0.0, "Wzy-dependent O-antigen (control)",
  "wzm",   6.6,   2.8, "Wzy-dependent O-antigen (control)",
  "wzt",   2.9,   2.6, "Wzy-dependent O-antigen (control)") |>
  dplyr::mutate(l2 = log2((winter+0.1)/(summer+0.1)))

p <- ggplot(eps, aes(l2, reorder(gene, l2), fill=system)) +
  geom_vline(xintercept=0, color="grey60") +
  geom_col(width=.65) +
  scale_fill_manual(values=c(
    "Capsule export (Wza\u2013Wzc)"="#2a78d6",
    "Group-2 capsule (Kps)"="#1baf7a",
    "Wzy-dependent O-antigen (control)"="#b0aeb0"), name=NULL) +
  labs(x=expression(log[2]*"(winter / summer TPM)"), y=NULL,
       subtitle="Capsule export coordinately winter-enriched; Wzy-dependent route flat (specificity control)") +
  theme_pub() + theme(legend.position="top", plot.subtitle=element_text(size=9))

ggsave(file.path(DIR,"FigS4_EPS_export_gate.pdf"), p, width=8, height=5.5)
ggsave(file.path(DIR,"FigS4_EPS_export_gate.png"), p, width=8, height=5.5, dpi=300)
message("wrote FigS4")
