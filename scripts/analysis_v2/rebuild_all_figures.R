#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(tidyverse); library(patchwork)})
BASE <- "/data/habib/metagenome/biosurfactant"
DIR  <- file.path(BASE,"integration/wp1_compositional")
AN   <- file.path(BASE,"analysis_v2_community_abundance")
theme_pub <- function(b=12) theme_classic(b)+theme(axis.text=element_text(color="black"))
SEASON <- c(summer="#eda100", winter="#2a78d6")

ft <- read_tsv(file.path(DIR,"function_table_controlled.tsv"), show_col_types=FALSE) |>
      rename(func=`function`)
af <- read_tsv(file.path(DIR,"archaeal_fraction.tsv"), show_col_types=FALSE)
ca <- read_tsv(file.path(DIR,"carrier_table.tsv"), show_col_types=FALSE) |>
      rename(fn=`function`)
cat("SANITY ectoine WS_bact =", ft$WS_bact[ft$func=="ectoine"],
    " (expect ~1.5 after the doe fix)\n")

# log2 ratio for a function on a chosen basis (whole or bacterial)
L2 <- function(name, basis="bact"){
  r <- ft[ft$func==name,]; if(nrow(r)==0) return(NA_real_)
  log2(if(basis=="bact") r$WS_bact else r$WS_whole)
}

## ---------------- Fig 3B: osmoadaptation functional panel ----------------
f3 <- tibble(
  lab  = c("Ectoine (biosynthesis)","Trehalose","Glycine betaine","Rhodopsin"),
  func = c("ectoine","trehalose","glycine_betaine","rhodopsin"),
  guild= c("Salt-out (winter)","Salt-out (winter)","Salt-out (winter)","Salt-in (summer)"),
  basis= c("bact","bact","bact","whole")) |>
  rowwise() |> mutate(l2=L2(func,basis)) |> ungroup()
p3b <- ggplot(f3, aes(l2, reorder(lab,l2), fill=guild)) +
  geom_vline(xintercept=0,color="grey60")+geom_col(width=.6)+
  scale_fill_manual(values=c("Salt-out (winter)"="#2a78d6","Salt-in (summer)"="#eda100"),name=NULL)+
  labs(x=expression(log[2]*"(winter / summer)"),y=NULL,
       subtitle="Osmoadaptation (solutes: bacterial compartment; rhodopsin: whole community; K-uptake omitted)")+
  theme_pub()+theme(legend.position="top",plot.subtitle=element_text(size=9))
ggsave(file.path(AN,"Fig3B_osmoadaptation_panel.pdf"),p3b,width=6.5,height=3.4)
ggsave(file.path(AN,"Fig3B_osmoadaptation_panel.png"),p3b,width=6.5,height=3.4,dpi=300)

## ---------------- Fig 4: surface-active response (A|B / C) ----------------
labA <- c(phospholipaseA="Phospholipase A", lysophospholipase="Lysophospholipase*",
  FA_desaturase="Fatty-acid desaturase", EPS_wza="Capsule export (wza)",
  EPS_wzc="Capsule export (wzc)", EPS_wzb="Capsule export (wzb)",
  EPS_kps_group2="Group-2 capsule (kps)*", EPS_wzxwzy_ctrl="O-antigen (control)")
d4a <- ft |> filter(func %in% names(labA)) |>
  mutate(lab=labA[func], l2=log2(WS_bact),
         grp=case_when(tier=="control"~"O-antigen control",
                       class=="Membrane"~"Membrane remodeling", TRUE~"Secreted EPS/capsule"),
         a=ifelse(tier=="supporting",.55,1))
p4a <- ggplot(d4a, aes(l2, reorder(lab,l2), fill=grp, alpha=I(a))) +
  geom_vline(xintercept=0,color="grey60")+geom_col(width=.62)+
  scale_fill_manual(values=c("Membrane remodeling"="#2a78d6","Secreted EPS/capsule"="#1baf7a",
                             "O-antigen control"="#b0aeb0"),name=NULL)+
  labs(x=expression(log[2]*"(winter / summer, bacterial compartment)"),y=NULL,
       subtitle="Surface-active response (bacterial-compartment controlled; * supporting)")+
  theme_pub()+theme(legend.position="top",plot.subtitle=element_text(size=9))

czp <- file.path(AN,"cazy_family_tpm.tsv")
if(!file.exists(czp)) czp <- file.path(BASE,"integration/cazyme/cazy_family_tpm.tsv")
cz <- read_tsv(czp, show_col_types=FALSE) |>
  mutate(class=str_extract(family,"^(GH|GT|PL|CE|CBM|AA)")) |> filter(!is.na(class)) |>
  group_by(class) |> summarise(summer=sum(summer),winter=sum(winter),.groups="drop") |>
  mutate(l2=log2(winter/summer))
p4b <- ggplot(cz, aes(l2, reorder(class,l2)))+geom_vline(xintercept=0,color="grey60")+
  geom_segment(aes(x=0,xend=l2,yend=class),color="#2a78d6",linewidth=1)+
  geom_point(color="#2a78d6",size=4)+
  labs(x=expression(log[2]*"(winter / summer TPM)"),y="CAZy class",
       subtitle="Carbohydrate-active enzymes")+theme_pub()

grid <- expand_grid(genus=unique(ca$genus), fn=unique(ca$fn)) |>
  left_join(ca, by=c("genus","fn")) |> mutate(genus=str_replace(genus,"_[A-Z]$",""))
p4c <- ggplot(grid, aes(fn, genus, fill=winter_TPM))+geom_tile(color="white")+
  geom_text(aes(label=ifelse(is.na(winter_TPM),"n.d.",sprintf("%.1f",winter_TPM))),size=3)+
  scale_fill_gradient(low="#eef4fb",high="#2a78d6",na.value="grey92",name="Winter TPM")+
  labs(x=NULL,y=NULL,subtitle="Verified carriers (recomputed winter TPM)")+
  theme_pub()+theme(axis.text.x=element_text(angle=25,hjust=1))
fig4 <- (p4a|p4b)/p4c + plot_layout(heights=c(1,.85)) + plot_annotation(tag_levels="A")
ggsave(file.path(AN,"Fig4_amphiphile_response.pdf"),fig4,width=11,height=8)
ggsave(file.path(AN,"Fig4_amphiphile_response.png"),fig4,width=11,height=8,dpi=300)

## ---------------- Fig 6: energy & element cycling ----------------
f6 <- tibble(
  lab  = c("Rhodopsin","AAP (bacteriochl.)","N fixation","Denitrification","Sulfur oxidation"),
  func = c("rhodopsin","AAP_bacteriochl","N_fixation","denitrification","sulfur_ox_sox"),
  cat  = c("Phototrophy","Phototrophy","Nitrogen","Nitrogen","Sulfur"),
  basis= c("whole","bact","whole","whole","bact")) |>
  rowwise() |> mutate(l2=L2(func,basis)) |> ungroup() |>
  mutate(season=ifelse(l2>0,"winter","summer"), cat=factor(cat,levels=c("Phototrophy","Nitrogen","Sulfur")))
p6 <- ggplot(f6, aes(l2, reorder(lab,l2), fill=season))+
  geom_vline(xintercept=0,color="grey60")+geom_col(width=.6)+
  facet_grid(cat~., scales="free_y", space="free_y")+
  scale_fill_manual(values=SEASON,name=NULL)+
  labs(x=expression(log[2]*"(winter / summer)"),y=NULL,
       subtitle="Energy & element cycling (AAP, sulfur oxidation: bacterial compartment)")+
  theme_pub()+theme(legend.position="top",strip.text.y=element_text(angle=0),plot.subtitle=element_text(size=9))
ggsave(file.path(AN,"Fig6_energy_element_cycling.pdf"),p6,width=6.5,height=5)
ggsave(file.path(AN,"Fig6_energy_element_cycling.png"),p6,width=6.5,height=5,dpi=300)

## ---------------- Table S3 ----------------
tS3 <- ft |> transmute(Function=func, Class=class, Tier=tier,
  `Winter/Summer (whole)`=round(WS_whole,2), `Winter/Summer (bacterial)`=round(WS_bact,2))
write_tsv(tS3, file.path(DIR,"TableS3_whole_vs_bacterial.tsv"))

## ---------------- Fig S1: archaeal read fraction ----------------
pS1 <- ggplot(af, aes(sample, pct_archaea, fill=season))+geom_col(width=.6)+
  scale_fill_manual(values=SEASON)+
  labs(x=NULL,y="Archaeal share of mapped reads (%)",
       subtitle="Seasonal archaeal collapse (summer ~58% -> winter ~21%)")+theme_pub()
ggsave(file.path(DIR,"FigS1_archaeal_fraction.pdf"),pS1,width=5,height=4)
ggsave(file.path(DIR,"FigS1_archaeal_fraction.png"),pS1,width=5,height=4,dpi=300)

message("DONE: Fig3B, Fig4, Fig6 (in ",AN,"); TableS3, FigS1 (in ",DIR,")")
