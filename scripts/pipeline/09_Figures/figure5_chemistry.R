#!/usr/bin/env Rscript
# Figure 5 — chemistry: novelty (with fragmentation caveat), C-domain subtypes, flagship architectures
suppressWarnings(suppressMessages({
  need <- c("ggplot2","patchwork","dplyr","scales","grid")
  miss <- need[!need %in% rownames(installed.packages())]
  if (length(miss)) install.packages(miss, repos="https://cloud.r-project.org")
  library(ggplot2); library(patchwork); library(dplyr); library(scales); library(grid)
}))
B <- "/data/habib/metagenome/biosurfactant"
geq <- "\u2265"
th <- theme_classic(base_size=8, base_family="sans") +
  theme(axis.text=element_text(colour="#444444", size=7.3), axis.title=element_text(colour="#111111", size=8.3),
        axis.line=element_line(colour="#444444", linewidth=0.4), axis.ticks=element_line(colour="#444444", linewidth=0.4),
        plot.tag=element_text(face="bold", size=12),
        plot.title=element_text(size=8.3, face="bold", colour="#111111"),
        plot.subtitle=element_text(size=6.5, colour="#666666"),
        legend.position="none", plot.margin=margin(6,8,4,6))

# a: novelty
nov <- data.frame(cat=factor(c("Novel","MIBiG"),levels=c("MIBiG","Novel")), n=c(9,1459))
pa <- ggplot(nov, aes(1, n, fill=cat)) +
  geom_col(width=0.55, colour="white", linewidth=0.4) +
  annotate("text", x=1.0, y=780, label="99.4%", colour="white", fontface="bold", size=3.9) +
  annotate("text", x=1, y=650, label="novel", colour="white", size=2.8) +
  annotate("text", x=1.55, y=1250, label="9 BGCs\nMIBiG-similar\n(0.6%)", size=2.4, colour="#C0392B", hjust=0, lineheight=1.0) +
  annotate("segment", x=1.30, xend=1.52, y=1455, yend=1330, colour="#C0392B", linewidth=0.4) +
  scale_fill_manual(values=c("MIBiG"="#C0392B","Novel"="#2166AC")) +
  scale_x_continuous(limits=c(0.6,2.3)) + scale_y_continuous(expand=expansion(mult=c(0,0.05))) +
  labs(title="BGC novelty vs MIBiG", subtitle="1,468 BGCs, distance < 0.50",
       x=NULL, y="Biosynthetic gene clusters", tag="a") +
  th + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank())

# b: fragmentation caveat
cav <- data.frame(bin=factor(c("<5 kb",paste0(geq,"5 kb")),levels=c("<5 kb",paste0(geq,"5 kb"))), n=c(1186,282))
pb <- ggplot(cav, aes(bin,n,fill=bin)) +
  geom_col(width=0.6, colour="white", linewidth=0.4) +
  geom_text(aes(label=paste0(round(100*n/1468),"%")), vjust=-0.4, size=2.7, colour="#333333") +
  annotate("text", x=1.5, y=1130, label="median 1,977 bp", size=2.5, colour="#B06A00", fontface="italic") +
  scale_fill_manual(values=setNames(c("#E69F00","#999999"), c("<5 kb",paste0(geq,"5 kb")))) +
  scale_y_continuous(expand=expansion(mult=c(0,0.22))) +
  labs(title="Fragmentation caveat", x=NULL, y="BGCs", tag="b") + th

# c: C-domain subtypes (VERIFIED counts from antiSMASH: LCL 150, Epi 60, DCL 48, Dual 22, Starter 12)
cd <- data.frame(dom=c("LCL","Epimerization","DCL","Dual","Starter"), n=c(150,60,48,22,12))
cd$dom <- factor(cd$dom, levels=cd$dom[order(cd$n)])
cd$grp <- ifelse(cd$dom=="LCL","Generic L-peptide","D-amino-acid / lipopeptide")
pc <- ggplot(cd, aes(n, dom, fill=grp)) +
  geom_col(width=0.7, colour="white", linewidth=0.4) +
  geom_text(aes(label=n), hjust=-0.35, size=2.7, colour="#333333") +
  scale_fill_manual(values=c("Generic L-peptide"="#BFBFBF","D-amino-acid / lipopeptide"="#009E73"), name=NULL) +
  scale_x_continuous(limits=c(0,178), expand=expansion(mult=c(0,0))) +
  labs(title="C-domain subtypes", x="Domain count", y=NULL, tag="c") + guides(fill=guide_legend(nrow=2)) +
  th + theme(legend.position="bottom", legend.text=element_text(size=5.8),
             legend.key.size=unit(0.24,"cm"), legend.margin=margin(0,0,0,0), legend.box.spacing=unit(2,"pt"))

# d: flagship domain architectures
dom <- rbind(
  data.frame(y=2, start=c(20,360,770,910), end=c(340,750,880,1020),
             lab=c("C","A","A_C","PCP"), pf=c("PF00668","PF00501","PF13193","PF00550"),
             col=c("#0072B2","#D55E00","#E69F00","#009E73")),
  data.frame(y=1, start=c(20,250), end=c(240,440),
             lab=c("WES_acyltransf","WS_DGAT_C"), pf=c("PF03007","PF06974"),
             col=c("#CC79A7","#56B4E9")))
backbone <- data.frame(y=c(1,2), L=c(452,1057))
titles <- data.frame(y=c(2.62,1.48),
  lab=c("Halomonas_B \u2014 surfactin-type NRPS module (1,057 aa)",
        "Moraxella_A \u2014 wax-ester synthase / WSD (452 aa)"))
pd <- ggplot() +
  geom_segment(data=backbone, aes(x=0, xend=L, y=y, yend=y), colour="#CCCCCC", linewidth=0.7) +
  geom_rect(data=dom, aes(xmin=start, xmax=end, ymin=y-0.20, ymax=y+0.20, fill=I(col)), colour="white", linewidth=0.5) +
  geom_text(data=dom, aes(x=(start+end)/2, y=y, label=lab), size=2.3, colour="white", fontface="bold") +
  geom_text(data=dom, aes(x=(start+end)/2, y=y-0.30, label=pf), size=2.0, colour="#666666") +
  geom_text(data=titles, aes(x=0, y=y, label=lab), hjust=0, size=2.7, fontface="italic", colour="#222222") +
  scale_y_continuous(limits=c(0.55,2.8)) + scale_x_continuous(limits=c(-15,1090)) +
  labs(title="Domain-verified flagship biosurfactant enzymes (hmmscan vs Pfam)",
       x="Protein position (aa)", y=NULL, tag="d") +
  th + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank())

layout <- "
ABC
DDD
"
fig <- pa + pb + pc + pd + plot_layout(design=layout, heights=c(1, 0.75))
dir.create(file.path(B,"figures"), showWarnings=FALSE)
ggsave(file.path(B,"figures/Figure5_chemistry.pdf"), fig, width=185, height=130, units="mm", device=cairo_pdf)
ggsave(file.path(B,"figures/Figure5_chemistry.png"), fig, width=185, height=130, units="mm", dpi=400)
cat("written: figures/Figure5_chemistry.pdf (+ .png)\n")

