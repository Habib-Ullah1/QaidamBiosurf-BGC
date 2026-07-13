#!/usr/bin/env Rscript
suppressMessages(library(vegan))
BASE<-"/data/habib/metagenome/biosurfactant"; OUT<-file.path(BASE,"integration/wp2_16s_power")
dir.create(OUT,showWarnings=FALSE,recursive=TRUE)
raw<-read.csv(file.path(BASE,"16S_data/level-6.csv"),header=TRUE,check.names=FALSE,row.names=1)
istax<-grepl(";|__",colnames(raw)); tax<-as.matrix(raw[,istax,drop=FALSE]); mode(tax)<-"numeric"
tax[is.na(tax)]<-0; tax<-tax[rowSums(tax)>0,colSums(tax)>0,drop=FALSE]
season<-ifelse(grepl("^C",rownames(tax),ignore.case=TRUE),"summer","winter")
cat("samples:",nrow(tax)," genera:",ncol(tax),"\n"); print(table(season))
rel<-sweep(tax,1,rowSums(tax),"/"); bc<-vegdist(rel,"bray")
set.seed(1); pm<-adonis2(bc~season,permutations=999); cat("\n== PERMANOVA ==\n"); print(pm)
bd<-betadisper(bc,factor(season)); cat("\n== betadisper ==\n"); print(permutest(bd,permutations=999))
sm<-colMeans(rel[season=="summer",,drop=FALSE]); wi<-colMeans(rel[season=="winter",,drop=FALSE])
g<-data.frame(genus=sub(".*g__","",colnames(rel)),full=colnames(rel),
   mean_summer=sm,mean_winter=wi,log2WS=log2((wi+1e-6)/(sm+1e-6)))
g$direction<-ifelse(g$log2WS>0,"winter","summer")
if(requireNamespace("indicspecies",quietly=TRUE)){
  library(indicspecies); set.seed(1)
  ctrl<-if(requireNamespace("permute",quietly=TRUE)) permute::how(nperm=999) else 999
  iv<-multipatt(as.data.frame(rel),season,func="IndVal.g",control=ctrl)
  g$IndVal<-iv$sign$stat[match(g$full,rownames(iv$sign))]
  g$IndVal_p<-iv$sign$p.value[match(g$full,rownames(iv$sign))]
} else cat("\n[!] indicspecies missing\n")
g<-g[order(-abs(g$log2WS)),]
write.table(g,file.path(OUT,"genus_seasonal_fidelity.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
cat("\n== top seasonal genera ==\n"); print(head(g[,c("genus","mean_summer","mean_winter","direction",intersect("IndVal_p",colnames(g)))],15),row.names=FALSE)
carr<-c("Moraxella","Loktanella","Roseovarius","Gillisia","Halomonas","Billgrantia",
 "Franzmannia","Persicimonas","Haloarcula","Halobaculum","Haloferax","Salinibacter")
cat("\n== MG carrier/marker genera in 17-sample amplicon ==\n")
h<-g[g$genus %in% carr | grepl(paste(carr,collapse="|"),g$full),]
print(h[,c("genus","mean_summer","mean_winter","direction",intersect("IndVal_p",colnames(g)))],row.names=FALSE)
write.table(h,file.path(OUT,"carrier_genera_concordance.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
cat("\n[DONE]",OUT,"\n")
