#!/usr/bin/env Rscript
# WP1 part2 — ALDEx2 CLR effect sizes; whole vs bacterial-only.
suppressMessages(library(ALDEx2))
BASE<-"/data/habib/metagenome/biosurfactant"; OUT<-file.path(BASE,"integration/wp1_compositional")
S<-c("C3","C6","H2","H6"); G<-c("summer","summer","winter","winter")  # winter>summer alphabetically => +effect = winter
FOCUS<-list(ectoine=c("K06718","K06720","K00836"),glycine_betaine=c("K00108","K00130"),
 trehalose=c("K00697","K01087"),K_uptake=c("K03498","K03499","K01546","K01547","K01548"),
 rhodopsin="K04641",wza="K01991",phospholipaseA=c("K01058","K01059"),
 lysophospholipase=c("K01048","K14676","K06999"))
run<-function(path,tag){
  m<-read.delim(path,row.names=1,check.names=FALSE)[,S]; m<-m[rowSums(m)>=4,]
  cat("\n====",tag,"====  features:",nrow(m),"\n"); set.seed(1)
  x<-aldex.clr(round(m),conds=G,mc.samples=128,denom="all",verbose=FALSE)
  ef<-aldex.effect(x); tt<-aldex.ttest(x)
  r<-data.frame(KO=rownames(ef),diff=ef$diff.btw,effect=ef$effect,we.ep=tt$we.ep)
  r<-r[order(-r$effect),]; write.table(r,file.path(OUT,paste0("aldex2_",tag,".tsv")),sep="\t",quote=FALSE,row.names=FALSE)
  cat("top winter (effect>0):\n"); print(head(r,6),row.names=FALSE)
  cat("top summer (effect<0):\n"); print(head(r[order(r$effect),],6),row.names=FALSE)
  r
}
rw<-run(file.path(OUT,"ko_counts_wholecommunity.tsv"),"wholecommunity")
rb<-run(file.path(OUT,"ko_counts_bacterialonly.tsv"),"bacterialonly")
cat("\n==== FOCUS: effect size whole vs bacterial-only (+ = winter) ====\n")
cat(sprintf("%-18s %-9s %8s %8s\n","function","KO","whole","bact"))
for(f in names(FOCUS)) for(ko in FOCUS[[f]]){
  ew<-rw$effect[match(ko,rw$KO)]; eb<-rb$effect[match(ko,rb$KO)]
  if(!is.na(ew)||!is.na(eb)) cat(sprintf("%-18s %-9s %8s %8s\n",f,ko,
     ifelse(is.na(ew),"-",sprintf("%+.2f",ew)),ifelse(is.na(eb),"-",sprintf("%+.2f",eb))))
}
cat("\nA winter function whose bact effect stays positive = genuine within-bacterial\nenrichment; one that collapses to ~0 was the archaeal-fraction artifact.\n")
