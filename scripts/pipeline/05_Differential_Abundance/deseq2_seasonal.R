library(DESeq2)
library(ggplot2)

PROJECT <- "/data/habib/metagenome"
OUTDIR  <- file.path(PROJECT, "biosurfactant/14_deseq2")
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

cat("Loading CoverM coverage data...\n")

samples <- c("C3","C6","H2","H6")
cov_list <- list()

for(s in samples){
  f <- file.path(PROJECT,"biosurfactant",s,"13_coverage/coverM",
                 paste0(s,"_contig_coverage.tsv"))
  df <- read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  colnames(df)[1] <- "contig"
  cov_col <- grep("Trimmed.Mean|trimmed_mean", colnames(df), value=TRUE)[1]
  cov_list[[s]] <- df[,c("contig",cov_col)]
  colnames(cov_list[[s]])[2] <- s
  cat(sprintf("  %s: %d contigs\n", s, nrow(cov_list[[s]])))
}

cat("Merging coverage matrix...\n")
cov_mat <- Reduce(function(a,b) merge(a,b,by="contig",all=TRUE), cov_list)
rownames(cov_mat) <- cov_mat$contig
cov_mat$contig <- NULL
cov_mat[is.na(cov_mat)] <- 0
cat(sprintf("Coverage matrix: %d contigs x %d samples\n",
            nrow(cov_mat), ncol(cov_mat)))

cat("Loading biosurfactant gene contigs...\n")
bs_contigs <- c()
for(s in samples){
  f <- file.path(PROJECT,"biosurfactant",s,"03_hmmer",
                 paste0(s,"_biosurfactant_hmmer_hits.tsv"))
  df <- read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  contigs <- sub("_[0-9]+$","", df$orf_id)
  bs_contigs <- union(bs_contigs, contigs)
}
cat(sprintf("  Unique biosurfactant contigs: %d\n", length(bs_contigs)))

bs_mat <- cov_mat[rownames(cov_mat) %in% bs_contigs,]
cat(sprintf("  In coverage matrix: %d\n", nrow(bs_mat)))

count_mat <- round(bs_mat * 1000)
count_mat <- count_mat[rowSums(count_mat) > 0,]
cat(sprintf("  Non-zero contigs: %d\n", nrow(count_mat)))

col_data <- data.frame(
  sample = samples,
  season = c("summer","summer","winter","winter"),
  row.names = samples
)
col_data$season <- factor(col_data$season, levels=c("summer","winter"))

cat("Running DESeq2...\n")
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = col_data,
  design    = ~ season
)
dds <- dds[rowSums(counts(dds) >= 10) >= 2,]
cat(sprintf("  Contigs after low-count filter: %d\n", nrow(dds)))
dds <- DESeq(dds)

res <- results(dds,
               contrast=c("season","winter","summer"),
               alpha=0.05)
res_df <- as.data.frame(res)
res_df$contig <- rownames(res_df)
res_df <- res_df[order(res_df$padj, na.last=TRUE),]

write.table(res_df,
            file.path(OUTDIR,"deseq2_winter_vs_summer_all.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

sig <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05,]
sig_summer <- sig[sig$log2FoldChange < 0,]
sig_winter  <- sig[sig$log2FoldChange > 0,]

write.table(sig_summer,
            file.path(OUTDIR,"deseq2_summer_enriched.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(sig_winter,
            file.path(OUTDIR,"deseq2_winter_enriched.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

cat(sprintf("\n=== DESeq2 Results ===\n"))
cat(sprintf("Total contigs tested:    %d\n", nrow(res_df)))
cat(sprintf("Significant (padj<0.05): %d\n", nrow(sig)))
cat(sprintf("  Summer enriched:       %d\n", nrow(sig_summer)))
cat(sprintf("  Winter enriched:       %d\n", nrow(sig_winter)))

# Volcano plot
plot_df <- res_df[!is.na(res_df$padj),]
plot_df$sig <- ifelse(plot_df$padj<0.05 & plot_df$log2FoldChange < -1,
                      "Summer enriched",
               ifelse(plot_df$padj<0.05 & plot_df$log2FoldChange > 1,
                      "Winter enriched","Not significant"))

p <- ggplot(plot_df, aes(x=log2FoldChange, y=-log10(padj), color=sig)) +
  geom_point(alpha=0.6, size=1.5) +
  scale_color_manual(values=c("Summer enriched"="#E8593C",
                               "Winter enriched"="#3B8BD4",
                               "Not significant"="#AAAAAA")) +
  geom_vline(xintercept=c(-1,1), linetype="dashed", alpha=0.5) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", alpha=0.5) +
  labs(title="Biosurfactant gene contigs: Winter vs Summer",
       x="log2 Fold Change (Winter/Summer)",
       y="-log10(adjusted p-value)",
       color="") +
  theme_bw(base_size=12) +
  theme(legend.position="bottom")

ggsave(file.path(OUTDIR,"deseq2_volcano.pdf"), p, width=8, height=6)
ggsave(file.path(OUTDIR,"deseq2_volcano.png"), p, width=8, height=6, dpi=300)

cat("\n=== Top 10 summer-enriched biosurfactant contigs ===\n")
print(head(sig_summer[,c("contig","log2FoldChange","padj","baseMean")],10))
cat("\n=== Top 10 winter-enriched biosurfactant contigs ===\n")
print(head(sig_winter[,c("contig","log2FoldChange","padj","baseMean")],10))

cat(sprintf("\nOutput: %s\n", OUTDIR))
cat("=== DESeq2 complete ===\n")
