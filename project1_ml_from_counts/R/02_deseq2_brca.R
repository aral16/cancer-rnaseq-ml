suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

counts <- read.csv("data/processed/brca_counts.csv", row.names=1, check.names=FALSE)
meta   <- read.csv("data/processed/brca_metadata.csv", stringsAsFactors = FALSE)

# Ensure sample order matches
counts <- counts[, meta$sample_id]

# DESeq2 requires integer counts
counts_mat <- round(as.matrix(counts))
meta$condition <- factor(meta$condition, levels=c("normal", "tumor"))
rownames(meta) <- meta$sample_id

dds <- DESeqDataSetFromMatrix(countData = counts_mat, colData = meta, design = ~ condition)
dds <- DESeq(dds)

res <- results(dds, contrast=c("condition","tumor","normal"))
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]

dir.create("results", showWarnings = FALSE)

write.csv(res_df, "results/deseq2_all_genes.csv", row.names=FALSE)

sig <- subset(res_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(sig, "results/deseq2_sig_genes.csv", row.names=FALSE)

# Volcano plot
#res_df$neglog10padj <- -log10(res_df$padj)
res_df$neglog10padj <- -log10(pmax(res_df$padj, 1e-300))

res_df$signif <- ifelse(!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "sig", "ns")

p <- ggplot(res_df, aes(x=log2FoldChange, y=neglog10padj)) +
  geom_point(alpha=0.4, size=0.8) +
  geom_point(data=subset(res_df, signif=="sig"), alpha=0.7, size=0.9) +
  theme_minimal() +
  labs(title="BRCA Tumor vs Normal (DESeq2)", x="log2 Fold Change", y="-log10 adj p-value")

ggsave("results/volcano.png", plot=p, width=7, height=5, dpi=200)

# Heatmap: top 50 DE genes by padj
top50 <- head(res_df$gene[!is.na(res_df$padj)], 50)
vsd <- vst(dds, blind=TRUE)
mat <- assay(vsd)[top50, ]

# z-score by gene for visualization
mat_z <- t(scale(t(mat)))

ann <- data.frame(condition=meta$condition)
rownames(ann) <- rownames(meta)

pheatmap(mat_z,
         annotation_col = ann,
         show_colnames = FALSE,
         fontsize_row = 7,
         main = "Top 50 DE genes (VST, z-scored)")
ggsave("results/heatmap_top50.png", width=8, height=10, dpi=200)

cat("DE results written to results/ and plots generated.\n")
cat("Significant genes:", nrow(sig), "\n")
