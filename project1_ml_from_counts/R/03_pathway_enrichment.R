suppressPackageStartupMessages({
  library(data.table)
  library(msigdbr)
  library(fgsea)
})

# Load DE results
res <- fread("results/deseq2_all_genes.csv")

# Keep valid rows only
res <- res[!is.na(gene) & !is.na(log2FoldChange) & !is.na(padj)]

# Rank by log2FC (tumor vs normal)
stats <- res$log2FoldChange
names(stats) <- res$gene
stats <- sort(stats, decreasing = TRUE)

# Hallmark gene sets (human)
m_df <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(m_df$gene_symbol, m_df$gs_name)

# Run fgsea
fg <- fgsea(pathways = pathways, stats = stats, nperm = 10000)
fg <- fg[order(padj)]

dir.create("results", showWarnings = FALSE)
fwrite(fg, "results/pathways_hallmark.csv")

# Save top enriched on each side
top_up <- fg[NES > 0][order(padj)][1:10, .(pathway, NES, padj, size)]
top_down <- fg[NES < 0][order(padj)][1:10, .(pathway, NES, padj, size)]

fwrite(top_up, "results/top10_pathways_up.csv")
fwrite(top_down, "results/top10_pathways_down.csv")

print(top_up)
print(top_down)
