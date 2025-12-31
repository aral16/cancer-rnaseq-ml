import pandas as pd

# Inputs
counts_path = "data/processed/brca_counts.csv"         # genes x samples
meta_path   = "data/processed/brca_metadata.csv"       # sample_id, condition
de_path     = "results/deseq2_all_genes.csv"           # gene, padj, log2FoldChange

# Output
out_path = "data/processed/ml_top200_features.csv"

# Load
counts = pd.read_csv(counts_path, index_col=0)   # genes x samples
meta   = pd.read_csv(meta_path)
de     = pd.read_csv(de_path)

# Choose top 200 genes by adjusted p-value (exclude NAs)
de = de.dropna(subset=["padj"])
top200 = de.sort_values("padj").head(200)["gene"].tolist()

# Ensure genes exist in counts matrix
top200_present = [g for g in top200 if g in counts.index]
if len(top200_present) < 150:
    raise ValueError(f"Only {len(top200_present)} of top200 genes found in counts. Check gene IDs.")

# Build sample x gene matrix
X = counts.loc[top200_present].T  # samples x genes

# Ensure sample order matches metadata
X = X.loc[meta["sample_id"]]

# Save
X.to_csv(out_path)
print(f"Saved ML feature matrix: {out_path} with shape {X.shape}")
print(f"Missing genes from top200: {len(top200) - len(top200_present)}")
