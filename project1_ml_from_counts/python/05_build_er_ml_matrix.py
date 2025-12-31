import pandas as pd

counts = pd.read_csv("data/processed/brca_counts.csv", index_col=0)
meta   = pd.read_csv("data/processed/brca_metadata_er.csv")
de     = pd.read_csv("results/deseq2_all_genes.csv")

# Feature selection
de = de.dropna(subset=["padj"])
top_genes = de.sort_values("padj").head(300)["gene"].tolist()

# Remove direct hormone receptors (no leakage)
leakage = {"ESR1", "ESR2", "PGR"}
top_genes = [g for g in top_genes if g not in leakage and g in counts.index]

# Build matrix
X = counts.loc[top_genes].T
X = X.loc[meta["sample_id"]]

X.to_csv("data/processed/ml_er_features.csv")

print(f"Feature matrix shape: {X.shape}")
print("Saved: data/processed/ml_er_features.csv")
