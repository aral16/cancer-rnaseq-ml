suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------
# Paths
# -----------------------
tumor_counts_path  <- "data/extracted/GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_FeatureCounts.txt.gz"
normal_counts_path <- "data/extracted/GSM1697009_06_01_15_TCGA_24.normal_Rsubread_FeatureCounts.txt.gz"

tumor_map_path  <- "data/raw/GSE62944_06_01_15_TCGA_24_CancerType_Samples.txt.gz"
normal_map_path <- "data/raw/GSE62944_06_01_15_TCGA_24_Normal_CancerType_Samples.txt.gz"

out_counts_path <- "data/processed/brca_counts.csv"
out_meta_path   <- "data/processed/brca_metadata.csv"

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)

# -----------------------
# 1) Read sample mappings (NO HEADERS)
# -----------------------
# These mapping files are 2-column TSVs:
#   V1 = TCGA sample ID
#   V2 = cancer type (BRCA, LUAD, etc.)
tum_map <- fread(tumor_map_path, header = FALSE)
nor_map <- fread(normal_map_path, header = FALSE)

setnames(tum_map, c("sample_id", "cancer_type"))
setnames(nor_map, c("sample_id", "cancer_type"))

cat("Tumor map head:\n"); print(head(tum_map))
cat("Normal map head:\n"); print(head(nor_map))

brca_tumor_samples  <- tum_map[cancer_type == "BRCA", unique(sample_id)]
brca_normal_samples <- nor_map[cancer_type == "BRCA", unique(sample_id)]

cat("BRCA tumor samples:", length(brca_tumor_samples), "\n")
cat("BRCA normal samples:", length(brca_normal_samples), "\n")

if (length(brca_tumor_samples) == 0 || length(brca_normal_samples) == 0) {
  stop("No BRCA samples found in mapping files. Check cancer_type labels.")
}

# -----------------------
# 2) Read counts matrices
# -----------------------
tum_counts <- fread(tumor_counts_path)
nor_counts <- fread(normal_counts_path)

# first column should be gene id
setnames(tum_counts, names(tum_counts)[1], "gene")
setnames(nor_counts, names(nor_counts)[1], "gene")

cat("Tumor counts ncol:", ncol(tum_counts), "nrow:", nrow(tum_counts), "\n")
cat("Normal counts ncol:", ncol(nor_counts), "nrow:", nrow(nor_counts), "\n")

# -----------------------
# 3) Match sample columns
# -----------------------
select_cols_by_samples <- function(df, sample_ids) {
  cols <- names(df)
  matched <- intersect(sample_ids, cols)

  # fallback: prefix match (rarely needed but safe)
  if (length(matched) < length(sample_ids) * 0.5) {
    matched2 <- c()
    for (s in sample_ids) {
      hit <- cols[startsWith(cols, s)]
      if (length(hit) == 1) matched2 <- c(matched2, hit)
    }
    matched <- unique(c(matched, matched2))
  }
  matched
}

tum_keep <- select_cols_by_samples(tum_counts, brca_tumor_samples)
nor_keep <- select_cols_by_samples(nor_counts, brca_normal_samples)

cat("Matched BRCA tumor columns:", length(tum_keep), "\n")
cat("Matched BRCA normal columns:", length(nor_keep), "\n")

if (length(tum_keep) == 0 || length(nor_keep) == 0) {
  stop("No BRCA samples matched counts matrix columns. Inspect column names in the FeatureCounts files.")
}

tum_sub <- tum_counts[, c("gene", tum_keep), with=FALSE]
nor_sub <- nor_counts[, c("gene", nor_keep), with=FALSE]

# -----------------------
# 4) Merge tumor + normal by gene
# -----------------------
merged <- merge(tum_sub, nor_sub, by="gene", all=TRUE)

# Replace NA with 0
for (j in 2:ncol(merged)) {
  na_idx <- which(is.na(merged[[j]]))
  if (length(na_idx) > 0) set(merged, na_idx, j, 0)
}

# -----------------------
# 5) Write processed outputs
# -----------------------
counts_df <- as.data.frame(merged)
rownames(counts_df) <- counts_df$gene
counts_df$gene <- NULL

write.csv(counts_df, out_counts_path)

all_samples <- colnames(counts_df)
meta <- data.frame(
  sample_id = all_samples,
  condition = ifelse(all_samples %in% tum_keep, "tumor", "normal"),
  stringsAsFactors = FALSE
)
write.csv(meta, out_meta_path, row.names = FALSE)

cat("Wrote:\n", out_counts_path, "\n", out_meta_path, "\n")
