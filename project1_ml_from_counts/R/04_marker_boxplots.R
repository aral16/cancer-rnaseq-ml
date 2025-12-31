# R/04_marker_boxplots.R
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(dplyr)
})

# Inputs from your pipeline
counts <- read.csv("data/processed/brca_counts.csv", row.names = 1, check.names = FALSE)
meta   <- read.csv("data/processed/brca_metadata.csv", stringsAsFactors = FALSE)

# Ensure sample order matches metadata
counts <- counts[, meta$sample_id, drop = FALSE]

# Markers to plot (skip any missing genes gracefully)
markers <- c("MKI67", "TOP2A", "CDC20", "ESR1", "PGR")
present <- intersect(markers, rownames(counts))
missing <- setdiff(markers, present)
if (length(missing) > 0) message("Missing genes not found in matrix: ", paste(missing, collapse = ", "))

# Simple log transform for visualization (counts can be huge)
df <- counts[present, , drop = FALSE] %>%
  t() %>%
  as.data.frame()
df$sample_id <- rownames(df)

df <- df %>%
  left_join(meta, by = "sample_id") %>%
  pivot_longer(cols = all_of(present), names_to = "gene", values_to = "count") %>%
  mutate(log2_count = log2(count + 1),
         condition = factor(condition, levels = c("normal", "tumor")))

dir.create("results", showWarnings = FALSE)

p <- ggplot(df, aes(x = condition, y = log2_count)) +
  geom_boxplot(outlier.alpha = 0.2) +
  geom_jitter(width = 0.15, alpha = 0.15, size = 0.7) +
  facet_wrap(~ gene, scales = "free_y", ncol = 3) +
  theme_minimal(base_size = 12) +
  labs(
    title = "BRCA marker expression (log2(count+1))",
    x = NULL,
    y = "log2(count + 1)"
  )

ggsave("results/marker_boxplots.png", plot = p, width = 10, height = 6, dpi = 200)
message("Saved: results/marker_boxplots.png")
