###############################################################################
## MetaboAnalyst-Safe Preprocessing + PCA (R)
## - Validates and cleans a feature-by-sample table
## - Harmonizes sample IDs between feature table and metadata
## - Applies zero-handling, log2 transform, and Pareto scaling
## - Runs PCA and generates a publication-ready PC1 vs PC2 plot
## - Exports MetaboAnalyst-ready feature and metadata files
##
## Author: Lee Chua
## Last updated: 2026-01-31
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(tibble)
})

###############################################################################
## 0) Configuration
###############################################################################
# Input files (rename freely; keep them generic and shareable)
input_feature_table <- "metabolomics_feature_table.csv"
input_metadata      <- "sample_metadata.csv"

# Output files
output_feature_table <- "metaboanalyst_PCA_ready_features.csv"
output_metadata      <- "metaboanalyst_PCA_ready_metadata.csv"
output_pca_plot      <- "PCA_PC1_PC2.png"

# Plot settings
ellipse_level <- 0.95
point_size    <- 3
base_fontsize <- 14

###############################################################################
## 1) Helpers
###############################################################################
stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
}

sanitize_ids <- function(x) {
  # Replace non-alphanumeric characters with underscore, then make syntactically valid
  x %>%
    str_replace_all("[^A-Za-z0-9]", "_") %>%
    make.names(unique = TRUE)
}

pareto_scale <- function(x) {
  # Pareto scaling: mean-center and divide by sqrt(sd)
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / sqrt(s)
}

###############################################################################
## 2) Load data
###############################################################################
stop_if_missing(input_feature_table)
stop_if_missing(input_metadata)

feature_raw <- read.csv(
  input_feature_table,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

meta_raw <- read.csv(
  input_metadata,
  stringsAsFactors = FALSE
)

if (ncol(feature_raw) < 2) stop("Feature table must have >= 2 columns (Feature + samples).")
if (!("Sample" %in% colnames(meta_raw))) {
  # Allow first column to be Sample if not explicitly named
  colnames(meta_raw)[1] <- "Sample"
}

###############################################################################
## 3) Remove embedded non-numeric rows from feature table
###############################################################################
# Assumes first column is feature IDs and remaining columns should be numeric.
numeric_rows <- apply(feature_raw[, -1, drop = FALSE], 1, function(x) {
  all(!is.na(suppressWarnings(as.numeric(x))))
})

feature_clean <- feature_raw[numeric_rows, , drop = FALSE]

###############################################################################
## 4) Standardize feature column
###############################################################################
colnames(feature_clean)[1] <- "Feature"
feature_clean$Feature <- make.names(feature_clean$Feature, unique = TRUE)

###############################################################################
## 5) Sanitize sample names in feature table (CRITICAL)
###############################################################################
sample_names_clean <- sanitize_ids(colnames(feature_clean)[-1])
colnames(feature_clean)[-1] <- sample_names_clean

###############################################################################
## 6) Sanitize sample names in metadata
###############################################################################
meta_clean <- meta_raw %>%
  mutate(Sample = sanitize_ids(Sample))

###############################################################################
## 7) Keep only common samples (feature table ∩ metadata)
###############################################################################
common_samples <- intersect(colnames(feature_clean)[-1], meta_clean$Sample)

if (length(common_samples) < 3) {
  stop(
    "Too few overlapping samples between feature table and metadata after sanitization.\n",
    "Check that sample IDs match between inputs."
  )
}

feature_clean <- feature_clean %>%
  select(Feature, all_of(common_samples))

meta_clean <- meta_clean %>%
  filter(Sample %in% common_samples)

# Ensure metadata rows align with the sample order used in PCA
meta_clean <- meta_clean %>%
  mutate(Sample = factor(Sample, levels = common_samples)) %>%
  arrange(Sample) %>%
  mutate(Sample = as.character(Sample))

###############################################################################
## 8) Convert to numeric matrix
###############################################################################
data_mat <- as.matrix(feature_clean[, -1, drop = FALSE])
storage.mode(data_mat) <- "numeric"
rownames(data_mat) <- feature_clean$Feature

###############################################################################
## 9) Replace zeros (avoid log / LoD failures)
###############################################################################
min_positive <- suppressWarnings(min(data_mat[data_mat > 0], na.rm = TRUE))
if (!is.finite(min_positive)) {
  stop("No positive values found in the data matrix; cannot replace zeros or log-transform.")
}
data_mat[data_mat == 0] <- min_positive / 2

###############################################################################
## 10) Remove non-finite features
###############################################################################
finite_rows <- apply(data_mat, 1, function(x) all(is.finite(x)))
data_mat <- data_mat[finite_rows, , drop = FALSE]

###############################################################################
## 11) Remove zero-variance features
###############################################################################
feature_sd <- apply(data_mat, 1, sd)
data_mat <- data_mat[feature_sd > 0, , drop = FALSE]

###############################################################################
## 12) Remove zero-variance samples
###############################################################################
sample_sd <- apply(data_mat, 2, sd)
data_mat <- data_mat[, sample_sd > 0, drop = FALSE]

# If any samples were removed, update metadata accordingly
kept_samples <- colnames(data_mat)
meta_clean <- meta_clean %>% filter(Sample %in% kept_samples)

if (length(kept_samples) < 3) stop("Too few samples remain after variance filtering.")

###############################################################################
## 13) Log2 transform
###############################################################################
data_mat <- log2(data_mat)

###############################################################################
## 14) Pareto scaling (PCA-stable)
###############################################################################
data_mat <- apply(data_mat, 2, pareto_scale)

# Remove any features that became NA due to scaling edge cases
data_mat <- data_mat[apply(data_mat, 1, function(x) all(is.finite(x))), , drop = FALSE]

###############################################################################
## 15) PCA
###############################################################################
pca_res <- prcomp(t(data_mat), center = FALSE, scale. = FALSE)

pca_scores <- as.data.frame(pca_res$x) %>%
  rownames_to_column("Sample") %>%
  left_join(meta_clean, by = "Sample")

variance <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100

###############################################################################
## 16) PCA plot (colored by first metadata factor after Sample)
###############################################################################
if (ncol(meta_clean) < 2) {
  stop("Metadata must contain at least one grouping column in addition to 'Sample'.")
}
group_var <- colnames(meta_clean)[2]

pca_plot <- ggplot(
  pca_scores,
  aes(x = PC1, y = PC2, color = .data[[group_var]])
) +
  geom_point(size = point_size, alpha = 0.85) +
  stat_ellipse(level = ellipse_level, linewidth = 0.8) +
  labs(
    title = "PCA of Metabolomics Data",
    x = paste0("PC1 (", round(variance[1], 1), "%)"),
    y = paste0("PC2 (", round(variance[2], 1), "%)"),
    color = group_var
  ) +
  theme_classic(base_size = base_fontsize)

print(pca_plot)

# Save plot
ggsave(filename = output_pca_plot, plot = pca_plot, width = 8.5, height = 6, dpi = 300)

###############################################################################
## 17) Export MetaboAnalyst-ready tables
###############################################################################
final_feature_table <- data.frame(
  Feature = rownames(data_mat),
  data_mat,
  check.names = FALSE
)

write.csv(final_feature_table, output_feature_table, row.names = FALSE)
write.csv(meta_clean, output_metadata, row.names = FALSE)

###############################################################################
## Script sign-off
###############################################################################
# — Lee Chua

###############################################################################
## Citations
###############################################################################
# Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G,
#   Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K,
#   Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K,
#   Yutani H (2019). “Welcome to the tidyverse.” Journal of Open Source Software,
#   4(43), 1686. https://doi.org/10.21105/joss.01686.
#
# Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag
#   New York. https://ggplot2.tidyverse.org.
#
# R Core Team (2025). R: A Language and Environment for Statistical Computing.
#   R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/.
###############################################################################
