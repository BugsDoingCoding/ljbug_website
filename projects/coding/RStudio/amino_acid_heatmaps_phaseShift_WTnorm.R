###############################################################################
## Amino Acid Heatmaps: Phase Shift + WT-Normalized Comparisons (R)
## - Computes mean amino acid abundance per Group within each medium
## - Derives within-strain phase shift: log2(Stationary / Log)
## - Derives WT-normalized shifts per phase: log2(Strain / WT)
## - Produces publication-ready heatmaps and combined multi-panel figures
##
## Author: Lee Chua
## Last updated: 2026-01-31
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
})

###############################################################################
## 0) Configuration
###############################################################################
# Input
input_csv <- "amino_acid_abundance_table.csv"

# Expected column conventions
# - Media: growth medium label
# - Group: composite key encoding strain and phase (e.g., "StrainX_log")
# - Numeric columns: amino acid measurements
media_col <- "Media"
group_col <- "Group"

# Define media names as they appear in your dataset
media_A <- "YMD"
media_B <- "YMD+WYFM"

# Define WT matching pattern (case-insensitive)
wt_pattern <- "WildType"

# Output
out_dir_fig <- "figures"
dir.create(out_dir_fig, showWarnings = FALSE, recursive = TRUE)

# Heatmap label text
label_phase_shift <- "log2(Stationary / Log)"
label_vs_wt        <- "log2(Strain / WT)"

###############################################################################
## 1) Helpers
###############################################################################
stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
}

mean_by_group <- function(df) {
  # Mean across numeric columns per composite Group, then split Group -> Strain, Phase
  df %>%
    group_by(.data[[group_col]]) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    separate(
      col  = all_of(group_col),
      into = c("Strain", "Phase"),
      sep  = "_",
      remove = FALSE,
      extra = "merge",
      fill  = "right"
    )
}

compute_log2fc_stationary_over_log <- function(mean_df) {
  # Wide -> long -> wide by phase to compute log2(stationary/log)
  mean_df %>%
    pivot_longer(
      cols      = where(is.numeric),
      names_to  = "AminoAcid",
      values_to = "MeanAbundance"
    ) %>%
    select(Strain, Phase, AminoAcid, MeanAbundance) %>%
    pivot_wider(names_from = Phase, values_from = MeanAbundance) %>%
    mutate(log2FC = log2(stationary / log))
}

compute_vs_wt <- function(mean_df, wt_pattern = "WildType") {
  # Compute log2(strain / WT) within each Phase x AminoAcid
  mean_long <- mean_df %>%
    pivot_longer(
      cols      = where(is.numeric),
      names_to  = "AminoAcid",
      values_to = "MeanAbundance"
    )
  
  wt_values <- mean_long %>%
    filter(grepl(wt_pattern, Strain, ignore.case = TRUE)) %>%
    select(Phase, AminoAcid, WT_value = MeanAbundance)
  
  mean_long %>%
    left_join(wt_values, by = c("Phase", "AminoAcid")) %>%
    mutate(log2FC_vs_WT = log2(MeanAbundance / WT_value))
}

make_strain_levels <- function(df, wt_pattern = "WildType") {
  strains <- df %>% distinct(Strain) %>% pull(Strain)
  c(
    sort(strains[grepl(wt_pattern, strains, ignore.case = TRUE)]),
    sort(strains[!grepl(wt_pattern, strains, ignore.case = TRUE)])
  )
}

make_aa_levels <- function(df) {
  df %>%
    distinct(AminoAcid) %>%
    pull(AminoAcid) %>%
    sort(decreasing = TRUE)
}

apply_factors <- function(df, strain_levels, aa_levels) {
  df %>%
    mutate(
      Strain    = factor(Strain, levels = strain_levels),
      AminoAcid = factor(AminoAcid, levels = aa_levels)
    )
}

plot_heatmap <- function(df, fill_col, title_text, fill_label) {
  ggplot(df, aes(x = Strain, y = AminoAcid, fill = .data[[fill_col]])) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient2(
      low = "blue",
      mid = "grey",
      high = "red",
      midpoint = 0,
      na.value = "grey90",
      name = fill_label
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      panel.grid  = element_blank()
    ) +
    labs(title = title_text, x = "Strain", y = "Amino Acid")
}

save_plot <- function(p, filename, width = 10, height = 6, dpi = 300) {
  ggsave(
    filename = file.path(out_dir_fig, filename),
    plot     = p,
    width    = width,
    height   = height,
    dpi      = dpi
  )
}

###############################################################################
## 2) Load data
###############################################################################
stop_if_missing(input_csv)

df <- read.csv(input_csv, stringsAsFactors = FALSE)

required_cols <- c(media_col, group_col)
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

###############################################################################
## 3) Split by medium
###############################################################################
df_A <- df %>% filter(.data[[media_col]] == media_A)
df_B <- df %>% filter(.data[[media_col]] == media_B)

if (nrow(df_A) == 0) stop("No rows found for medium: ", media_A)
if (nrow(df_B) == 0) stop("No rows found for medium: ", media_B)

###############################################################################
## 4) Mean abundance per Group (Strain_Phase) within each medium
###############################################################################
mean_A <- mean_by_group(df_A)
mean_B <- mean_by_group(df_B)

###############################################################################
## 5) Phase shift: log2(Stationary / Log)
###############################################################################
log2fc_A <- compute_log2fc_stationary_over_log(mean_A)
log2fc_B <- compute_log2fc_stationary_over_log(mean_B)

heatmap_phase_A <- log2fc_A %>% select(Strain, AminoAcid, log2FC)
heatmap_phase_B <- log2fc_B %>% select(Strain, AminoAcid, log2FC)

strain_levels_A <- make_strain_levels(heatmap_phase_A, wt_pattern)
strain_levels_B <- make_strain_levels(heatmap_phase_B, wt_pattern)

aa_levels_A <- make_aa_levels(heatmap_phase_A)
aa_levels_B <- make_aa_levels(heatmap_phase_B)

heatmap_phase_A <- apply_factors(heatmap_phase_A, strain_levels_A, aa_levels_A)
heatmap_phase_B <- apply_factors(heatmap_phase_B, strain_levels_B, aa_levels_B)

p_phase_A <- plot_heatmap(heatmap_phase_A, "log2FC", title_text = media_A, fill_label = label_phase_shift)
p_phase_B <- plot_heatmap(heatmap_phase_B, "log2FC", title_text = "WYFM", fill_label = label_phase_shift)

combined_phase_shift <- p_phase_A + p_phase_B +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Amino Acid Abundance: Log vs Stationary (Within-Strain)",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

print(combined_phase_shift)
save_plot(combined_phase_shift, "AA_phase_shift_heatmaps.png", width = 12, height = 6)

###############################################################################
## 6) WT-normalized: log2(Strain / WT) within each phase (Log, Stationary)
###############################################################################
vsWT_A <- compute_vs_wt(mean_A, wt_pattern)
vsWT_B <- compute_vs_wt(mean_B, wt_pattern)

# Split by phase
heat_A_log <- vsWT_A %>% filter(Phase == "log") %>% select(Strain, AminoAcid, log2FC_vs_WT)
heat_A_sta <- vsWT_A %>% filter(Phase == "stationary") %>% select(Strain, AminoAcid, log2FC_vs_WT)

heat_B_log <- vsWT_B %>% filter(Phase == "log") %>% select(Strain, AminoAcid, log2FC_vs_WT)
heat_B_sta <- vsWT_B %>% filter(Phase == "stationary") %>% select(Strain, AminoAcid, log2FC_vs_WT)

# Apply same factor ordering as phase-shift plots for consistency
heat_A_log <- apply_factors(heat_A_log, strain_levels_A, aa_levels_A)
heat_A_sta <- apply_factors(heat_A_sta, strain_levels_A, aa_levels_A)
heat_B_log <- apply_factors(heat_B_log, strain_levels_B, aa_levels_B)
heat_B_sta <- apply_factors(heat_B_sta, strain_levels_B, aa_levels_B)

# Remove WT after normalization (used only as baseline)
heat_A_log <- heat_A_log %>% filter(!grepl(wt_pattern, Strain, ignore.case = TRUE))
heat_A_sta <- heat_A_sta %>% filter(!grepl(wt_pattern, Strain, ignore.case = TRUE))
heat_B_log <- heat_B_log %>% filter(!grepl(wt_pattern, Strain, ignore.case = TRUE))
heat_B_sta <- heat_B_sta %>% filter(!grepl(wt_pattern, Strain, ignore.case = TRUE))

p_A_log <- plot_heatmap(heat_A_log, "log2FC_vs_WT", title_text = paste0(media_A, " – Log"), fill_label = label_vs_wt)
p_A_sta <- plot_heatmap(heat_A_sta, "log2FC_vs_WT", title_text = paste0(media_A, " – Stationary"), fill_label = label_vs_wt)
p_B_log <- plot_heatmap(heat_B_log, "log2FC_vs_WT", title_text = "WYFM – Log", fill_label = label_vs_wt)
p_B_sta <- plot_heatmap(heat_B_sta, "log2FC_vs_WT", title_text = "WYFM – Stationary", fill_label = label_vs_wt)

combined_vs_WT <- (p_A_log + p_A_sta) / (p_B_log + p_B_sta) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Amino Acid Abundance: WT-Normalized (Strain vs WT)",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

print(combined_vs_WT)
save_plot(combined_vs_WT, "AA_WT_normalized_heatmaps.png", width = 12, height = 10)

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
# Pedersen TL (2025). patchwork: The Composer of Plots. R package.
#   https://patchwork.data-imaginist.com/.
###############################################################################
