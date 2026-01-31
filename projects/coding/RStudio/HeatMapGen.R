# 1. Libraries and WD --------------------------------------------------------
# Load tidyverse for data manipulation and ggplot2 for plotting
library(tidyverse)
library(ggplot2)
library(patchwork)   # Added ONLY to combine plots for publication

# Set working directory to project folder
setwd("C:/your/file/directory")


# 2. Create Data Frames ------------------------------------------------------
# Read in metabolomics amino acid dataset

df <- read.csv("Metabolomics_AA_conc.csv")

# Split dataset by growth medium
df_WYFM <- df %>%
  filter(Media == "YMD+WYFM")

df_YMD <- df %>%
  filter(Media == "YMD")


# 3. Data Processing ---------------------------------------------------------
# Compute mean amino acid peak areas per Group (strain_phase)
# Then split Group into Strain and Phase metadata columns

mean_YMD <- df_YMD %>%
  group_by(Group) %>%
  summarise(
    across(
      where(is.numeric),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  separate(
    Group,
    into = c("Strain", "Phase"),
    sep = "_",
    remove = FALSE
  )

mean_WYFM <- df_WYFM %>%
  group_by(Group) %>%
  summarise(
    across(
      where(is.numeric),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  separate(
    Group,
    into = c("Strain", "Phase"),
    sep = "_",
    remove = FALSE
  )


# 4. Calculate Log Fold Change -----------------------------------------------
# Helper function calculates log2(stationary / log) per strain
# Converts wide AA data to long format, then back to wide by phase

compute_log2fc <- function(mean_df) {
  mean_df %>%
    pivot_longer(
      cols = where(is.numeric),
      names_to = "AminoAcid",
      values_to = "MeanPeakArea"
    ) %>%
    select(Strain, Phase, AminoAcid, MeanPeakArea) %>%
    pivot_wider(
      names_from = Phase,
      values_from = MeanPeakArea
    ) %>%
    mutate(
      log2FC = log2(stationary / log)
    )
}

# Apply phase fold-change calculation per medium
log2fc_YMD  <- compute_log2fc(mean_YMD)
log2fc_WYFM <- compute_log2fc(mean_WYFM)


# 5. Final Processing --------------------------------------------------------
# Reduce to only columns needed for heatmap plotting

heatmap_YMD <- log2fc_YMD %>%
  select(Strain, AminoAcid, log2FC)

heatmap_WYFM <- log2fc_WYFM %>%
  select(Strain, AminoAcid, log2FC)

# Define strain ordering: Wild Type first, then mutants alphabetically
strain_levels_YMD <- heatmap_YMD %>%
  distinct(Strain) %>%
  pull(Strain) %>%
  {c(
    sort(.[grepl("WildType", .)]),
    sort(.[!grepl("WildType", .)])
  )}

strain_levels_WYFM <- heatmap_WYFM %>%
  distinct(Strain) %>%
  pull(Strain) %>%
  {c(
    sort(.[grepl("WildType", .)]),
    sort(.[!grepl("WildType", .)])
  )}

# Alphabetical amino acid ordering (top → bottom in heatmap)
aa_levels_YMD <- heatmap_YMD %>%
  distinct(AminoAcid) %>%
  pull(AminoAcid) %>%
  sort(decreasing = TRUE)

aa_levels_WYFM <- heatmap_WYFM %>%
  distinct(AminoAcid) %>%
  pull(AminoAcid) %>%
  sort(decreasing = TRUE)

# Apply factor levels for consistent plotting
heatmap_YMD <- heatmap_YMD %>%
  mutate(AminoAcid = factor(AminoAcid, levels = aa_levels_YMD)) %>%
  mutate(Strain = factor(Strain, levels = strain_levels_YMD))

heatmap_WYFM <- heatmap_WYFM %>%
  mutate(AminoAcid = factor(AminoAcid, levels = aa_levels_WYFM)) %>%
  mutate(Strain = factor(Strain, levels = strain_levels_WYFM))


# 6. Generate Log v Stat Heatmaps -------------------------------------------
# These plots visualize phase-dependent changes within each strain

plot1 <- ggplot(heatmap_YMD, aes(x = Strain, y = AminoAcid, fill = log2FC)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "blue",
    mid = "grey",
    high = "red",
    midpoint = 0,
    na.value = "grey90",
    name = "log2(Stationary / Log)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    title = "YMD",
    x = "Strain",
    y = "Amino Acid"
  )

plot2 <- ggplot(heatmap_WYFM, aes(x = Strain, y = AminoAcid, fill = log2FC)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "blue",
    mid = "grey",
    high = "red",
    midpoint = 0,
    na.value = "grey90",
    name = "log2(Stationary / Log)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    title = "WYFM",
    x = "Strain",
    y = "Amino Acid"
  )


# 7. Compare Strains to WT in Log & Stationary ------------------------------
# Compute log2(strain / WT) within each phase and medium

compute_vs_WT <- function(mean_df, wt_pattern = "WildType") {
  
  mean_long <- mean_df %>%
    pivot_longer(
      cols = where(is.numeric),
      names_to = "AminoAcid",
      values_to = "MeanPeakArea"
    )
  
  wt_values <- mean_long %>%
    filter(grepl(wt_pattern, Strain)) %>%
    select(Phase, AminoAcid, WT_value = MeanPeakArea)
  
  mean_long %>%
    left_join(
      wt_values,
      by = c("Phase", "AminoAcid")
    ) %>%
    mutate(
      log2FC_vs_WT = log2(MeanPeakArea / WT_value)
    )
}

# Compute WT-normalized values per medium
vsWT_YMD  <- compute_vs_WT(mean_YMD)
vsWT_WYFM <- compute_vs_WT(mean_WYFM)


# Prepare Heatmap Data ----------------------------------------------------
# Split WT-normalized data by phase and medium

heatmap_YMD_log <- vsWT_YMD %>%
  filter(Phase == "log") %>%
  select(Strain, AminoAcid, log2FC_vs_WT)

heatmap_YMD_stationary <- vsWT_YMD %>%
  filter(Phase == "stationary") %>%
  select(Strain, AminoAcid, log2FC_vs_WT)

heatmap_WYFM_log <- vsWT_WYFM %>%
  filter(Phase == "log") %>%
  select(Strain, AminoAcid, log2FC_vs_WT)

heatmap_WYFM_stationary <- vsWT_WYFM %>%
  filter(Phase == "stationary") %>%
  select(Strain, AminoAcid, log2FC_vs_WT)

# Reapply factor ordering for consistent axes
apply_factors <- function(df, strain_levels, aa_levels) {
  df %>%
    mutate(
      Strain = factor(Strain, levels = strain_levels),
      AminoAcid = factor(AminoAcid, levels = aa_levels)
    )
}

heatmap_YMD_log        <- apply_factors(heatmap_YMD_log, strain_levels_YMD, aa_levels_YMD)
heatmap_YMD_stationary <- apply_factors(heatmap_YMD_stationary, strain_levels_YMD, aa_levels_YMD)

heatmap_WYFM_log        <- apply_factors(heatmap_WYFM_log, strain_levels_WYFM, aa_levels_WYFM)
heatmap_WYFM_stationary <- apply_factors(heatmap_WYFM_stationary, strain_levels_WYFM, aa_levels_WYFM)

# Remove Wild Type AFTER normalization
# WT is used as the reference but excluded from visualization

heatmap_YMD_log <- heatmap_YMD_log %>%
  filter(!grepl("WildType", Strain))

heatmap_YMD_stationary <- heatmap_YMD_stationary %>%
  filter(!grepl("WildType", Strain))

heatmap_WYFM_log <- heatmap_WYFM_log %>%
  filter(!grepl("WildType", Strain))

heatmap_WYFM_stationary <- heatmap_WYFM_stationary %>%
  filter(!grepl("WildType", Strain))


# Plot Heatmaps -----------------------------------------------------------
# WT-normalized strain comparisons

plot3 <- ggplot(heatmap_YMD_log, aes(x = Strain, y = AminoAcid, fill = log2FC_vs_WT)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "blue",
    mid = "grey",
    high = "red",
    midpoint = 0,
    na.value = "grey90",
    name = "log2(Strain / WT)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    title = "YMD – Log",
    x = "Strain",
    y = "Amino Acid"
  )

plot4 <- ggplot(heatmap_YMD_stationary, aes(x = Strain, y = AminoAcid, fill = log2FC_vs_WT)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "blue",
    mid = "grey",
    high = "red",
    midpoint = 0,
    na.value = "grey90",
    name = "log2(Strain / WT)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    title = "YMD – Stationary",
    x = "Strain",
    y = "Amino Acid"
  )

plot5 <- ggplot(heatmap_WYFM_log, aes(x = Strain, y = AminoAcid, fill = log2FC_vs_WT)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "blue",
    mid = "grey",
    high = "red",
    midpoint = 0,
    na.value = "grey90",
    name = "log2(Strain / WT)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    title = "WYFM – Log",
    x = "Strain",
    y = "Amino Acid"
  )

plot6 <- ggplot(heatmap_WYFM_stationary, aes(x = Strain, y = AminoAcid, fill = log2FC_vs_WT)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "blue",
    mid = "grey",
    high = "red",
    midpoint = 0,
    na.value = "grey90",
    name = "log2(Strain / WT)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    title = "WYFM – Stationary",
    x = "Strain",
    y = "Amino Acid"
  )


# 8. Combine Plots for Publication -----------------------------------------
# Arrange all heatmaps into a single multi-panel figure
# Shared legend improves readability and consistency

combined_phase_shift <- plot1 + plot2 +
  plot_layout(guides = "collect") + 
  plot_annotation(
    title = "Amino Acid Abundance Log v Stationary",
    theme = theme(
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5
      )))
  

combined_vs_WT <- (plot3 + plot4) /
  (plot5 + plot6) +
  plot_layout(guides = "collect") + 
  plot_annotation(
    title = "Amino Acid Abundance WT v Mu",
    theme = theme(
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5
      )))


# Display combined figures
combined_phase_shift
combined_vs_WT

