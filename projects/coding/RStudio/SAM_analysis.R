# ============================================================
# SAM METABOLOMICS ANALYSIS
# Author: You
# Purpose:
#   - Compare SAM levels across growth phase
#   - Normalize strains to matched WT
#   - Compute log2 fold-change (Log → Stationary)
#   - REMOVE all non-existent factor levels in plots
# ============================================================


# ------------------------------------------------------------
# 1. Load required libraries
# ------------------------------------------------------------
library(tidyverse)
library(ggpubr)


# ------------------------------------------------------------
# 2. Import data
# ------------------------------------------------------------
df <- read.csv("SAM_levels.csv")

# Quick structure check
glimpse(df)


# ------------------------------------------------------------
# 3. Define ordered factors
#    (Wild Type first, others alphabetical)
# ------------------------------------------------------------

group_levels <- df %>%
  distinct(Group) %>%
  arrange(Group) %>%
  pull(Group)

group_levels <- c("Wild Type", setdiff(group_levels, "Wild Type"))

df <- df %>%
  mutate(
    Group = factor(Group, levels = group_levels),
    Phase = factor(Phase, levels = c("Log", "Stationary"))
  )


# ------------------------------------------------------------
# 4. Summary statistics (mean + SE)
# ------------------------------------------------------------
df_summary <- df %>%
  group_by(Group, Media, Phase) %>%
  summarise(
    mean_SAM = mean(SAM_level, na.rm = TRUE),
    se_SAM   = sd(SAM_level, na.rm = TRUE) / sqrt(n()),
    n        = n(),
    .groups  = "drop"
  )


# ============================================================
# GRAPH 1 — GROUP × PHASE COMPARISON
# (Unused groups removed per Media automatically)
# ============================================================

p_phase <- ggplot(
  df_summary,
  aes(x = Group, y = mean_SAM, fill = Phase)
) +
  geom_col(
    position = position_dodge(width = 0.8)
  ) +
  geom_errorbar(
    aes(ymin = mean_SAM - se_SAM,
        ymax = mean_SAM + se_SAM),
    width = 0.2,
    position = position_dodge(width = 0.8)
  ) +
  facet_wrap(
    ~ Media,
    scales = "free_x"   # <<< KEY: drops unused groups
  ) +
  scale_x_discrete(drop = TRUE) +  # <<< removes empty space
  labs(
    title = "SAM Levels Across Growth Phases",
    x = "Group",
    y = "Mean SAM Level ± SE"
  ) +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_phase


# ============================================================
# GRAPH 2 — STRAINS VS MATCHED WILD TYPE
# ============================================================

# Extract WT reference per Media + Phase
wt_reference <- df_summary %>%
  filter(Group == "Wild Type") %>%
  select(Media, Phase, wt_mean = mean_SAM)

# Normalize strains to matched WT
df_vs_wt <- df_summary %>%
  left_join(wt_reference, by = c("Media", "Phase")) %>%
  mutate(
    ratio_to_WT = mean_SAM / wt_mean
  )

p_vs_wt <- ggplot(
  df_vs_wt,
  aes(x = Group, y = ratio_to_WT, fill = Phase)
) +
  geom_col(position = position_dodge(width = 0.8)) +
  facet_wrap(
    ~ Media,
    scales = "free_x"   # <<< removes empty group slots
  ) +
  scale_x_discrete(drop = TRUE) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    title = "SAM Levels Relative to Matched Wild Type",
    x = "Group",
    y = "Fold Change vs WT"
  ) +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_vs_wt


# ============================================================
# GRAPH 3 — LOG2 FOLD CHANGE (LOG → STATIONARY) VS WT
# ============================================================

log2fc <- df_summary %>%
  select(Group, Media, Phase, mean_SAM) %>%
  pivot_wider(names_from = Phase, values_from = mean_SAM) %>%
  filter(!is.na(Log), !is.na(Stationary)) %>%  # <<< remove non-existent pairs
  mutate(
    log2FC = log2(Stationary / Log)
  )

# WT baseline
wt_log2fc <- log2fc %>%
  filter(Group == "Wild Type") %>%
  select(Media, wt_log2FC = log2FC)

log2fc_norm <- log2fc %>%
  left_join(wt_log2fc, by = "Media") %>%
  mutate(
    delta_log2FC = log2FC - wt_log2FC
  )

p_log2fc <- ggplot(
  log2fc_norm,
  aes(x = Group, y = delta_log2FC)
) +
  geom_col(fill = "steelblue") +
  facet_wrap(
    ~ Media,
    scales = "free_x"
  ) +
  scale_x_discrete(drop = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Log2 Fold Change (Stationary vs Log) Relative to WT",
    x = "Group",
    y = expression(Delta~log[2]~Fold~Change)
  ) +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_log2fc


# ============================================================
# GRAPH 4 — RAW SAM DISTRIBUTIONS (QC)
# ============================================================

ggplot(
  df,
  aes(x = Group, y = SAM_level)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1) +
  facet_grid(
    Media ~ Phase,
    scales = "free_x"   # <<< critical
  ) +
  scale_x_discrete(drop = TRUE) +
  labs(
    title = "Raw SAM Distributions (QC View)",
    x = "Group",
    y = "SAM Level"
  ) +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

