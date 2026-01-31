###############################################################################
## Targeted IAA Analysis (WT-normalized Log2FC + Paired WT Phase Test + Summaries)
## - Cleans Excel-based IAA intensity table
## - Standardizes Condition / Group / Phase labels
## - Splits into two conditions (YMD+M and YMD+AAA+M)
## - Computes WT phase-specific references and per-sample log2FC vs WT
## - Generates condition-specific log2FC boxplots with statistical comparisons
## - Performs paired WT Stationary vs Log tests within each condition (Wilcoxon)
## - Produces mean ± SE bar plots (per condition, combined, and WT vs ego3Δ focus)
##
## Author: Lee Chua
## Last updated: 2026-01-31
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(ggsci)
  library(ggpubr)
  library(patchwork)
})

###############################################################################
## 0) Configuration
###############################################################################
# Input
input_excel <- "MassSpec.xlsx"
input_sheet <- "Fulldata"

# Output
out_dir_fig <- "figures"
dir.create(out_dir_fig, showWarnings = FALSE, recursive = TRUE)

# Standardize conditions
condition_map <- c("YMD+WYFM" = "YMD+AAA+M")
conditions_keep <- c("YMD+M", "YMD+AAA+M")

# Parsing controls
blank_label <- "Blank"
value_col   <- "IAA"

# Group normalization rules
# - Convert "...ko" -> "Δ..."
# - Normalize WT naming to "WT"
normalize_group <- function(x) {
  x %>%
    str_replace_all("ko", "Δ") %>%
    str_replace_all("\\bWt\\b", "WT")
}

# Phase normalization rules (accepts L/S and full strings)
normalize_phase <- function(x) {
  x <- toupper(trimws(x))
  case_when(
    x %in% c("L", "LOG") ~ "Log",
    x %in% c("S", "STATIONARY") ~ "Stationary",
    TRUE ~ x
  )
}

# Condition-specific suffix stripping
strip_suffix_by_condition <- function(group, condition) {
  case_when(
    condition == "YMD+M"      ~ str_replace_all(group, "a", ""),
    condition == "YMD+AAA+M"  ~ str_replace_all(group, "b", ""),
    TRUE ~ group
  )
}

# Optional strain relabeling rules (preserve your intent, but keep explicit)
relabel_special_strains <- function(group) {
  group %>%
    str_replace_all("aro1Δ", "aro1Δ_EGO3S288c") %>%
    str_replace_all("ptr3Δ", "ptr3Δ_EGO3S288c")
}

# Predefined plotting orders (edit as needed)
strain_order_M <- c("WT", "ego3Δ")
strain_order_AAA_M <- c("WT", "EGO3S288c", "aro1Δ_EGO3S288c", "ptr3Δ_EGO3S288c")

comparisons_M <- list(c("WT", "ego3Δ"))
comparisons_AAA <- list(
  c("WT", "EGO3S288c"),
  c("WT", "aro1Δ_EGO3S288c"),
  c("WT", "ptr3Δ_EGO3S288c")
)

###############################################################################
## 1) Load and clean data
###############################################################################
stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
}
stop_if_missing(input_excel)

raw <- read_excel(input_excel, sheet = input_sheet, col_names = FALSE)
if (ncol(raw) < 3) stop("Expected >= 3 columns: Condition, Sample, IAA.")

dat <- raw[, 1:3]
colnames(dat) <- c("Condition", "Sample", value_col)

fulldata <- dat %>%
  drop_na(Sample, .data[[value_col]]) %>%
  fill(Condition) %>%
  mutate(
    # robust numeric parsing
    !!value_col := parse_number(as.character(.data[[value_col]])),
    Group_raw   = str_extract(Sample, "^[^_]+"),
    Phase_raw   = str_extract(Sample, "_[LSls]|_[A-Za-z]+") %>%
      str_remove("^_") %>%
      str_remove("^_") %>%
      str_remove("_") %>%
      toupper() %>%
      trimws()
  ) %>%
  filter(!is.na(Group_raw), Group_raw != blank_label) %>%
  mutate(
    Condition = recode(Condition, !!!condition_map, .default = Condition),
    Group     = normalize_group(Group_raw),
    Phase     = normalize_phase(Phase_raw)
  ) %>%
  filter(Condition %in% conditions_keep, !is.na(.data[[value_col]]))

# Apply condition-specific group cleanup and preserve your special renaming
fulldata <- fulldata %>%
  mutate(
    Group = strip_suffix_by_condition(Group, Condition),
    Group = relabel_special_strains(Group),
    Phase = factor(Phase, levels = c("Log", "Stationary")),
    Condition = factor(Condition, levels = conditions_keep)
  )

###############################################################################
## 2) Split into conditions
###############################################################################
Fulldata_M <- fulldata %>% filter(Condition == "YMD+M")
Fulldata_AAA_M <- fulldata %>% filter(Condition == "YMD+AAA+M")

###############################################################################
## 3) Compute WT references by Phase and add log2FC vs WT (per condition)
###############################################################################
compute_wt_reference <- function(df, value_col = "IAA") {
  df %>%
    filter(Group == "WT") %>%
    group_by(Phase) %>%
    summarise(Reference = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop") %>%
    filter(!is.na(Reference), Reference > 0)
}

add_log2fc_vs_reference <- function(df, ref_df, value_col = "IAA") {
  df %>%
    left_join(ref_df, by = "Phase") %>%
    filter(!is.na(Reference), Reference > 0, .data[[value_col]] > 0) %>%
    mutate(Log2FC = log2(.data[[value_col]] / Reference)) %>%
    filter(is.finite(Log2FC))
}

ref_WT_M <- compute_wt_reference(Fulldata_M, value_col)
ref_WT_AAA_M <- compute_wt_reference(Fulldata_AAA_M, value_col)

Fulldata_M_fc <- add_log2fc_vs_reference(Fulldata_M, ref_WT_M, value_col)
Fulldata_AAA_M_fc <- add_log2fc_vs_reference(Fulldata_AAA_M, ref_WT_AAA_M, value_col)

###############################################################################
## 4) Condition-specific log2FC plots with t-tests (ggpubr)
###############################################################################
p_M <- ggplot(
  Fulldata_M_fc,
  aes(x = factor(Group, levels = strain_order_M), y = Log2FC, fill = Group)
) +
  geom_boxplot(outlier.shape = 16, outlier.size = 2) +
  facet_wrap(~Phase, scales = "free_y") +
  labs(title = "YMD+M", x = "Strain", y = "log2(FC vs WT)") +
  scale_fill_startrek() +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  stat_compare_means(
    comparisons = comparisons_M,
    label = "p.signif",
    method = "t.test",
    hide.ns = TRUE
  )

p_AAA_M <- ggplot(
  Fulldata_AAA_M_fc,
  aes(x = factor(Group, levels = strain_order_AAA_M), y = Log2FC, fill = Group)
) +
  geom_boxplot(outlier.shape = 16, outlier.size = 2) +
  facet_wrap(~Phase, scales = "free_y") +
  labs(title = "YMD+AAA+M", x = "Strain", y = "log2(FC vs WT)") +
  scale_fill_d3() +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  stat_compare_means(
    comparisons = comparisons_AAA,
    label = "p.signif",
    method = "t.test",
    hide.ns = TRUE
  )

###############################################################################
## 5) Paired WT Stationary vs Log test within each condition (Wilcoxon)
###############################################################################
# We rely on a numeric replicate identifier embedded in Sample.
# If your Sample IDs encode replicates differently, adjust SampleNum regex.
extract_sample_num <- function(x) str_extract(x, "[0-9]+")

paired_wt_phase_fc <- function(df_wt, value_col = "IAA") {
  df_wt %>%
    mutate(
      SampleNum = extract_sample_num(Sample),
      Phase = factor(Phase, levels = c("Log", "Stationary"))
    ) %>%
    filter(!is.na(SampleNum)) %>%
    select(SampleNum, Phase, !!sym(value_col)) %>%
    pivot_wider(names_from = Phase, values_from = !!sym(value_col)) %>%
    drop_na(Log, Stationary) %>%
    mutate(Log2FC = log2(Stationary / Log))
}

WT_M <- Fulldata_M %>% filter(Group == "WT")
WT_AAA <- Fulldata_AAA_M %>% filter(Group == "WT")

df_M_wide <- paired_wt_phase_fc(WT_M, value_col)
df_AAA_wide <- paired_wt_phase_fc(WT_AAA, value_col)

test_M <- if (nrow(df_M_wide) >= 2) wilcox.test(df_M_wide$Stationary, df_M_wide$Log, paired = TRUE) else NULL
test_AAA <- if (nrow(df_AAA_wide) >= 2) wilcox.test(df_AAA_wide$Stationary, df_AAA_wide$Log, paired = TRUE) else NULL

df_combined <- bind_rows(
  df_AAA_wide %>% mutate(Condition = "YMD+AAA+M"),
  df_M_wide %>% mutate(Condition = "YMD+M")
) %>%
  mutate(Condition = factor(Condition, levels = conditions_keep))

p_paired <- ggplot(df_combined, aes(x = Condition, y = Log2FC, fill = Condition)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 2, width = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  stat_compare_means(method = "wilcox.test", paired = TRUE, label = "p.format") +
  scale_fill_startrek() +
  labs(
    title = "WT: Log2 Fold Change in IAA (Stationary vs Log Phase)",
    y = "Log2 Fold Change (Stationary / Log)",
    x = "Condition"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

###############################################################################
## 6) Mean ± SE bar plots (per condition, combined)
###############################################################################
summary_se <- function(data, measurevar, groupvars) {
  data %>%
    group_by(across(all_of(groupvars))) %>%
    summarise(
      mean = mean(.data[[measurevar]], na.rm = TRUE),
      se   = sd(.data[[measurevar]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[measurevar]]))),
      n    = sum(!is.na(.data[[measurevar]])),
      .groups = "drop"
    )
}

IAA_summary_M <- summary_se(Fulldata_M, value_col, c("Group", "Phase")) %>%
  mutate(Group = factor(Group, levels = strain_order_M))

IAA_summary_AAA_M <- summary_se(Fulldata_AAA_M, value_col, c("Group", "Phase")) %>%
  mutate(Group = factor(Group, levels = strain_order_AAA_M))

bar_M <- ggplot(IAA_summary_M, aes(x = Group, y = mean, fill = Group)) +
  geom_col(position = position_dodge(0.9), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.2, position = position_dodge(0.9)
  ) +
  facet_wrap(~Phase, scales = "free_y") +
  scale_fill_startrek() +
  labs(
    title = "Average IAA Intensity (YMD+M)",
    y = "Mean IAA Intensity ± SE",
    x = "Strain"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

bar_AAA_M <- ggplot(IAA_summary_AAA_M, aes(x = Group, y = mean, fill = Group)) +
  geom_col(position = position_dodge(0.9), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.2, position = position_dodge(0.9)
  ) +
  facet_wrap(~Phase, scales = "free_y") +
  scale_fill_d3() +
  labs(
    title = "Average IAA Intensity (YMD+AAA+M)",
    y = "Mean IAA Intensity ± SE",
    x = "Strain"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

IAA_summary_combined <- bind_rows(
  summary_se(Fulldata_M, value_col, c("Group", "Phase")) %>% mutate(Condition = "YMD+M"),
  summary_se(Fulldata_AAA_M, value_col, c("Group", "Phase")) %>% mutate(Condition = "YMD+AAA+M")
) %>%
  mutate(
    Group = factor(Group, levels = unique(c(strain_order_M, strain_order_AAA_M))),
    Condition = factor(Condition, levels = conditions_keep)
  )

bar_combined <- ggplot(
  IAA_summary_combined,
  aes(x = Group, y = mean, fill = Condition)
) +
  geom_col(position = position_dodge(0.9), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    position = position_dodge(0.9),
    width = 0.2
  ) +
  facet_wrap(~Phase, scales = "free_y") +
  scale_fill_startrek() +
  labs(
    title = "Average IAA Intensity ± SE Across Conditions",
    x = "Strain",
    y = "Mean IAA Intensity ± SE",
    fill = "Condition"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

###############################################################################
## 7) Focus plot: WT vs ego3Δ mean ± SE by Phase with significance
###############################################################################
IAA_two <- bind_rows(
  Fulldata_M %>% mutate(Condition = "YMD+M"),
  Fulldata_AAA_M %>% mutate(Condition = "YMD+AAA+M")
) %>%
  filter(Group %in% c("WT", "ego3Δ"))

IAA_summary_focus <- IAA_two %>%
  group_by(Phase, Group) %>%
  summarise(
    mean = mean(.data[[value_col]], na.rm = TRUE),
    se   = sd(.data[[value_col]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[value_col]]))),
    .groups = "drop"
  ) %>%
  mutate(
    Phase = factor(Phase, levels = c("Log", "Stationary")),
    Group = factor(Group, levels = c("WT", "ego3Δ"))
  )

p_values <- IAA_two %>%
  group_by(Phase) %>%
  summarise(
    p = if (n_distinct(Group) == 2) t.test(.data[[value_col]] ~ Group)$p.value else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(sig = case_when(
    is.na(p)     ~ "ns",
    p <= 0.001   ~ "***",
    p <= 0.01    ~ "**",
    p <= 0.05    ~ "*",
    TRUE         ~ "ns"
  ))

# label height per phase (phase-specific maxima + padding)
y_lims <- IAA_summary_focus %>%
  group_by(Phase) %>%
  summarise(y = max(mean + se, na.rm = TRUE) * 1.07, .groups = "drop")

sig_labels <- p_values %>%
  left_join(y_lims, by = "Phase")

bar_WT_vs_ego3 <- ggplot(IAA_summary_focus, aes(x = Phase, y = mean, fill = Group)) +
  geom_col(position = position_dodge(0.8), width = 0.6) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.2, position = position_dodge(0.8)
  ) +
  geom_text(
    data = sig_labels,
    aes(x = Phase, y = y, label = sig),
    inherit.aes = FALSE,
    size = 6,
    vjust = 0
  ) +
  scale_fill_startrek(name = "Strain") +
  labs(
    title = "WT vs ego3Δ: IAA Intensity by Phase",
    x = "Phase",
    y = "IAA Intensity"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

###############################################################################
## 8) Print plots and tests
###############################################################################
print(p_M)
print(p_AAA_M)
print(bar_M)
print(bar_AAA_M)
print(p_paired)
print(bar_combined)
print(bar_WT_vs_ego3)

cat("Wilcoxon paired test for WT (YMD+AAA+M):\n")
print(test_AAA)

cat("\nWilcoxon paired test for WT (YMD+M):\n")
print(test_M)

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
# Wickham H, Bryan J (2025). readxl: Read Excel Files. R package.
#   https://readxl.tidyverse.org.
#
# Xiao N (2025). ggsci: Scientific Journal and Sci-Fi Themed Color Palettes for
#   'ggplot2'. R package. https://cran.r-project.org/package=ggsci
#
# Kassambara A (2025). ggpubr: 'ggplot2' Based Publication Ready Plots.
#   R package. https://cran.r-project.org/package=ggpubr
#
# Pedersen TL (2025). patchwork: The Composer of Plots. R package.
#   https://patchwork.data-imaginist.com/
###############################################################################
