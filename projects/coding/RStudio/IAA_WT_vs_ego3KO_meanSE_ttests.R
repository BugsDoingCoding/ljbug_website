###############################################################################
## Targeted IAA Analysis: WT vs ego3 Knockout(s) with Mean ± SE and Significance
## - Cleans Excel-based intensity table
## - Standardizes group/phase labels
## - Filters to WT + ego3-related genotypes detected in the dataset
## - Computes mean ± SE per Group × Condition × Phase
## - Performs pairwise two-sided t-tests: (KO vs WT) within each Condition × Phase
## - Produces a faceted bar plot with significance annotations
##
## Author: Lee Chua
## Last updated: 2026-01-31
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(ggsci)
  library(ggpubr)
})

###############################################################################
## 0) Configuration
###############################################################################
# Inputs
input_excel <- "targeted_metabolite_table.xlsx"  # e.g., MassSpec.xlsx
input_sheet <- "Fulldata"

# Column semantics
value_col   <- "IAA"      # targeted metabolite column name (as created after import)
blank_label <- "Blank"    # samples/groups to exclude

# Condition normalization (optional)
# Map any incoming condition names to standardized labels
condition_map <- c(
  "YMD+WYFM" = "YMD+AAA+M"
)

# Group normalization (optional)
# - Convert "...ko" -> "Δ..."
# - Normalize WT naming (e.g., "Wt" -> "WT")
normalize_group <- function(x) {
  x %>%
    str_replace_all("ko", "Δ") %>%
    str_replace_all("\\bWt\\b", "WT")
}

# Detect ego3-related groups using this pattern
ego3_regex <- "(?i)ego3"

# Plot output
out_dir_fig <- "figures"
dir.create(out_dir_fig, showWarnings = FALSE, recursive = TRUE)
out_plot    <- file.path(out_dir_fig, "IAA_WT_vs_ego3_mean_SE_ttests.png")

###############################################################################
## 1) Load and clean data
###############################################################################
stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
}

stop_if_missing(input_excel)

raw <- read_excel(input_excel, sheet = input_sheet, col_names = FALSE)

if (ncol(raw) < 3) {
  stop("Expected >= 3 columns in the Excel sheet: Condition, Sample, Value.")
}

# Standardize to exactly 3 columns (first 3), then name them
dat <- raw[, 1:3]
colnames(dat) <- c("Condition", "Sample", value_col)

clean <- dat %>%
  drop_na(Sample, .data[[value_col]]) %>%
  fill(Condition) %>%
  mutate(
    # robust numeric parsing (handles strings like "123,456" or "IAA=123")
    !!value_col := parse_number(as.character(.data[[value_col]])),
    Group_raw   = str_extract(Sample, "^[^_]+"),
    Phase_raw   = str_extract(Sample, "_[LSls]") %>%
      str_remove("_") %>%
      toupper() %>%
      trimws()
  ) %>%
  filter(!is.na(Group_raw), Group_raw != blank_label) %>%
  mutate(
    # Standardize condition labels
    Condition = recode(Condition, !!!condition_map, .default = Condition),
    
    # Standardize group and phase
    Group = normalize_group(Group_raw),
    Phase = case_when(
      Phase_raw %in% c("L", "LOG") ~ "Log",
      Phase_raw %in% c("S", "STATIONARY") ~ "Stationary",
      TRUE ~ Phase_raw
    ),
    Phase = factor(Phase, levels = c("Log", "Stationary"))
  ) %>%
  filter(!is.na(.data[[value_col]]))

###############################################################################
## 2) Optional: condition-specific cleanup (kept explicit, not hidden)
###############################################################################
# If your raw data encodes condition-specific suffix letters ("a"/"b") in Group,
# you can standardize them here. This matches the intent of your original code
# but keeps it transparent.
#
# Example (uncomment if applicable):
# clean <- clean %>%
#   mutate(
#     Group = case_when(
#       Condition == "YMD+M"       ~ str_replace_all(Group, "a", ""),
#       Condition == "YMD+AAA+M"   ~ str_replace_all(Group, "b", ""),
#       TRUE ~ Group
#     )
#   )

# If you have specific strain relabeling rules you want preserved, keep them here:
# clean <- clean %>%
#   mutate(
#     Group = str_replace_all(Group, "aro1Δ", "aro1Δ_EGO3S288c"),
#     Group = str_replace_all(Group, "ptr3Δ", "ptr3Δ_EGO3S288c")
#   )

# Ensure Condition is a factor with stable ordering (edit as needed)
clean <- clean %>%
  mutate(
    Condition = factor(Condition, levels = sort(unique(as.character(Condition))))
  )

###############################################################################
## 3) Identify WT and ego3 groups to plot
###############################################################################
groups_present <- unique(clean$Group)

ego3_groups <- unique(groups_present[grepl(ego3_regex, groups_present)])
groups_to_plot <- unique(c("WT", ego3_groups))
groups_to_plot <- groups_to_plot[groups_to_plot %in% groups_present]

if (length(groups_to_plot) < 2) {
  stop(
    "No ego3-related groups detected (pattern: ", ego3_regex, ").\n",
    "Check your group naming conventions or update ego3_regex."
  )
}

plot_data <- clean %>%
  filter(Group %in% groups_to_plot)

###############################################################################
## 4) Summary stats: mean ± SE per Group × Condition × Phase
###############################################################################
summary_df <- plot_data %>%
  group_by(Group, Condition, Phase) %>%
  summarise(
    mean = mean(.data[[value_col]], na.rm = TRUE),
    se   = sd(.data[[value_col]], na.rm = TRUE) /
      sqrt(sum(!is.na(.data[[value_col]]))),
    n    = sum(!is.na(.data[[value_col]])),
    .groups = "drop"
  )

###############################################################################
## 5) Pairwise t-tests: each KO vs WT within Condition × Phase
###############################################################################
wt_data <- plot_data %>% filter(Group == "WT")

stat_tests <- expand_grid(
  Group     = setdiff(groups_to_plot, "WT"),
  Condition = levels(plot_data$Condition),
  Phase     = levels(plot_data$Phase)
) %>%
  rowwise() %>%
  mutate(
    wt_vals = list(wt_data %>% filter(Condition == Condition, Phase == Phase) %>% pull(.data[[value_col]])),
    ko_vals = list(plot_data %>% filter(Group == Group, Condition == Condition, Phase == Phase) %>% pull(.data[[value_col]])),
    p.value = {
      if (length(wt_vals) >= 2 && length(ko_vals) >= 2) {
        tt <- try(t.test(ko_vals, wt_vals, alternative = "two.sided"), silent = TRUE)
        if (inherits(tt, "try-error")) NA_real_ else tt$p.value
      } else {
        NA_real_
      }
    },
    p.label = case_when(
      is.na(p.value)     ~ NA_character_,
      p.value < 0.001    ~ "***",
      p.value < 0.01     ~ "**",
      p.value < 0.05     ~ "*",
      TRUE               ~ "ns"
    )
  ) %>%
  ungroup() %>%
  select(Group, Condition, Phase, p.value, p.label)

# Determine y-positions for significance labels (placed above KO mean + SE)
# Use a single padding value derived from the overall dynamic range.
global_range <- range(summary_df$mean, finite = TRUE)
pad <- 0.05 * diff(global_range)

stat_tests <- stat_tests %>%
  left_join(summary_df, by = c("Group", "Condition", "Phase")) %>%
  mutate(ypos = mean + se + pad)

###############################################################################
## 6) Plot: Faceted bar plot (mean ± SE) + significance labels
###############################################################################
plot_df <- summary_df %>% filter(Group %in% groups_to_plot)

p <- ggplot(plot_df, aes(x = Phase, y = mean, fill = Condition)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.6,
    color = "black"
  ) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  facet_wrap(~ Group, nrow = 1) +
  scale_fill_startrek(name = "Condition") +
  labs(
    title = "Mean IAA Intensity ± SE — WT and ego3 Knockout(s)",
    x = "Growth Phase",
    y = "Mean IAA Intensity (± SE)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 0, vjust = 0.5),
    strip.text = element_text(face = "bold")
  ) +
  geom_text(
    data = stat_tests %>% filter(!is.na(p.label)),
    aes(x = Phase, y = ypos, label = p.label),
    inherit.aes = FALSE,
    position = position_dodge(width = 0.8),
    size = 4
  )

print(p)

ggsave(filename = out_plot, plot = p, width = 12, height = 5, dpi = 300)

###############################################################################
## 7) Optional: print test table
###############################################################################
cat("\nPairwise t-tests (KO vs WT) per Condition × Phase:\n")
print(stat_tests %>% arrange(Group, Condition, Phase))

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
###############################################################################
