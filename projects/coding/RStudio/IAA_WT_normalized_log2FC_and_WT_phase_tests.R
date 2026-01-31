###############################################################################
## Targeted Metabolite Summary + Fold-Change Workflow (Trp + IAA)
## - Mean ± SE by Group
## - log2 fold-change vs phase-specific WT reference(s) (Log and Stationary)
## - Within-group Stationary vs Log contrast (|log2(S/L)|)
##
## Author: Lee Chua
## Last updated: 2026-01-31
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(ggsci)
})

###############################################################################
## 0) Configuration
###############################################################################
# You can keep this as ONE script and run both analytes (default),
# or set `analytes <- analytes["Trp"]` to run only one.

analytes <- list(
  Trp = list(
    input_file   = "MassSpec2.xlsx",
    sheet        = "Fulldata",
    value_col    = "Trp",
    plot_prefix  = "Trp",
    # WT references expected to exist in the data:
    wt_groups    = c("Wta", "Wtb"),
    # Group -> which WT to use as reference (edit as needed):
    reference_map = c(
      "ego3ko" = "Wta",
      "Wtb"    = "Wtb",
      "EGO3GSY" = "Wtb",
      "aro1ko"  = "Wtb",
      "ptr3ko"  = "Wtb",
      "Wta"     = "Wta"
    )
  ),
  IAA = list(
    input_file   = "MassSpec.xlsx",
    sheet        = "Fulldata",
    value_col    = "IAA",
    plot_prefix  = "IAA",
    wt_groups    = c("Wta", "Wtb"),
    reference_map = c(
      "ego3ko" = "Wta",
      "Wtb"    = "Wtb",
      "EGO3GSY" = "Wtb",
      "aro1ko"  = "Wtb",
      "ptr3ko"  = "Wtb",
      "Wta"     = "Wta"
    )
  )
)

# Output controls
out_dir_fig <- "figures"
dir.create(out_dir_fig, showWarnings = FALSE, recursive = TRUE)

###############################################################################
## 1) Helpers
###############################################################################
read_and_clean <- function(input_file, sheet, value_col) {
  # Expect 3 columns: Condition, Sample, {value_col}
  raw <- read_excel(input_file, sheet = sheet)
  
  if (ncol(raw) < 3) {
    stop("Expected at least 3 columns in the Excel sheet: Condition, Sample, Value.")
  }
  
  # Standardize to exactly 3 columns (first 3), and rename
  dat <- raw[, 1:3]
  colnames(dat) <- c("Condition", "Sample", value_col)
  
  dat %>%
    drop_na(Sample, .data[[value_col]]) %>%
    fill(Condition) %>%
    mutate(
      Group = str_extract(Sample, "^[^_]+"),
      Phase = str_extract(Sample, "_[LSls]") %>%
        str_remove("_") %>%
        toupper()
    ) %>%
    filter(!is.na(Group), Group != "Blank") %>%
    mutate(
      Phase = case_when(
        Phase == "L" ~ "Log",
        Phase == "S" ~ "Stationary",
        TRUE ~ Phase
      ),
      Phase = factor(Phase, levels = c("Log", "Stationary"))
    )
}

summarise_mean_se <- function(dat, value_col) {
  dat %>%
    group_by(Group) %>%
    summarise(
      mean = mean(.data[[value_col]], na.rm = TRUE),
      sd   = sd(.data[[value_col]], na.rm = TRUE),
      n    = sum(!is.na(.data[[value_col]])),
      se   = sd / sqrt(n),
      .groups = "drop"
    ) %>%
    arrange(desc(mean))
}

plot_mean_se <- function(summary_df, title, y_lab) {
  ggplot(summary_df, aes(x = reorder(Group, -mean), y = mean, fill = Group)) +
    geom_col(width = 0.6) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    labs(title = title, x = "Group", y = y_lab) +
    theme_minimal() +
    scale_fill_uchicago() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
}

compute_log2fc_vs_reference <- function(dat, value_col, phase_value, wt_groups, reference_map) {
  # Filter phase
  df_phase <- dat %>% filter(Phase == phase_value)
  
  # Compute WT means for the phase (reference values)
  ref_means <- df_phase %>%
    group_by(Group) %>%
    summarise(ref_mean = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop") %>%
    filter(Group %in% wt_groups) %>%
    deframe()
  
  # Map each Group to a WT label then to WT mean
  df_phase %>%
    mutate(
      ref_group = unname(reference_map[Group]),
      Reference = ref_means[ref_group]
    ) %>%
    filter(!is.na(Reference), Reference > 0, .data[[value_col]] > 0) %>%
    mutate(Log2FC = log2(.data[[value_col]] / Reference))
}

plot_log2fc_box <- function(df_fc, title) {
  ggplot(df_fc, aes(x = Group, y = Log2FC, fill = Group)) +
    geom_boxplot(outlier.shape = 16, outlier.size = 2) +
    labs(title = title, x = "Group", y = "log2(Fold Change)") +
    scale_fill_uchicago() +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
}

compute_abs_log2fc_stationary_vs_log <- function(dat, value_col) {
  # NOTE:
  # This uses a join by Group only, which creates all S×L pairwise combinations
  # within a group (cartesian product). If you have replicate IDs embedded in Sample,
  # you can refine this to true pairing (recommended when possible).
  log_df  <- dat %>% filter(Phase == "Log") %>% select(Group, !!sym(value_col))
  stat_df <- dat %>% filter(Phase == "Stationary") %>% select(Group, !!sym(value_col))
  
  colnames(log_df)[2]  <- "value_L"
  colnames(stat_df)[2] <- "value_S"
  
  inner_join(stat_df, log_df, by = "Group") %>%
    filter(value_S > 0, value_L > 0) %>%
    mutate(abs_Log2FC_S_vs_L = abs(log2(value_S / value_L)))
}

plot_abs_phase_fc <- function(df_phase_fc, title) {
  ggplot(df_phase_fc, aes(x = Group, y = abs_Log2FC_S_vs_L, fill = Group)) +
    geom_boxplot(outlier.shape = 16, outlier.size = 2) +
    labs(title = title, x = "Group", y = "|log2(Stationary / Log)|") +
    scale_fill_uchicago() +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
}

save_plot <- function(p, filename, width = 9, height = 5, dpi = 300) {
  ggsave(
    filename = file.path(out_dir_fig, filename),
    plot = p,
    width = width,
    height = height,
    dpi = dpi
  )
}

###############################################################################
## 2) Run workflow for each analyte
###############################################################################
for (analyte_name in names(analytes)) {
  
  cfg <- analytes[[analyte_name]]
  message("Running analyte: ", analyte_name)
  
  # 2.1 Read + clean
  dat <- read_and_clean(
    input_file = cfg$input_file,
    sheet      = cfg$sheet,
    value_col  = cfg$value_col
  )
  
  # 2.2 Mean ± SE
  sum_df <- summarise_mean_se(dat, cfg$value_col)
  p_mean <- plot_mean_se(
    sum_df,
    title = paste0("Mean ", cfg$plot_prefix, " by Group (± SE)"),
    y_lab = paste0(cfg$plot_prefix, " intensity (mean ± SE)")
  )
  print(p_mean)
  save_plot(p_mean, paste0(cfg$plot_prefix, "_mean_SE.png"))
  
  # 2.3 log2FC vs WT reference (Log)
  df_fc_log <- compute_log2fc_vs_reference(
    dat           = dat,
    value_col     = cfg$value_col,
    phase_value   = "Log",
    wt_groups     = cfg$wt_groups,
    reference_map = cfg$reference_map
  )
  p_fc_log <- plot_log2fc_box(
    df_fc_log,
    title = paste0(cfg$plot_prefix, ": Log phase log2 fold change vs reference")
  )
  print(p_fc_log)
  save_plot(p_fc_log, paste0(cfg$plot_prefix, "_log2FC_Log.png"))
  
  # 2.4 log2FC vs WT reference (Stationary)
  df_fc_stat <- compute_log2fc_vs_reference(
    dat           = dat,
    value_col     = cfg$value_col,
    phase_value   = "Stationary",
    wt_groups     = cfg$wt_groups,
    reference_map = cfg$reference_map
  )
  p_fc_stat <- plot_log2fc_box(
    df_fc_stat,
    title = paste0(cfg$plot_prefix, ": Stationary phase log2 fold change vs reference")
  )
  print(p_fc_stat)
  save_plot(p_fc_stat, paste0(cfg$plot_prefix, "_log2FC_Stationary.png"))
  
  # 2.5 |log2(S/L)| within-group phase change
  df_phase_fc <- compute_abs_log2fc_stationary_vs_log(dat, cfg$value_col)
  p_phase_fc <- plot_abs_phase_fc(
    df_phase_fc,
    title = paste0(cfg$plot_prefix, ": Stationary vs Log (|log2(S/L)|)")
  )
  print(p_phase_fc)
  save_plot(p_phase_fc, paste0(cfg$plot_prefix, "_absLog2FC_Stationary_vs_Log.png"))
}

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
#   4(43), 1686. https://doi.org/10.21105/joss.01686. :contentReference[oaicite:0]{index=0}
#
# Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag
#   New York. https://ggplot2.tidyverse.org. :contentReference[oaicite:1]{index=1}
#
# Wickham H, Bryan J (2025). readxl: Read Excel Files. R package version 1.4.5.
#   https://readxl.tidyverse.org. :contentReference[oaicite:2]{index=2}
#
# Xiao N (2025). ggsci: Scientific Journal and Sci-Fi Themed Color Palettes for
#   'ggplot2'. R package (see CRAN metadata). :contentReference[oaicite:3]{index=3}
