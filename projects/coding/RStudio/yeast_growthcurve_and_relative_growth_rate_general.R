###############################################################################
## General Growth Curve + Relative Growth Rate Workflow (CSV or Excel)
## - Accepts a wide-format table: first column = strain/group ID, remaining columns
##   = timepoints (e.g., X0, X2, X24) with OD measurements (replicates as rows)
## - Produces:
##   (1) Mean OD over time with SE ribbon
##   (2) Mean relative growth rate over time with SD ribbon
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
# ---- Input ----
# Works for .csv, .xlsx, .xls
input_file  <- "YMD_GC.xlsx"
excel_sheet <- 1          # ignored for CSV; for Excel can be index or sheet name

# Name to assign to the first column (group/strain ID)
id_col_name <- "Yeast_Strain"

# Plot titles (edit as desired)
title_growth_curve <- "Growth Curve"
title_growth_rate  <- "Relative Growth Rate"

# Summary behavior
# - Use SE ribbon for OD means (typical)
# - Use SD ribbon for relative growth rates (typical)
od_ribbon <- c("se", "sd")[1]          # "se" or "sd"
min_points_for_sd <- 2                # guardrail for sd/se calculations

# Output (optional): set to TRUE to save figures
save_figures <- FALSE
out_dir_fig  <- "figures"
dir.create(out_dir_fig, showWarnings = FALSE, recursive = TRUE)

###############################################################################
## 1) Helpers
###############################################################################
stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
}

read_wide_table <- function(path, sheet = 1) {
  stop_if_missing(path)
  
  ext <- tolower(tools::file_ext(path))
  
  if (ext %in% c("xlsx", "xls")) {
    df <- readxl::read_excel(path, sheet = sheet)
    return(as.data.frame(df))
  }
  
  if (ext == "csv") {
    df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    return(df)
  }
  
  stop("Unsupported file type: .", ext, " (supported: csv, xlsx, xls)")
}

parse_time_from_colname <- function(x) {
  # Handles "X0", "X2", "0", "2", "Time_24" (extracts first number-like token)
  # Returns numeric vector (NA if no number found).
  out <- suppressWarnings(as.numeric(str_extract(x, "-?\\d+\\.?\\d*")))
  out
}

compute_sd <- function(x) {
  if (sum(!is.na(x)) < min_points_for_sd) return(NA_real_)
  sd(x, na.rm = TRUE)
}

compute_se <- function(x) {
  n <- sum(!is.na(x))
  if (n < min_points_for_sd) return(NA_real_)
  sd(x, na.rm = TRUE) / sqrt(n)
}

###############################################################################
## 2) Load and reshape to long format
###############################################################################
wide <- read_wide_table(input_file, sheet = excel_sheet)

if (ncol(wide) < 2) stop("Input must have at least 2 columns: ID + >=1 timepoint.")

# Standardize ID column name
colnames(wide)[1] <- id_col_name

# Pivot longer: each row becomes one observation at one timepoint
long <- wide %>%
  pivot_longer(
    cols = -all_of(id_col_name),
    names_to = "Time_raw",
    values_to = "OD",
    names_repair = "unique"
  ) %>%
  mutate(
    Time = parse_time_from_colname(Time_raw),
    OD   = as.numeric(OD)
  ) %>%
  filter(!is.na(Time), !is.na(OD))

if (nrow(long) == 0) stop("No valid (Time, OD) observations after cleaning/time parsing.")

# Ensure time ordering is sensible
long <- long %>% arrange(.data[[id_col_name]], Time)

###############################################################################
## 3) Growth curve summary: mean OD ± (SE or SD)
###############################################################################
od_summary <- long %>%
  group_by(.data[[id_col_name]], Time) %>%
  summarise(
    mean_OD = mean(OD, na.rm = TRUE),
    sd_OD   = compute_sd(OD),
    se_OD   = compute_se(OD),
    n       = sum(!is.na(OD)),
    .groups = "drop"
  ) %>%
  mutate(
    ribbon = if_else(od_ribbon == "se", se_OD, sd_OD)
  )

p_growth <- ggplot(
  od_summary,
  aes(x = Time, y = mean_OD, color = .data[[id_col_name]], group = .data[[id_col_name]])
) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(
    aes(
      ymin = mean_OD - ribbon,
      ymax = mean_OD + ribbon,
      fill = .data[[id_col_name]]
    ),
    alpha = 0.2,
    color = NA
  ) +
  labs(
    title = title_growth_curve,
    x = "Time (hours)",
    y = "OD"
  ) +
  theme_minimal(base_size = 16) +
  scale_fill_tron()

print(p_growth)

###############################################################################
## 4) Relative growth rate: (OD_t - OD_{t-1}) / OD_{t-1}
###############################################################################
# Compute per replicate row (i.e., per original row in wide table)
# NOTE: This assumes each row in the wide file is a replicate timecourse for a strain.
# The growth rate is computed within each replicate trajectory.
# If your file structure differs (e.g., replicates encoded differently), tell me.

# Add a replicate identifier based on original row number within each strain.
# This helps compute lag() within a replicate rather than across replicates.
long_with_rep <- long %>%
  group_by(.data[[id_col_name]]) %>%
  mutate(.replicate_id = row_number()) %>%
  ungroup()

# The above replicate logic is a best-effort default; if your input has true replicate IDs,
# swap .replicate_id to that column.

growth_long <- long_with_rep %>%
  group_by(.data[[id_col_name]], .replicate_id) %>%
  arrange(Time, .by_group = TRUE) %>%
  mutate(
    relative_growth_rate = (OD - lag(OD)) / lag(OD)
  ) %>%
  ungroup() %>%
  filter(is.finite(relative_growth_rate))

growth_summary <- growth_long %>%
  group_by(.data[[id_col_name]], Time) %>%
  summarise(
    mean_growth_rate = mean(relative_growth_rate, na.rm = TRUE),
    sd_growth_rate   = compute_sd(relative_growth_rate),
    n                = sum(!is.na(relative_growth_rate)),
    .groups = "drop"
  ) %>%
  # Drop the earliest timepoint per strain (lag undefined)
  group_by(.data[[id_col_name]]) %>%
  filter(Time != min(Time, na.rm = TRUE)) %>%
  ungroup()

p_growth_rate <- ggplot(
  growth_summary,
  aes(x = Time, y = mean_growth_rate, color = .data[[id_col_name]], group = .data[[id_col_name]])
) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(
    aes(
      ymin = mean_growth_rate - sd_growth_rate,
      ymax = mean_growth_rate + sd_growth_rate,
      fill = .data[[id_col_name]]
    ),
    alpha = 0.2,
    color = NA
  ) +
  labs(
    title = title_growth_rate,
    x = "Time (hours)",
    y = "Relative Growth Rate"
  ) +
  theme_minimal(base_size = 15) +
  scale_fill_tron()

print(p_growth_rate)

###############################################################################
## 5) Optional: save figures
###############################################################################
if (isTRUE(save_figures)) {
  ggsave(file.path(out_dir_fig, "growth_curve.png"), p_growth, width = 10, height = 6, dpi = 300)
  ggsave(file.path(out_dir_fig, "relative_growth_rate.png"), p_growth_rate, width = 10, height = 6, dpi = 300)
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
###############################################################################
