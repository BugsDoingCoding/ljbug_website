###############################################################################
## Yeast Growth Curves: Endpoint + Midpoint OD Summaries (Bar Plots)
##
## - Reads a wide OD timecourse CSV (columns like X0, X2, X24 ...)
## - Converts to long format
## - Produces:
##   (1) Endpoint mean ± SD at a specified timepoint (default: 24)
##   (2) Midpoint mean ± SD (closest observed time to max_time/2 per strain)
##
## Author: Lee Chua
## Last updated: 2026-01-31
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
})

###############################################################################
## 0) Configuration
###############################################################################
input_csv <- "Fake_data_set.csv"
endpoint_time <- 24  # change as needed

###############################################################################
## 1) Load and reshape data
###############################################################################
if (!file.exists(input_csv)) stop("File not found: ", input_csv)

data_wide <- read.csv(input_csv, check.names = FALSE, stringsAsFactors = FALSE)

# Ensure the first column is the strain identifier
colnames(data_wide)[1] <- "Yeast_Strain"

data_long <- data_wide %>%
  pivot_longer(
    cols = -Yeast_Strain,
    names_to = "Time",
    values_to = "OD",
    names_repair = "unique"
  ) %>%
  mutate(
    Time = as.numeric(gsub("^X", "", Time)),
    OD = as.numeric(OD)
  ) %>%
  filter(!is.na(Time), !is.na(OD))

if (nrow(data_long) == 0) stop("No valid (Time, OD) observations after cleaning.")

###############################################################################
## 2) Endpoint summary (mean ± SD at endpoint_time)
###############################################################################
endpoint_summary <- data_long %>%
  filter(Time == endpoint_time) %>%
  group_by(Yeast_Strain) %>%
  summarise(
    mean_OD = mean(OD, na.rm = TRUE),
    sd_OD   = sd(OD, na.rm = TRUE),
    n       = sum(!is.na(OD)),
    .groups = "drop"
  )

if (nrow(endpoint_summary) == 0) {
  warning("No observations found at endpoint_time = ", endpoint_time, ". Endpoint plot will be empty.")
}

p_endpoint <- ggplot(endpoint_summary, aes(x = Yeast_Strain, y = mean_OD, fill = Yeast_Strain)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_OD - sd_OD, ymax = mean_OD + sd_OD), width = 0.2) +
  labs(
    title = paste0("Mean Endpoint ODs (Time = ", endpoint_time, ")"),
    x = "Yeast Strain",
    y = "Mean OD (± SD)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

print(p_endpoint)

###############################################################################
## 3) Midpoint summary (mean ± SD at closest observed time to max_time/2)
###############################################################################
max_time_per_strain <- data_long %>%
  group_by(Yeast_Strain) %>%
  summarise(
    max_time = max(Time, na.rm = TRUE),
    mid_time = max_time / 2,
    .groups = "drop"
  )

midpoint_data <- data_long %>%
  left_join(max_time_per_strain, by = "Yeast_Strain") %>%
  group_by(Yeast_Strain) %>%
  filter(abs(Time - mid_time) == min(abs(Time - mid_time), na.rm = TRUE)) %>%
  ungroup()

midpoint_summary <- midpoint_data %>%
  group_by(Yeast_Strain) %>%
  summarise(
    mean_OD = mean(OD, na.rm = TRUE),
    sd_OD   = sd(OD, na.rm = TRUE),
    n       = sum(!is.na(OD)),
    .groups = "drop"
  )

p_midpoint <- ggplot(midpoint_summary, aes(x = Yeast_Strain, y = mean_OD, fill = Yeast_Strain)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_OD - sd_OD, ymax = mean_OD + sd_OD), width = 0.2) +
  labs(
    title = "Mean Midpoint ODs (Closest Observed Time to max/2)",
    x = "Yeast Strain",
    y = "Mean OD (± SD)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

print(p_midpoint)

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
###############################################################################
