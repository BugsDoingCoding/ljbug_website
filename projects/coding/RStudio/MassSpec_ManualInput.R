# === Load packages ===
library(tidyverse)
library(readxl)
library(ggsci)

# === 1. Load and prepare data ===
fulldata <- read_excel("MassSpec.xlsx", sheet = "Fulldata")
colnames(fulldata) <- c("Condition", "Sample", "IAA")

fulldata <- fulldata %>%
  drop_na(Sample, IAA) %>%
  fill(Condition) %>%
  mutate(
    Group = str_extract(Sample, "^[^_]+"),
    Phase = str_extract(Sample, "_[LSls]") %>% str_remove("_") %>% toupper()
  ) %>%
  filter(Group != "Blank")

# === 2. Mean plot with SE error bars ===
summary_stats <- fulldata %>%
  group_by(Group) %>%
  summarise(
    Mean = mean(IAA),
    SD = sd(IAA),
    N = n(),
    SE = SD / sqrt(N),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean))

ggplot(summary_stats, aes(x = reorder(Group, -Mean), y = Mean, fill = Group)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  labs(
    title = "Mean IAA by Group with Standard Error",
    x = "Group",
    y = "IAA Intensity"
  ) +
  theme_minimal() +
  scale_fill_uchicago() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

# === 3. Log2 Fold Change Plot (Log Phase) ===
log_data <- fulldata %>% filter(Phase == "L")

ref_vals_log <- log_data %>%
  group_by(Group) %>%
  summarise(Mean = mean(IAA), .groups = "drop") %>%
  filter(Group %in% c("Wta", "Wtb")) %>%
  deframe()

log_data <- log_data %>%
  mutate(
    Reference = case_when(
      Group == "ego3ko" ~ ref_vals_log["Wta"],
      Group %in% c("Wtb", "EGO3GSY", "aro1ko", "ptr3ko") ~ ref_vals_log["Wtb"],
      Group == "Wta" ~ ref_vals_log["Wta"],
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(Reference) & Reference > 0 & IAA > 0) %>%
  mutate(Log2FC = log2(IAA / Reference))

ggplot(log_data, aes(x = Group, y = Log2FC, fill = Group)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 2) +
  labs(title = "Plot 1: Log Phase log2 Fold Changes", 
       x = "Strain", 
       y = "log2(FC)") +
  scale_fill_uchicago() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

# === 4. Log2 Fold Change Plot (Stationary Phase) ===
stationary_data <- fulldata %>% filter(Phase == "S")

ref_vals_stat <- stationary_data %>%
  group_by(Group) %>%
  summarise(Mean = mean(IAA), .groups = "drop") %>%
  filter(Group %in% c("Wta", "Wtb")) %>%
  deframe()

stationary_data <- stationary_data %>%
  mutate(
    Reference = case_when(
      Group == "ego3ko" ~ ref_vals_stat["Wta"],
      Group %in% c("Wtb", "EGO3GSY", "aro1ko", "ptr3ko") ~ ref_vals_stat["Wtb"],
      Group == "Wta" ~ ref_vals_stat["Wta"],
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(Reference) & Reference > 0 & IAA > 0) %>%
  mutate(Log2FC = log2(IAA / Reference))

ggplot(stationary_data, aes(x = Group, y = Log2FC, fill = Group)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 2) +
  labs(title = "Plot 2: Stationary Phase log2 Fold Changes", 
       x = "Strain", 
       y = "log2(FC)") +
  scale_fill_jco() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

# === 5. S vs L Fold Change Plot ===

log_data <- fulldata %>% filter(Phase == "L")
stat_data <- fulldata %>% filter(Phase == "S")

# Inner join by Group only — not Sample
paired_data <- inner_join(stat_data, log_data, by = "Group", suffix = c("_S", "_L")) %>%
  mutate(Log2FC_S_vs_L = abs(log2(IAA_S / IAA_L)))  # Or just log2(...) for signed

# Now plot
ggplot(paired_data, aes(x = Group, y = Log2FC_S_vs_L, fill = Group)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 2) +
  labs(
    title = "Stationary vs Log Phase",
    x = "Group",
    y = "|log2(Stationary / Log)|"
  ) +
  scale_fill_jco() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

