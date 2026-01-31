# Load necessary libraries
library(tidyverse)

# Read the CSV file
data <- read.csv("Fake_data_set.csv")

# Ensure the column names are unique
names(data)[1] <- "Yeast_Strain"

# Convert the data to long format for easier manipulation
data_long <- pivot_longer(data, cols = -Yeast_Strain, 
                          names_to = "Time", values_to = "OD", 
                          names_repair = "unique")

# Convert Time to numeric
data_long$Time <- as.numeric(gsub("X", "", data_long$Time))

# Filter out time points with no data
data_long <- data_long %>% filter(!is.na(OD))

# Determine the maximum Time (endpoint) for each strain
max_time_per_strain <- data_long %>%
  group_by(Yeast_Strain) %>%
  summarize(max_time = max(Time, na.rm = TRUE))

# Calculate the midpoint time for each strain
midpoints <- max_time_per_strain %>%
  mutate(mid_time = max_time / 2)

# Add the midpoint time information to the data
data_long <- data_long %>%
  left_join(midpoints, by = "Yeast_Strain")

# For each strain, find the OD closest to the calculated midpoint
midpoint_data <- data_long %>%
  group_by(Yeast_Strain) %>%
  filter(abs(Time - mid_time) == min(abs(Time - mid_time))) %>%
  ungroup()

# Calculate mean and standard deviation of the ODs at the midpoint for each strain
midpoint_summary <- midpoint_data %>%
  group_by(Yeast_Strain) %>%
  summarize(mean_OD = mean(OD), sd_OD = sd(OD))

# Create a bar graph comparing mean midpoint ODs for all strains
ggplot(midpoint_summary, aes(x = Yeast_Strain, y = mean_OD, fill = Yeast_Strain)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_OD - sd_OD, ymax = mean_OD + sd_OD), width = 0.2) +
  labs(
    title = "Mean Midpoint ODs for All Yeast Strains",
    x = "Yeast Strain",
    y = "Mean OD (± SD)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

