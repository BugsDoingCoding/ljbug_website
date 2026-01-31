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

# Calculate mean and standard deviation for each strain at each time point
data_summary <- data_long %>%
  group_by(Yeast_Strain, Time) %>%
  summarize(mean_OD = mean(OD), sd_OD = sd(OD))

# Filter data for the endpoint (Time = 24)
endpoint_summary <- data_summary %>% filter(Time == 24)

# Create a bar graph comparing mean endpoint ODs for all strains
ggplot(endpoint_summary, aes(x = Yeast_Strain, y = mean_OD, fill = Yeast_Strain)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_OD - sd_OD, ymax = mean_OD + sd_OD), width = 0.2) +
  labs(
    title = "Mean Endpoint ODs for All Yeast Strains",
    x = "Yeast Strain",
    y = "Mean OD (± SD)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

