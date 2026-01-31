# Load necessary libraries
library(tidyverse)
library(ggplot2)
library(readxl)
library(ggsci)

# Read the CSV file
data <- read_excel("YMD_GC.xlsx")

# Ensure the column names are unique
names(data)[1] <- "Yeast_Strain"

# Convert the data to long format for easier manipulation
data_long <- pivot_longer(data, cols = -Yeast_Strain, 
                          names_to = "Time", 
                          values_to = "OD", 
                          names_repair = "unique")

# Convert Time to numeric
data_long$Time <- as.numeric(gsub("X", "", data_long$Time))

# Filter out time points with no data
data_long <- data_long %>% filter(!is.na(OD))

# Calculate mean and standard deviation for each strain at each time point
data_summary <- data_long %>%
  group_by(Yeast_Strain, Time) %>%
  summarize(mean_OD = mean(OD), sd_OD = sd(OD)/sqrt(3))

# Plot the time series
ggplot(data_summary, aes(x = Time, 
                         y = mean_OD, 
                         color = Yeast_Strain, 
                         group = Yeast_Strain)) +
  geom_line(size = 1.5) + # Make lines thicker
  geom_ribbon(aes(ymin = mean_OD - sd_OD, 
                  ymax = mean_OD + sd_OD, 
                  fill = Yeast_Strain),
                  alpha = 0.2) +
                  labs(title = "RM11 YMD+WYFM Growth Curve", #Name of Growth Curve Chart
                  x = "Time (hours)", 
                  y = "OD") +
  theme_minimal(base_size = 16)+
  scale_fill_tron()

# Calculate relative growth rate
data_long <- data_long %>%
  group_by(Yeast_Strain) %>%
  mutate(relative_growth_rate = (OD - lag(OD)) / lag(OD))

# Summarize the relative growth rate
growth_summary <- data_long %>%
  group_by(Yeast_Strain, Time) %>%
  summarize(mean_growth_rate = mean(relative_growth_rate, na.rm = TRUE),
            sd_growth_rate = sd(relative_growth_rate, na.rm = TRUE))

# Filter out the first 8 time points in growth_summary
growth_summary <- growth_summary %>%
  filter(Time > sort(unique(Time))[1])

# Plot the relative growth rate
ggplot(growth_summary, aes(x = Time,
                           y = mean_growth_rate, 
                           color = Yeast_Strain, 
                           group = Yeast_Strain)) +
  geom_line(size = 1.5) + # Make lines thicker
  geom_ribbon(aes(ymin = mean_growth_rate - sd_growth_rate,
                  ymax = mean_growth_rate + sd_growth_rate, 
                  fill = Yeast_Strain), 
                  alpha = 0.2) +
  labs(title = "RM11 ptr3∆ YPD vs Rapamycin Relative Growth Rates", 
      x = "Time (hours)", 
      y = "Relative Growth Rate") +
      theme_minimal(base_size = 15) # Increase base font size

