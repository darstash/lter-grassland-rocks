# TITLE:        LTER Grassland Rock: Create core figures 
# AUTHORS:      Ashley Darst
# COLLABORATORS:  
# DATA INPUT:   Data imported as csv files from shared Google drive L2 folder
# DATA OUTPUT:  Core figures
# PROJECT:      LTER Grassland Rock
# DATE:         October 2024

# Clear all existing data
rm(list=ls())

# Load packages
library(tidyverse)

# Set working directory 
L0_dir <- Sys.getenv("L0DIR")
L1_dir <- Sys.getenv("L1DIR")
L2_dir <- Sys.getenv("L2DIR")
list.files(L2_dir)

# Read in CSV files
plot <- read.csv(file.path(L2_dir, "plot_metrics_SPEI_diversity.csv"))
ece <- read.csv(file.path(L2_dir, "ece_resist_resil.csv"))

# Only keep distinct rows in ece
ece <- distinct(ece)

# Change ex_year to year
ece <- ece %>%
  rename(year = ex_year)

# Merge plot with resistance and resilience
plot_ece <- left_join(plot, ece)

# Make column with categories for high and low dominance
plot_ece <- plot_ece %>%
  mutate(dom_category = case_when(
    Berger_Parker >= 0.5 ~ "high",
    Berger_Parker < 0.5 ~ "low"
  ))

# Figure 1: resistance vs richness ----
# Why is resistance so large for some plots
plot_ece %>%
  ggplot(aes(x = Richness, y = resistance)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_y_log10() +
  annotation_logticks(sides = "l")

# Look at wet vs dry
plot_ece %>%
  drop_na(resistance) %>%
  ggplot(aes(x = Richness, y = resistance, col = spei6_category)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm") +
  scale_y_log10() +
  annotation_logticks(sides = "l")

plot_ece %>%
  drop_na(resistance) %>%
  ggplot(aes(x = log(resistance))) +
  geom_histogram() +
  facet_wrap(~site)
  
# Look at dominance (hard to see any pattern, low richness = high dominance)
plot_ece %>%
  drop_na(resistance) %>%
  ggplot(aes(x = Richness, y = resistance, col = Berger_Parker)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm", aes(linetype = dom_category)) +
  scale_y_log10() +
  annotation_logticks(sides = "l")

# Compare dominance vs richness
plot_ece %>%
  ggplot(aes(x = Richness, y = Berger_Parker)) +
  geom_point(alpha = 0.2) +
  geom_smooth()

# Figure 2: resilience vs richness ----
# Why are there some weird outliers
# No trend
plot_ece %>%
  ggplot(aes(x = Richness, y = resilience)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_y_log10() +
  annotation_logticks(sides = "l")

# Look at wet vs dry
plot_ece %>%
  drop_na(resistance) %>%
  ggplot(aes(x = Richness, y = resilience, col = spei6_category)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm") +
  scale_y_log10() +
  annotation_logticks(sides = "l")

# Look at dominance (hard to see any pattern, low richness = high dominance)
plot_ece %>%
  drop_na(resistance) %>%
  ggplot(aes(x = Richness, y = resilience, col = Berger_Parker)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm", aes(linetype = dom_category)) +
  scale_y_log10() +
  annotation_logticks(sides = "l")

