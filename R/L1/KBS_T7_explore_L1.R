# TITLE:         LTER Grassland Rock: KBS - calculating diversity & stability & exploring data
# AUTHOR:        Ashley Darst
# COLLABORATORS: LTER synthesis group
# DATA INPUT:    Data imported as csvs from shared Google Drive L1 folder
# DATA OUTPUT:
# PROJECT:       LTER Grassland Rock synthesis group
# DATE:          11/15/2023

# Set-up ----
# Clear all existing data
rm(list=ls())

#Load packages
library(tidyverse)
library(janitor)

# Set working directory 
L0_dir <- Sys.getenv("L0DIR")
L1_dir <- Sys.getenv("L1DIR")
list.files(L1_dir)

# Load data ----
anpp <- read.csv(file.path(L1_dir, "KBS_MCSE_T7_ANPP.csv"))
spcomp <- read.csv(file.path(L1_dir, "KBS_MCSE_T7_SpComp.csv"))
anpp <- clean_names(anpp)
spcomp <- clean_names(spcomp)

# Clean data ----
# Remove unsorted data (493 unsorted!!)
spcomp <- spcomp %>%
  filter(species != "UnSorted")

# Summarize data ----
# Calculate ANPP at the plot level
# Don't really need this, less power
anpp_plot <- anpp %>%
  group_by(year, treatment, replicate, meantemp, annualprecip) %>%
  summarize(mean_anpp = mean(anpp_g_m2))

# Calculate richness
spcomp <- spcomp %>%
  group_by(year, treatment, station, replicate, disturbed_microplot, fertilized_microplot) %>%
  mutate(richness = n())

# Calculate Shannon diversity
spcomp <- spcomp %>%
  mutate(proportions = (pseudo_perc_cover/100)*(log(pseudo_perc_cover/100)))
spcomp$proportions[is.na(spcomp$proportions)] <- 0  # Fix Nas to zeros
spcomp <- spcomp %>%
  group_by(year, treatment, station, replicate, disturbed_microplot, fertilized_microplot) %>%
  mutate(shannon = -1*sum(proportions))

# Calculate stability of each sample across ALL years
stability <- spcomp %>%
  group_by(treatment, station, replicate, disturbed_microplot, fertilized_microplot) %>%
  summarize(mean_anpp = mean(anpp_total),
            sd_anpp = sd(anpp_total),
            mean_richness = mean(richness),
            mean_shannon = mean(shannon))
stability <- stability %>%
  mutate(stability = mean_anpp/sd_anpp)

# Calculate stability every 10 years (2019-2022 only 4 years)
# NOT all sites sampled every year (5 years in not enough to get SD)
stability10 <- spcomp %>%
  mutate(period = cut(year, breaks = c(1988, 1999, 2009, 2022), dig.lab = 4)) %>%
  group_by(treatment, station, replicate, disturbed_microplot, fertilized_microplot, period) %>%
  summarize(mean_anpp = mean(anpp_total),
            sd_anpp = sd(anpp_total),
            mean_richness = mean(richness),
            mean_shannon = mean(shannon))
stability10 <- stability10 %>%
  mutate(stability = mean_anpp/sd_anpp)

# Plot data ----
# Does mean temperature predict ANPP? (not really, maybe a peak at mid-temps for normal stations)
anpp %>%
  ggplot(aes(x = meantemp, y = anpp_g_m2, col = fertilized_microplot)) +
  geom_point() +
  geom_smooth()
anpp %>%
  ggplot(aes(x = meantemp, y = anpp_g_m2, col = fertilized_microplot)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~disturbed_microplot)

# Does annual precipitation predict ANPP? (increases with precip, especially normal stations)
# Stronger increasing trend in undisturbed than disturbed microplots)
anpp %>%
  ggplot(aes(x = annualprecip, y = anpp_g_m2, col = fertilized_microplot)) +
  geom_point() +
  geom_smooth(method = "lm")
anpp %>%
  ggplot(aes(x = annualprecip, y = anpp_g_m2, col = fertilized_microplot)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~disturbed_microplot)

# Does ANPP change over time? (started lower in stations)
anpp %>%
  ggplot(aes(x = year, y = anpp_g_m2, col = fertilized_microplot)) +
  geom_point() + 
  geom_smooth() +
  facet_wrap(~disturbed_microplot)

# Does richness change over time? (The microplots have a unimodel richness pattern, strange)
spcomp %>%
  group_by(year, treatment, station, replicate, disturbed_microplot, fertilized_microplot) %>%
  summarize(richness = mean(richness)) %>%
  ggplot(aes(x = year, y = richness, col = fertilized_microplot)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~disturbed_microplot, scales = "free_x")
spcomp %>%
  group_by(year, treatment, station, replicate, disturbed_microplot, fertilized_microplot) %>%
  summarize(richness = mean(richness)) %>%
  ggplot(aes(x = year, y = richness, col = replicate)) +
  geom_point() +
  facet_wrap(~disturbed_microplot, scales = "free_x")

# Do samples with greater richness have greater ANPP? (increase in undisturbed microplots and stations)
spcomp %>%
  group_by(year, treatment, station, replicate, disturbed_microplot, fertilized_microplot) %>%
  summarize(richness = mean(richness), anpp_total = mean(anpp_total)) %>%
  ggplot(aes(x = richness, y = anpp_total, col = fertilized_microplot)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~disturbed_microplot, scales = "free_x")

# Do samples with greater diversity (H) have greater ANPP? (increase in undisturbed microplots)
spcomp %>%
  group_by(year, treatment, station, replicate, disturbed_microplot, fertilized_microplot) %>%
  summarize(shannon = mean(shannon), anpp_total = mean(anpp_total)) %>%
  ggplot(aes(x = shannon, y = anpp_total, col = fertilized_microplot)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~disturbed_microplot, scales = "free_x")

# Are samples with greater richness more stable (u/o)? (No trend for stations, maybe increase in undisturbed microplots) 
# All years
stability %>%
  ggplot(aes(x = mean_richness, y = stability, col = fertilized_microplot)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~disturbed_microplot, scales = "free_x")
stability %>%
  ggplot(aes(x = mean_richness, y = stability)) +
  geom_point() +
  geom_smooth(method = "lm")

# Are samples with greater richness more stable (u/o)? (maybe increase in stations and disturbed microplots) 
# Every ~ 10 years (heavily impacted by which ten years you use)
stability10 %>%
  ggplot(aes(x = mean_richness, y = stability, col = fertilized_microplot)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~disturbed_microplot, scales = "free_x")
stability10 %>%
  ggplot(aes(x = mean_richness, y = stability)) +
  geom_point() +
  geom_smooth(method = "lm")

# Are samples with greater shannon more stable (u/o)? (No trend for stations, maybe increase in undisturbed microplots) 
# All years
stability %>%
  ggplot(aes(x = mean_richness, y = stability, col = fertilized_microplot)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~disturbed_microplot, scales = "free_x")

# Are samples with greater shannon more stable (u/o)? (maybe increase in stations) 
# Every ~ 10 years (heavily impacted by which ten years you use)
stability10 %>%
  ggplot(aes(x = mean_shannon, y = stability, col = fertilized_microplot)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~disturbed_microplot, scales = "free_x")

# Analysis ----
# Are richer samples more stable (u/o)? (No)
stability.lm <- lm(stability ~ fertilized_microplot + disturbed_microplot + mean_richness, data = stability)
summary(stability.lm)

# Are richer samples more stable (u/o)? (every ~10 years) (No)
stability10.lm <- lm(stability ~ fertilized_microplot + disturbed_microplot + mean_richness, data = stability10)
summary(stability10.lm)


