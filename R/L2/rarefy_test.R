# Clear all existing data
rm(list=ls())

# Load packages
library(tidyverse)
library(lme4)
library(lmerTest)
library(janitor)
library(DHARMa)
library(bbmle)
library(ggeffects)
library(MuMIn)
library(vegan)

# Set working directory 
L0_dir <- Sys.getenv("L0DIR")
L1_dir <- Sys.getenv("L1DIR")
L2_dir <- Sys.getenv("L2DIR")
list.files(L2_dir)

# Read in CSV files
specabun <- read.csv(file.path(L2_dir, "species_abundance_SPEI.csv"))
plot <- read.csv(file.path(L2_dir, "plot_metrics.csv"))

specabun_meta <- left_join(specabun, plot)

specabun_wide_filter <- specabun_meta %>%
  rowid_to_column() %>% # add due to duplicates
  mutate(relative_abundance = as.integer(relative_abundance*10^6)) %>%
  filter(rowid < 201 | rowid > 272800) %>%
  pivot_wider(names_from = species,
              names_sep = "_",
              values_from = relative_abundance,
              values_fill = 0)


plot_area_vector <- specabun_wide_filter$measurement_scale_cover
plot_area_vector <- replace_na(plot_area_vector, 1)
plot_area_vector <- plot_area_vector*10

rarefied_data <- rarefy(specabun_wide_filter[20:144], sample = plot_area_vector)

data(BCI)
S <- specnumber(BCI) # observed number of species
(raremax <- min(rowSums(BCI)))
Srare <- rarefy(BCI, raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(BCI, step = 20, sample = raremax, col = "blue", cex = 0.6)
