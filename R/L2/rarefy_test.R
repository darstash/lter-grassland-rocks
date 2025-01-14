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

specabun_wide <- specabun_meta %>%
  rowid_to_column() %>% # add due to duplicates
  pivot_wider(names_from = species,
              names_sep = "_",
              values_from = relative_abundance,
              values_fill = list(Value = 0))
specabun_wide <- specabun_wide %>% mutate_all( ~replace(., lengths(.)==0, 0))
specabun_wide[19:1475] <- sapply(specabun_wide[19:1475],as.numeric)

plot_area_vector <- specabun_wide$measurement_scale_cover

rarefied_data <- rarefy(specabun_wide[19:1475], sample = plot_area_vector)
