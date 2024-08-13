# TITLE:        LTER Grassland Rock: Merge CDR, KBS, and KNZ datasets
# AUTHORS:      
# COLLABORATORS:  
# DATA INPUT:   Data imported as cleaned site-specific csv files from shared Google drive L1 folder
# DATA OUTPUT:    
# PROJECT:      LTER Grassland Rock
# DATE:         July 2024

# this code is to harmonize plant biomass and plant composition datasets from CDR, KBS, and KNZ

# Clear all existing data
rm(list=ls())

# Load packages
library(tidyverse)
library(gtools)
library(lubridate)

# Set working directory 
L0_dir <- Sys.getenv("L0DIR")
L1_dir <- Sys.getenv("L1DIR")
L2_dir <- Sys.getenv("L2DIR")
list.files(L1_dir)

# Read in CSV files ----
# CDR
CDR_metadata <- read.csv(file.path(L1_dir, "CDR_metadata.csv"))
CDR_species_abundance <- read.csv(file.path(L1_dir, "CDR_specieslevel_abundance.csv"))
CDR_plot_metrics <- read.csv(file.path(L1_dir, "CDR_plotlevel_metrics.csv"))

# KBS
KBS_metadata <- read.csv(file.path(L1_dir, "KBS_metadata.csv"))
KBS_species_abundance <- read.csv(file.path(L1_dir, "KBS_species_level_abundance.csv"))
KBS_plot_metrics <- read.csv(file.path(L1_dir, "KBS_plot_level_metrics.csv"))

# KNZ 
KNZ_metadata <- read.csv(file.path(L1_dir, "KNZ_metadata.csv"))
KNZ_species_abundance <- read.csv(file.path(L1_dir, "KNZ_specieslevel_abundance.csv"))
KNZ_plot_metrics <- read.csv(file.path(L1_dir, "KNZ_plotlevel_metrics.csv"))

# Merge datasets ----
# Plot metrics

# Make plot metric datasets compatible for merge
CDR_plot_metrics$measurement_scale_biomass <- gsub("\\m.*", "", CDR_plot_metrics$measurement_scale_biomass)
CDR_plot_metrics$measurement_scale_biomass <- as.double(CDR_plot_metrics$measurement_scale_biomass)
CDR_plot_metrics$measurement_scale_cover <- gsub("\\m.*", "", CDR_plot_metrics$measurement_scale_cover)
CDR_plot_metrics$measurement_scale_cover <- as.double(CDR_plot_metrics$measurement_scale_cover)
KNZ_plot_metrics <- rename(KNZ_plot_metrics, shannon = Shannon)
CDR_plot_metrics$plot <- as.character(CDR_plot_metrics$plot)
KNZ_plot_metrics$plot <- as.character(KNZ_plot_metrics$plot)

plot_metrics <- full_join(CDR_plot_metrics, KBS_plot_metrics)
plot_metrics <- full_join(plot_metrics, KNZ_plot_metrics)

# write.csv(plot_metrics, file.path(L2_dir, "./plot_metrics.csv"), row.names=F)

# Metadata

# Make metadata datasets compatible for merge
CDR_metadata$plot <- as.character(CDR_metadata$plot)
KNZ_metadata$plot <- as.character(KNZ_metadata$plot)

metadata <- full_join(CDR_metadata, KBS_metadata)
metadata <- full_join(metadata, KNZ_metadata) # Two different nutrient added columns for KNZ (what is happening?)

# write.csv(metadata, file.path(L2_dir, "./metadata.csv"), row.names=F)

# Species abundance

# Make species abundance datasets compatible for merge
CDR_species_abundance$plot <- as.character(CDR_species_abundance$plot)
KNZ_species_abundance$plot <- as.character(KNZ_species_abundance$plot)

species_abundance <- full_join(CDR_species_abundance, KBS_species_abundance)
species_abundance <- full_join(species_abundance, KNZ_species_abundance)

# write.csv(species_abundance, file.path(L2_dir, "./species_abundance.csv"), row.names=F)

