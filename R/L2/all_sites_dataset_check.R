# TITLE:        LTER Grassland Rock: Data checking
# AUTHORS:      Ashley Darst, Joshua Ajowele, Seraina Cappelli
# COLLABORATORS:  
# DATA INPUT:   Data imported as cleaned site-specific csv files from shared Google drive L2 folder
# DATA OUTPUT:    
# PROJECT:      LTER Grassland Rock
# DATE:         October 2024 (Working group meeting)

# this code is to check and validate the combined data sets, before approving the use of them for analysis
# datasets to check:
# - species abundance data
# - plot metric data
# - metadata

# Clear all existing data
rm(list=ls())

# Load packages
library(tidyverse)

# Set working directory 
L0_dir <- Sys.getenv("L0DIR")
L1_dir <- Sys.getenv("L1DIR")
L2_dir <- Sys.getenv("L2DIR")
list.files(L2_dir)

# Get data
plot    <- read.csv(file.path(L2_dir, "plot_metrics.csv"))      %>% select(! "X")
species <- read.csv(file.path(L2_dir, "species_abundance.csv")) %>% select(! "X")
meta    <- read.csv(file.path(L2_dir, "metadata.csv"))          %>% select(! "X")

names(meta)
head(meta)
str(meta)
summary(meta)

# Metadata
# - growing_precipitation = growprecip
# - meantemp = temperature 
# - annualprecip = precipitation
# - konza: some experiments are missing the nutrients added -> set NA, "" and "none" to "no_fertilizer",
# - konza watershed (incorporated in higher order organization) -> delete from metadata
# - get rid of experiment column
# - grazing column: set "" to "ungrazed"
#
# - fire frequency and time since fire missing for some experiments (crucial?)
# - Konza and kbs have NA in some disutbance as NA (crucial?)
# - growtemp missing for konza (not crucial)


head(species)
str(species)
summary(species %>% select(!c("cover_method", "area_sampled_bio", "area_sampled_cover")))

# delete cover_method, area_sampled_bio, area_sampled_cover
# original_measuurement_unit = original_measurement_unit
# some relative abundance values are scaled to 100 instead of 1 in KBS and KNZ
# KNZ WAT -> add KNZ to site column!!!

# species names
#_______________
# knz nutnet ._. ?????
# sites: think about filtering categories (and export them for metadata)
#   genus_sp_in_biomass (list of things that are identified to genus level: genus sp.)
#   non_plant_things_in_biomass
#   maybe_plant_things_in_biomass (miscellaneous grasses etc.)
# unify species names -> capitalization, underscores, etc

head(plot)
str(plot)
summary(plot)



# ALL DATASETS - Konza: fix year issue






# things that should have happened in earier code
# - remove metadatastuff  from plotlevel data: treatment, source
# - remove speciesstuff from plotlevel data: orignial_measurement_unit
# - things that should be calculated based on the full species dataset (add to plotlevel data later): evenness, shannon, richness (?)





# plot level data: for all sites in some experiments the measurement_scale_cover is missing (not crucial)
# plot level data: Konza and KBS: in some experiments the measurement_scale_biomass is missing (not crucial)

