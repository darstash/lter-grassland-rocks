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
plot    <- read.csv(file.path(L1_dir, "plot_metrics_L1.csv"))      %>% 
  # select(! "X") %>%
  mutate(site                      = factor(site),
         original_measurement_unit = factor(original_measurement_unit) #,
         # dataset                   = factor(dataset),
         # source                    = factor(source)
         )

species <- read.csv(file.path(L1_dir, "species_abundance_L1.csv")) %>% 
  # select(! "X") %>%
  mutate(original_measurement_unit = factor(original_measurement_unit),
         cover_method              = factor(cover_method))

meta    <- read.csv(file.path(L1_dir, "metadata_L1.csv"))          %>% 
  # select(! "X") %>%
  mutate(site = factor(site),
         dataset = factor(dataset),
         treatment = factor(treatment),
         nutrients_added = factor(nutrients_added),
         disturbance = factor(disturbance),
         grazing = factor(grazing),
         source = factor(source),
         diversity_manipulated = factor(diversity_manipulated))

names(meta)
head(meta)
str(meta)
summary(meta)

# Metadata
## - growing_precipitation = growprecip (all site done!)
## - meantemp = temperature (all site done!)
## - annualprecip = precipitation (all site done!)
# - konza: some experiments are missing the nutrients added -> set NA, "" and "none" to "no_fertilizer" (KNZ individual: HQ), (done: CDR "" e247)
## - konza watershed (incorporated in higher order organization) -> delete from metadata (all site done)
## - get rid of experiment column (all site done)
# - grazing column: set "" to "ungrazed" (individual sites, KNZ HQ, Nutnet)
# - treatment: make sure control is always named control, CDR: add more treatment info back in -> make good metadata what which treatment variable means. (individual all, CDR and KNZ) CDR done
#
# - fire frequency and time since fire missing for some experiments (crucial? individual sites) -> what to do with 0 fire frequency when determining time since fire
# - Konza and kbs have NA in some disutbance as NA (crucial?, individual sites)

names(species)
head(species)
str(species)
summary(species %>% select(!c("cover_method", "area_sampled_bio", "area_sampled_cover")))

x <- species %>%
  group_by(year, site, dataset, plot, higher_order_organization, uniqueid, species, abundance) %>%
  summarize(n = n()) %>%
  filter(n > 1)

# - delete cover_method, area_sampled_cover, area_sampled_bio (all sites, )
## - original_measuurement_unit = original_measurement_unit (all sites, done)
## - some relative abundance values are scaled to 100 instead of 1 in KBS and KNZ (individual all, done)
## - KNZ WAT -> add KNZ to site column!!! (KNZ individual, done)

# species names
#_______________
## - knz nutnet ._. ????? (individual knz, done)
# - sites: think about filtering categories (and export them for metadata) like CDR
#     genus_sp_in_biomass (list of things that are identified to genus level: genus sp.)
#     non_plant_things_in_biomass
#     maybe_plant_things_in_biomass (miscellaneous grasses etc.)
# - unify species names -> capitalization, underscores, etc (all site)

#biomass was collected in 2020 for KNZ PPLot but not spp comp due to covid -> leave in, add comment?
#suggest removing ConsME data since we dont have up to five years of data (all data sites) -> inclusion criteria not fulfilled, remove!

names(plot)
head(plot)
str(plot)
summary(plot)

plot %>% 
  group_by(site, higher_order_organization, uniqueid) %>% 
  summarize(duration = max(year) - min(year)) %>%
  filter(duration < 5) %>%
  View()


# - KNZ & KBS: issues with year in uniqueid
## - remove metadatastuff  from plotlevel data: treatment, source, dataset (done)
# - remove speciesstuff from plotlevel data: orignial_measurement_unit
## - things that should be calculated based on the full species dataset (add to plotlevel data LATER -> remove now): evenness, shannon, richness (done)
# - for all sites in some experiments the measurement_scale_cover is missing! KNZ HQ nutnet wat, KBS glbrc, done: cdr e061 e245 e054 e001
# - Konza and KBS: in some experiments the measurement_scale_biomass is missing! KBS glbrc, KNZ HQ
## - Konza: 1NA in the biomass KNZ_VI_2c_3_22_2 2022 (OK)

## ALL DATASETS - Konza: fix year issue-Ramps (done)













