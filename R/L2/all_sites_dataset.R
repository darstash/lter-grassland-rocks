# TITLE:        LTER Grassland Rock: Merge CDR, KBS, and KNZ datasets
# AUTHORS:      Ashley Darst
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

plot_metrics <- plot_metrics %>%
  select(-c(source, shannon, evenness, richness, X, dataset))

# write.csv(plot_metrics, file.path(L2_dir, "./plot_metrics.csv"), row.names=F)

# Metadata

# Make metadata datasets compatible for merge
CDR_metadata$plot <- as.character(CDR_metadata$plot)
KNZ_metadata$plot <- as.character(KNZ_metadata$plot)

metadata <- full_join(CDR_metadata, KBS_metadata)
metadata <- full_join(metadata, KNZ_metadata) # Two different nutrient added columns for KNZ

# Merge temperature columns
metadata <- metadata %>%
  mutate(meantemp = coalesce(meantemp, temperature)) %>%
  select(-temperature)
metadata <- metadata %>%
  mutate(annualprecip = coalesce(annualprecip, precipitation)) %>%
  select(-c(precipitation, growtemp, experiment))
metadata <- metadata %>%
  mutate(growprecip = coalesce(growprecip, growing_precipitation)) %>%
  select(-c(growing_precipitation, X))
  
# write.csv(metadata, file.path(L2_dir, "./metadata.csv"), row.names=F)

# Species abundance

# Make species abundance datasets compatible for merge
CDR_species_abundance$plot <- as.character(CDR_species_abundance$plot)
KNZ_species_abundance$plot <- as.character(KNZ_species_abundance$plot)

species_abundance <- full_join(CDR_species_abundance, KBS_species_abundance)
species_abundance <- full_join(species_abundance, KNZ_species_abundance)

species_abundance <- species_abundance %>%
  select(-X)

######### UNFINISHED SPECIES NAMES FIXING #####################
# make species information dataset
# idea: at this stage, don't touch names in the dataset anymore, but have a 
# dataframe that assigns each weirdo and non-weirdo entry in the species column
# a cleaned up version of names 

## CDR list of things
genus_sp_in_biomass <- c("Acer sp.",
                         "Allium sp.",
                         "Alnus sp.",
                         "Apocynum sp.",
                         "Arabis sp.",
                         "Aristida sp.",
                         "Asclepias sp.",
                         "Aster sp.",
                         "Bromus sp.",
                         "Calamovilfa sp.",
                         "Carex sp.",
                         "Chenopodium sp.",
                         "Cirsium sp.",
                         "Cyperus sp.",
                         "Digitaria sp.",
                         "Digitaria sp.",
                         "Equisetum sp.",
                         "Eragrostis sp.",
                         "Erigeron sp.",
                         "Galium sp.",
                         "Helianthus sp.",
                         "Hieracium sp.",
                         "Juncus sp.",
                         "Lactuca sp.",
                         "Liatris sp.",
                         "Lithospermum sp.",
                         "Melilotus sp.",
                         "Oenothera sp.",
                         "Oxalis sp.",
                         "Panicum sp.",
                         "Parthenocissus sp.",
                         "Penstemon sp.",
                         "Pinus sp.",
                         "Plantago sp.",
                         "Poa sp.",
                         "Polygala sp.",
                         "Polygonatum sp.",
                         "Potentilla sp.",
                         "Prunus sp.",
                         "Quercus borealis-ellipsoidalis",
                         "Quercus sp.",
                         "Ranunculus sp.",
                         "Rhus sp.",
                         "Rudbeckia sp.",
                         "Rubus sp.",
                         "Rumex sp.",
                         "Salix sp.",
                         "Sedges",
                         "Senecio sp.",
                         "Setaria sp.",
                         "Silene sp.",
                         "Solidago sp.",
                         "Sporobolus sp.",
                         "Tradescantia sp.",
                         "Tragopogon sp.",
                         "Trifolium sp.",
                         "Viola sp.")

non_plant_things_in_biomass <- c("Corn litter", 
                                 "Fungi",
                                 "Fungi sp.",
                                 "Ground",
                                 "Miscellaneous litter",
                                 "Mosses",
                                 "Mosses & lichens",
                                 "Mosses & lichens 2",
                                 "moses & lichens",
                                 "Mosses and lichens",
                                 "Lichen",
                                 "Lichens",
                                 "Other",
                                 "Other animal diggings",
                                 "Other litter",
                                 "Pine litter",
                                 "Pine cones",
                                 "pine needles",
                                 "Pine needles",
                                 "Pine twigs",
                                 "woody debris",
                                 "Woody debris",
                                 "Leaves")

maybe_plant_things_in_biomass <- c("Grass seedlings",
                                   "1st year woody",
                                   "Bryophyte",
                                   "C3 grasses",
                                   "C4 grasses",
                                   "Forb",
                                   "Forb seedlings",
                                   "Forb sp.",
                                   "Forbes",
                                   "Legumes",
                                   "Miscellaneous forb",
                                   "Miscellaneous forbs",
                                   "Miscellaneous forb 1",
                                   "Miscellaneous forb 2",
                                   "Miscellaneous grass",
                                   "Miscellaneous grasses",
                                   "Miscellaneous grasses 2",
                                   "Miscellaneous herb",
                                   "Miscellaneous herbs",
                                   "Miscellaneous herbs 2",
                                   "Miscellaneous legumes",
                                   "Miscellaneous liter",
                                   "Miscellaneous litter",
                                   "Miscellaneous rushes",
                                   "Miscellaneous sedges",
                                   "Miscellaneous sp.",
                                   "Miscellaneous sp. 2",
                                   "Miscellaneous woody tree",
                                   "Miscellaneous  woody",
                                   "Miscellaneous Woody",
                                   "Miscellaneous woody 1",
                                   "Miscellaneous woody 2",
                                   "Miscellaneous woody plants",
                                   "Miscellaneous woody plants 1",
                                   "Miscellaneous woody plants 2",
                                   "Miscellaneous woody litter",
                                   "Unknown",
                                   "Unknown cupressaceae sp.",
                                   "Unknown fabaceae",
                                   "Unknown lamiaceae",
                                   "Unknown sp.",
                                   "Sedges",
                                   "Woody")



species_list <- cbind.data.frame(
  species = species_abundance$species %>%
    factor() %>%
    levels(),
  species_clean = species_abundance$species %>%
    factor() %>%
    levels() %>%
    str_trim() %>%
    str_squish() %>%
    str_to_sentence() %>%
    gsub(pattern = " \\(\\*\\)", replacement = "") %>%
    gsub(pattern = " \\(l\\.\\)", replacement = "_") %>%
    gsub(pattern = " l\\.", replacement = "_") %>%
    gsub(pattern = "Unk ", replacement = "Unknown ") %>%
    gsub(pattern = "Unk_", replacement = "Unknown ") %>%
    gsub(pattern = " ", replacement = "_") 
) %>%
  mutate(
    category = case_when(startsWith(as.character(species), "Unknown")       ~ "unidentified_plant_things",
                         startsWith(as.character(species), "Miscellaneous") ~ "unidentified_plant_things",
                         species_clean %in% genus_sp_in_biomass           | species %in% genus_sp_in_biomass ~"identified_to_genus",
                         species_clean %in% non_plant_things_in_biomass   | species %in% non_plant_things_in_biomass ~ "non_plant_things",
                         species_clean %in% maybe_plant_things_in_biomass | species %in% maybe_plant_things_in_biomass  ~ "unidentified_plant_things")
  ) %>%
  arrange(species_clean)

species_list %>% select(species_clean, category) %>% unique()

species_abundance %>% filter(species %in% (species_list %>% filter(str_detect(species, " \\(\\*\\)")))$species) %>% select(site, higher_order_organization) %>% unique()

# write.csv(species_abundance, file.path(L2_dir, "./species_abundance.csv"), row.names=F)

