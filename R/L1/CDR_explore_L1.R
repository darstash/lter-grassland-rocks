# TITLE:          LTER Grassland Rock: KBS ANPP biomass and plant composition cleanup
# AUTHORS:        Seraina Cappelli
# COLLABORATORS:  LTER synthesis group
# DATA INPUT:     Data imported as txt files from shared Google drive L0 folder*
# DATA OUTPUT:    
# PROJECT:        LTER Grassland Rock
# DATE:           11/28/2023

# R VERSION:      R-4.3.2, Rstudio: 2023.06.2 Build 561
# MUSIC:          Singuläre Frau - Jeffi Lou
#                 Não Ao Marco Temporal - Esperanza Spalding
#                 run>> - IKAN HYU

# * Species level Biomass (and cover data if available) from those CDR experiments:
# e001: https://cedarcreek.umn.edu/research/experiments/e001
#       Variables
#       exp:         experiment e001 -> 1
#       year:        1982 - 2018 (2012 missing)
#       field:       A: fencing given up after 2004, burned annually after 2005
#                    B: fencing given up after 2004, burned annually after 2005
#                    C: fencing given up after 2004, burned annually after 2005 (half of the plots split into e172, that continued to be fenced)
#                    D: burned after 1987 every 2-3 years, fencing given up after 2004
#       plot:        1 - 54 in each field
#       n_trt:       1 - 9 (1: control)
#       n_add:       amount of n added
#       nitr_add:    
#       n_atm_n_add: amount of n added plus atmospheric deposition (n_add + 1)
#       species:     unique species name
#       biomass:     
# e054: https://cedarcreek.umn.edu/research/experiments/e054
#       Variables
#       exp:         experiment e054 -> 54
#       year:        1988 - 2018
#       oldfield:    
#       plot:        
#       transect:   
#       yearab:      Year when the field was abandonned from agriculture (between 1927 and 2015)  
#       species:     unique species name
#       biomass:     
# e097: https://cedarcreek.umn.edu/research/experiments/e097
#       Variables
#       field:       
#       exp:        
#       plot:        
#       date:        
#       ntrt:        
#       nadd:       
#       ntrtreceived:
#       species:    
#       biomass:     
# e245: https://cedarcreek.umn.edu/research/experiments/e245
#       Variables
#       year:        2008 - 2021
#       date:        exact sampling date
#       plot:        1-48 (1-6 in block 1, 7-12 in block 2, ... 43-48 in block 3)
#       subplot:     after 2019 divided in half, west subplots received fertilizer treatment (see ferttrt)
#       treatment:   Control, FoliarFungicide, Insecticide, SoilDrenchFungicide, AllPesticides, Fenced
#       ferttrt:     fertility treatment y (yes), n (no)
#       species:     unique species name
#       biomass:     g.m-1 biomass collected in a 10cm x 1m strip (exact location varies between years)
#       strip:       in 2017 and 2018 two instead of just one 10cm x 1m strip was sampled for biomass

# Clear all existing data
rm(list=ls())

#Load packages
library(tidyverse) # tidyverse 2.0.0
# dplyr     1.1.3     # readr     2.1.4
# forcats   1.0.0     # stringr   1.5.1
# ggplot2   3.4.4     # tibble    3.2.1
# lubridate 1.9.3     # tidyr     1.3.0
# purrr     1.0.2     

# Set working directory 
L0_dir <- Sys.getenv("L0DIR")
list.files(L0_dir)


# Load data ####
e001_anpp <- 
  read.table(paste(L0_dir, "e001_Plant_aboveground_biomass_data.txt", sep = "/"), 
             sep  = "\t", 
             skip = 1
  ) %>%
  rename("exp"         = "V1", 
         "year"        = "V2", 
         "field"       = "V3",
         "plot"        = "V4",
         "n_trt"       = "V5", 
         "n_add"       = "V6", 
         "nitr_add"    = "V7", 
         "n_atm_n_add" = "V8",  
         "species"     = "V9", 
         "biomass"     = "V10"
  )

e054_anpp <- 
  read.table(paste(L0_dir, "e054_Plant_aboveground_biomass_data.txt", sep = "/"), 
             sep  = "\t",
             skip = 1
  ) %>%
  rename("exp"     = "V1", 
         "year"    = "V2", 
         "oldfield"= "V3",
         "plot"    = "V4",
         "transect"= "V5", 
         "yearab"  = "V6", 
         "species" = "V7", 
         "biomass" = "V8"
         )

e097_anpp <- 
  read.table(paste(L0_dir, "e097_Plant_aboveground_biomass_data.txt", sep = "/"), 
             sep  = "\t", 
             skip = 1
  ) %>%
  rename("field"       = "V1", 
         "exp"         = "V2", 
         "plot"        = "V3",
         "date"        = "V4",
         "ntrt"        = "V5", 
         "nadd"        = "V6", 
         "ntrtreceived"= "V7", 
         "species"     = "V8",  
         "biomass"     = "V9"
  )


e245_anpp <- 
  read.table(paste(L0_dir, "e245_Plant_aboveground_biomass_data.txt", sep = "/"), 
             sep  = "\t", 
             skip = 1
  ) %>%
  rename("year"     = "V1", 
         "date"     = "V2", 
         "plot"     = "V3",
         "subplot"  = "V4",
         "treatment"= "V5", 
         "ferttrt"  = "V6", 
         "species"  = "V7", 
         "biomass"  = "V8",  
         "strip"    = "V9"
  )


# Clean data ####
#e001 first#
str(e001_anpp)

e001_anpp <- e001_anpp %>%
  mutate(species = as.factor(species)) #to look at all the species and identify things to remove


levels(e001_anpp$species)

#remove none plant species data, is there a better way to do this?
e001_anpp <- e001_anpp %>%
  filter(species != "Corn litter") %>%
  filter(species != "Fungi") %>%
  filter(species != "Grass seedlings") %>%
  filter(species != "Leaves") %>%
  filter(species != "Lichens") %>%
  filter(species != "Miscellaneous liter") %>%
  filter(species != "Miscellaneous forb") %>%
  filter(species != "Miscellaneous litter") %>%
  filter(species != "Miscellaneous woody litter") %>%
  filter(species != "moses & lichens") %>%
  filter(species != "Mosses") %>%
  filter(species != "Mosses & lichens") %>%
  filter(species != "Mosses & lichens 2") %>%
  filter(species != "Pine litter") %>%
  filter(species != "pine needles") %>%
  filter(species != "Pine needles") %>%
  filter(species != "Pine twigs") %>%
  filter(species != "woody debris") %>%
  filter(species != "Woody debris") %>%
  filter(species != "Miscellaneous  woody") %>%
  filter(species != "Miscellaneous forb 1") %>%
  filter(species != "Miscellaneous forb 2") %>%
  filter(species != "Miscellaneous grass") %>%
  filter(species != "Miscellaneous grasses") %>%
  filter(species != "Miscellaneous herbs") %>%
  filter(species != "Miscellaneous herbs 2") %>%
  filter(species != "Miscellaneous sedges") %>%
  filter(species != "Miscellaneous sp.") %>%
  filter(species != "Miscellaneous woody 1") %>%
  filter(species != "Miscellaneous woody 2") %>%
  filter(species != "Miscellaneous woody plants") %>%
  filter(species != "Miscellaneous woody plants 1") %>%
  filter(species != "Miscellaneous woody plants 2") %>%
  filter(species != "Miscellaneous woody tree") %>%
  filter(species != "Woody debris")

#combine rows that have same species but different biomass - this would be due to error I assume (they measured biomass of a species and entered it, then had another of the same species and added that entry as well)
names(e001_anpp)
e001_anpp <- e001_anpp %>%
  group_by(exp, year, field, n_trt, n_add, plot, species) %>% # this removes nitr_add and n_atmn_n_add columns which we don't want for cleaned data
  summarize(biomass=sum(biomass))

# Summarize data ----
# Calculate ANPP and species richness for each plot in a given year
e001_anpp <- e001_anpp %>%
  group_by(exp, year, field, n_trt, n_add, plot) %>%
  summarize(anpp = sum(biomass), #anpp
            richness = n())      #species richness

View(e001_anpp) #note some plots have species richness of 1 which is not an error based on looking at raw data.

#next add other data (precip and temp?) needed for cleaned data. Also work on other datasets.
