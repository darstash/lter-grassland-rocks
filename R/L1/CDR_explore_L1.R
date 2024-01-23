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
#       year:        1982 - 2022 (2012 missing)
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
# e002: https://cedarcreek.umn.edu/research/experiments/e002
#       Variables
#       exp:         experiment e002 -> 2
#       year:        1982 - 2022 (2012 missing)
#       field:       A: fencing given up after 2004, burned annually after 2005
#                    B: fencing given up after 2004, burned annually after 2005
#                    C: fencing given up after 2004, burned annually after 2005 
#                    D: burned after 1987 every 2-3 years, fencing given up after 2004
#       plot:        1 - 54 in each field
#       n_trt:       1 - 9 (1: control)
#       n_add:       amount of n added
#       nitr_add:    
#       n_atm_n_add: amount of n added plus atmospheric deposition (n_add + 1)
#       species:     unique species name
#       biomass:
#       note: e001 and e002 data in same initial df. only difference experimental design wise is that e002 was disced (plowed) before establishing nutrient treatments
# e054: https://cedarcreek.umn.edu/research/experiments/e054
#       Variables
#       exp:         experiment e054 -> 54
#       year:        1988 - 2022
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

library(janitor)

# Set working directory 
L0_dir <- Sys.getenv("L0DIR")
list.files(L0_dir)

# Load data ####
species_list_CDR <- read.csv(paste(L0_dir, "CDR/cc_plant_species.csv", sep = "/"))

e001e002_anpp <- read.csv(paste(L0_dir, "CDR/E001 E002 Aboveground Biomass for ML through 2022.csv", sep = "/"))

#e001_anpp <-
#  read.table(paste(L0_dir, "CDR/e001_Plant_aboveground_biomass_data.txt", sep = "/"), 
#             sep  = "\t", 
#             skip = 1
#  ) %>%
#  rename("exp"         = "V1", 
#         "year"        = "V2", 
#         "field"       = "V3",
#         "plot"        = "V4",
#         "n_trt"       = "V5", 
#         "n_add"       = "V6", 
#         "nitr_add"    = "V7", 
#         "n_atm_n_add" = "V8",  
#         "species"     = "V9", 
#         "biomass"     = "V10"
#  )

e054_anpp <- read.csv(paste(L0_dir, "CDR/e54_biomass_1221_ML.csv", sep = "/")) 

#e054_anpp <- 
#  read.table(paste(L0_dir, "CDR/E001 E002 Aboveground Biomass for ML through 2022.csv", sep = "/"), 
#             sep  = "\t",
#             skip = 1
#  ) %>%
#  rename("exp"     = "V1", 
#         "year"    = "V2", 
#         "oldfield"= "V3",
#         "plot"    = "V4",
#         "transect"= "V5", 
#         "yearab"  = "V6", 
#         "species" = "V7", 
#         "biomass" = "V8"
#         )

e097_anpp <- 
  read.table(paste(L0_dir, "CDR/e097_Plant_aboveground_biomass_data.txt", sep = "/"), 
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


e245_anpp <- read.delim(paste(L0_dir, "CDR/e245_Plant_aboveground_biomass_data.txt", sep = "/")) %>%
  mutate(Year = as.integer(paste(Year)),
         Plot = as.integer(paste(Plot)),
         Treatment = factor(Treatment, 
                            levels = c("Control", "Insecticide", 
                                       "FoliarFungicide", 
                                       "SoilDrenchFungicide", 
                                       "AllPesticides", "Fenced")), 
         FertTrt = factor(FertTrt),
         Species = factor(Species),
         Species = str_squish(Species),
         Species = str_to_sentence(Species),
         Species = as.factor(Species)) %>%
  # the data of the first year of experiment (2008) was not sorted to species.
  filter(!Year %in% c(NA,1,2,3,4,5,6, 2008)) %>%
  # in 2017 and 2018 biomass was collected in two strips and values need to be 
  # averaged to get one measurement per plant and year
  group_by(Year, Plot, Subplot, Treatment, FertTrt, Species) %>%
  reframe(Mass.g.m.2. = ifelse(Year %in% c(2017, 2018), sum(Mass.g.m.2.)/2, sum(Mass.g.m.2.))) %>%
  # at some point, the plots were split in two subplots and the western subplot started to get fertilized. This has not been running for long enough to be considered here.
  filter(!Subplot %in% "West") %>%
  clean_names() 


# Clean data ####
##e001 and e002####
str(e001e002_anpp)
names(e001e002_anpp)

#adding study information and renaming columns to match master data sheet
e001e002_anpp <- e001e002_anpp %>%
  mutate(site="CDR",
         higher_order_organization = paste("Experiment", Exp, " field", Field), # note I added the field string, so that it is not a random A. Not sure this is needed
         species = as.factor(Species)) #to look at all the species and identify things to remove

#make new column that designates if fertilized or not
e001e002_anpp <- e001e002_anpp %>%
  mutate(nutrients_added = ifelse(NTrt %in% 9, "no_fertilizer", 
                                  ifelse(NTrt %in% 1, "PK+", "NPK+")),
         nitrogen_amount = NAdd)

#make new column designating fence treatment and burn treatment (based on cdr experimental design info from website - see flowchart for e001 and e002)
e001e002_anpp <- e001e002_anpp %>%
  mutate(grazing = case_when(Exp == 1 & Year < 2004 ~ "ungrazed",
                             Exp == 1 & Year >= 2004 ~ "grazed",
                             Exp == 2 & Year < 2004 ~ "ungrazed",
                             Exp == 2 & Year >= 2004 ~ "grazed"),
         fire_frequency = case_when(Exp == 1 & Field == "A" & Year < 2005 ~ 0,
                                    Exp == 1 & Field == "A" & Year >= 2005 ~ 1,
                                    Exp == 1 & Field == "B" & Year < 2005 ~ 0,
                                    Exp == 1 & Field == "B" & Year >= 2005 ~ 1,
                                    Exp == 1 & Field == "C" & Year < 2005 ~ 0,
                                    Exp == 1 & Field == "C" & Year >= 2005 ~ 1,
                                    Exp == 1 & Field == "D" & Year < 1987 ~ 0,
                                    Exp == 1 & Field == "D" & Year >= 1987 ~ 2.5,
                                    Exp == 2 & Field == "A" ~ 0,
                                    Exp == 2 & Field == "B" ~ 0,
                                    Exp == 2 & Field == "C" ~ 0,))

#remove none plant species data, is there a better way to do this? --> use 
# species list first and then create a second lists of things to kick out (list
# started below). These lists can then be used for all data sets.


# find all "species" names that are not an identified species -> fix those
merge(e001e002_anpp,
      species_list_CDR %>% 
        select(Species, ITISRecognizedName) %>%
        rename(species = Species),
      by = "species",
      all.x = T) %>% 
  select(species, ITISRecognizedName) %>%
  unique() %>%
  arrange(species) %>%
  filter(ITISRecognizedName %in% c(NA, ""))

merge(e245_anpp,
      species_list_CDR %>% 
        select(Species, ITISRecognizedName) %>%
        rename(species = Species),
      by = "species",
      all.x = T) %>% 
  select(species, ITISRecognizedName) %>%
  unique() %>%
  arrange(species) %>%
  filter(ITISRecognizedName %in% c(NA, ""))

genus_sp_in_biomass <- c("Allium sp.",
                         "Apocynum sp.",
                         "Arabis sp.",
                         "Aristida sp.",
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
                         "Poa sp.",
                         "Polygala sp.",
                         "Polygonatum sp.",
                         "Potentilla sp.",
                         "Prunus sp.",
                         "Quercus borealis-ellipsoidalis",
                         "Quercus sp.",
                         "Ranunculus sp.",
                         "Rhus sp.",
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
                         "Trifolium sp.",
                         "Viola sp.")

non_plant_things_in_biomass <- c("Corn litter", 
                                 "Fungi",
                                 "Miscellaneous litter",
                                 "Mosses",
                                 "Mosses & lichens",
                                 "Mosses & lichens 2",
                                 "moses & lichens",
                                 "Lichens",
                                 "Pine litter",
                                 "Pine cones",
                                 "pine needles",
                                 "Pine needles",
                                 "Pine twigs",
                                 "woody debris",
                                 "Woody debris",
                                 "Woody",
                                 "Leaves")

maybe_plant_things_in_biomass <- c("Grass seedlings",
                                   "1st year woody",
                                   "C3 grasses",
                                   "C4 grasses",
                                   "Forb seedlings",
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
                                   "Sedges")


# correct typos & kick out things that are in the kick out things on lists
e001e002_anpp = e001e002_anpp %>%
  merge(.,
        species_list_CDR %>% 
          select(Species, ITISRecognizedName) %>%
          rename(species = Species),
        by = "species",
        all.x = T) %>%
  mutate(species_fixed = species,
         species_fixed = ifelse(!ITISRecognizedName %in% c("", NA), ITISRecognizedName, species_fixed),
         species_fixed = gsub(species_fixed, pattern = "carex sp.",             replacement = "Carex sp."),
         species_fixed = gsub(species_fixed, pattern = "Chenopodium Album",     replacement = "Chenopodium album"),
         species_fixed = gsub(species_fixed, pattern = "cyperus sp.",           replacement = "Cyperus sp."),
         species_fixed = gsub(species_fixed, pattern = "Heliantus Laetifiorus", replacement = "Helianthus laetiflorus"),
         species_fixed = gsub(species_fixed, pattern = "Poa Pratensis",         replacement = "Poa pratensis"),
         species_fixed = gsub(species_fixed, pattern = "Pycnan Vir",            replacement = "Pycnanthemum virginianum"),
         species_fixed = gsub(species_fixed, pattern = "ranoncolos rhombiodeus",replacement = "Ranunculus rhomboideus"),	
         species_fixed = gsub(species_fixed, pattern = "Taraxicum officinalis", replacement = "Taraxacum officinalis")
         ) %>%
  # filter(!species %in% maybe_plant_things_in_biomass) %>%
  filter(!species %in% non_plant_things_in_biomass)

#metadata df#
e001e002_metadata <- e001e002_anpp %>%
  clean_names(.) %>%
  select(site, year, plot, higher_order_organization, nutrients_added, nitrogen_amount, grazing, fire_frequency)
  
#e001 still needs to add temp, precip, and other variables that the master datasheet will have

#combine rows that have same species but different biomass - this would be due to error I assume (they measured biomass of a species and entered it, then had another of the same species and added that entry as well)
e001e002_anpp <- e001e002_anpp %>%
  group_by(Year, site, Plot, higher_order_organization, species) %>% # this removes nitr_add and n_atmn_n_add columns which we don't want for cleaned data
  summarize(Biomass.g.m2.=sum(Biomass.g.m2.))

e001e002_anpp <- e001e002_anpp %>%
  group_by(Year, site, Plot, higher_order_organization) %>%
  mutate(relative_abundance = Biomass.g.m2./sum(Biomass.g.m2.),
         abundance = Biomass.g.m2.,
         original_measurement_unit = "biomass_g/m2") %>%
  clean_names(.)

e001e002_anpp <- e001e002_anpp %>%
  select(year, site, plot, higher_order_organization, species, abundance, relative_abundance, original_measurement_unit)

#e001 and e002 anpp data is now in correct format
#metadata still needs time since last burn and weather data


##e054####
str(e054_anpp)
names(e054_anpp)

#adding study information and renaming columns to match master data sheet
e054_anpp <- e054_anpp %>%
  mutate(site="CDR",
         higher_order_organization = paste("Experiment", Exp, " field", OldField), # note I added the field string, so that it is not a random A. Not sure this is needed
         species = as.factor(Species_names), #to look at all the species and identify things to remove
         Plot = case_when(Transect == "G" ~ 1, #sampling is done once in each transect within a field
                          Transect == "R" ~ 2, #so we are denoting the transect as 1-4 where biomass was collected in this plot column
                          Transect == "W" ~ 3,
                          Transect == "Y" ~ 4)) 

#make new column that designates if fertilized or not
e054_anpp <- e054_anpp %>%
  mutate(nutrients_added = "no_fertilizer",
         nitrogen_amount = 0)

#make new column designating fence treatment and burn treatment (based on cdr experimental design info from website - see flowchart for e001 and e002)
e054_anpp <- e054_anpp %>%
  mutate(grazing = "ungrazed",
         fire_frequency = case_when(Burnt_started2007 == "Unburned" ~ 0,
                                    Burnt_started2007 == "Burned" ~ 5,))




# correct typos & kick out things that are in the kick out things on lists
e054_anpp = e054_anpp %>%
  merge(.,
        species_list_CDR %>% 
          select(Species, ITISRecognizedName) %>%
          rename(species = Species),
        by = "species",
        all.x = T) %>%
  filter(!species %in% non_plant_things_in_biomass)

#metadata df#
e054_metadata <- e054_anpp %>%
  clean_names(.) %>%
  select(site, year, plot, higher_order_organization, nutrients_added, nitrogen_amount, grazing, fire_frequency)

#e054 still needs to add temp, precip, and other variables that the master datasheet will have

#combine rows that have same species but different biomass - this would be due to error I assume (they measured biomass of a species and entered it, then had another of the same species and added that entry as well)
e054_anpp <- e054_anpp %>%
  group_by(Year, site, Plot, higher_order_organization, species) %>% # this removes nitr_add and n_atmn_n_add columns which we don't want for cleaned data
  summarize(Biomass=sum(Biomass))

e054_anpp <- e054_anpp %>%
  group_by(Year, site, Plot, higher_order_organization) %>%
  mutate(relative_abundance = Biomass/sum(Biomass),
         abundance = Biomass,
         original_measurement_unit = "biomass_g/m2") %>%
  filter(abundance != 0) %>%
  clean_names(.)

e054_anpp <- e054_anpp %>%
  select(year, site, plot, higher_order_organization, species, abundance, relative_abundance, original_measurement_unit)


#e054_anpp matches data template - metadata still needs more info - precip, temp, time since last burn, etc..

##e245###
# fix species names
e245_anpp <- e245_anpp %>%
  mutate(species = recode_factor(species,
                                 "Achillea millefolium(lanulosa)" = "Achillea millefolium (lanulosa)",
                                 "Miscellaneous forb"             = "Miscellaneous forbs",
                                 "Miscellaneous grass"            = "Miscellaneous grasses",
                                 "Taraxicum officinalis"          = "Taraxacum officinale",
                                 "Viola petatifida"               = "Viola pedatifida")) %>%
  filter(!species %in% non_plant_things_in_biomass)

# make columns match with our goal dataset
e245_anpp <- e245_anpp %>%
  group_by(year, plot, treatment) %>%
  mutate(site                      = "CDR",
         higher_order_organization = case_when(plot %in% c(1:6)   ~ "e245_block1",
                                               plot %in% c(7:12)  ~ "e245_block2", 
                                               plot %in% c(13:18) ~ "e245_block3", 
                                               plot %in% c(19:24) ~ "e245_block4",
                                               plot %in% c(25:30) ~ "e245_block5",
                                               plot %in% c(31:36) ~ "e245_block6",
                                               plot %in% c(37:42) ~ "e245_block7",
                                               plot %in% c(43:48) ~ "e245_block8"),
         original_measurement_unit = "g/m2",
         relative_abundance        = mass_g_m_2/sum(mass_g_m_2)) %>%
  rename("abundance" = mass_g_m_2) %>%
  select(year, site, plot, higher_order_organization, species, abundance, relative_abundance, original_measurement_unit, treatment)

e245_metadata <- e245_anpp %>%
  select(year, site, plot, higher_order_organization, treatment) %>%
  unique()

e245_anpp <- e245_anpp %>%
  mutate(species = recode_factor(species,
                                 "Achillea millefolium(lanulosa)" = "Achillea millefolium (lanulosa)",
                                 "Miscellaneous forb"             = "Miscellaneous forbs",
                                 "Miscellaneous grass"            = "Miscellaneous grasses",
                                 "Taraxicum officinalis"          = "Taraxacum officinale",
                                 "Viola petatifida"               = "Viola pedatifida")) %>%
  filter(!species %in% non_plant_things_in_biomass)

e245_anpp <- e245_anpp %>%
  select(!treatment)
