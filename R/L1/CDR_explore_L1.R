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
#       Biomass Data Variables
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
#       Biomass Data Variables
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
#       Biomass Data Variables
#       exp:         experiment e054 -> 54
#       year:        1988 - 2022
#       oldfield:    
#       plot:        
#       transect:   
#       yearab:      Year when the field was abandonned from agriculture (between 1927 and 2015)  
#       species:     unique species name
#       biomass:     
# e097: https://cedarcreek.umn.edu/research/experiments/e097
#       Biomass Data Variables
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
#       Biomass Data Variables
#        year:        2008 - 2021
#        date:        exact sampling date
#        plot:        1-48 (1-6 in block 1, 7-12 in block 2, ... 43-48 in block 3)
#        subplot:     after 2019 divided in half, west subplots received fertilizer treatment (see ferttrt)
#        treatment:   Control, FoliarFungicide, Insecticide, SoilDrenchFungicide, AllPesticides, Fenced
#        ferttrt:     fertility treatment y (yes), n (no)
#        species:     unique species name
#        biomass:     g.m-1 biomass collected in a 10cm x 1m strip (exact location varies between years)
#        strip:       in 2017 and 2018 two instead of just one 10cm x 1m strip was sampled for biomass
# e061: https://cedarcreek.umn.edu/research/experiments/e061
#       Biomass Data Variables
#        Field.number.letter:  Letter B. Cedar Creek Ecosystem Science Reserve is subsetted into different fields. e061 is located in field B in the  macroplots of the experiments e004 and e060. Field burned in 1989, 1990, 1992, 1995
#        Experiment.number:    61 (experiment 61 = e061)
#        Macroplot.Number:     macroplots 2, 3, 4, 6, 7, and 10 of experiment e004 and e060 were used here
#        Plot.number:          in each plot there were 3 plots (1,2,3)
#        Bird.treatment:       In each plot there was a subplot with either birds exlcuded (BEX) or birds allowed (BAL)
#        Nitrogen.Treatment:   Nitrogen treatment I: control, G: 26g/m2 ammonium nitrate added (~8.9g/m^2 nitrogen) UNSURE IF ONLY IN 1990 OR ALL YEARS
#        Sampling.year:        Sampling year. Incomplete column. Probably a mistake column (There is a Year column, too)
#        Species.Name:         Species names
#        Species.Biomass:      Species biomass in g/m2
#        Sampling.date:        Sampling date in YYMMDD format
#       Year:                 Sampling year 1989, 1990, 1996-2003
# e247: https://cedarcreek.umn.edu/research/experiments/e247
#       Biomass Data Variables
#        year:                 sampling year
#        site_name:            "Cedar Creek LTER" site name within entire NutNet
#        site_code:            "cdr.us" site code identifier within entire NutNet
#        year_trt:             year of treatment = experiment duration at time of measurement. (Not all Nutrient Network sites started their treatments at the same year)
#        trt:                  Control, NP, NPK, K, N, PK, NPK+Fence etc. categorical treatment varialbe
#        block:                block within site (5 blocks in Cedar Creek LTER NutNet site)
#        plot:                 unique plot number within Cdar Creek LTER NutNet site (1-60)
#        subplot:              a random subplot (out of 4) is used for biomass sampling (always the same per plot)
#        N:                    0: no N added, 1: N added
#        P:                    0: no P added, 1: P added
#        K:                    0: no K added, 1: K added
#        Exclose:              0: no Fence, 1: Fence
#        Nitrogen.fertilizer:  Additional treatments within the N treatments: variable amount of N 0, 1, 5, 10 gN/m2/y
#        Aboveground.biomass:  dried mass (g/m2), Biomass collected in 2 0.1 x 1m strips and sorted to functional group.
#        category:             functional group, to which biomass was sorted.
#       Cover Data Variables
#        year:                 sampling year
#        site_name:            "Cedar Creek LTER" site name within entire NutNet
#        site_code:            "cdr.us" site code identifier within entire NutNet
#        year_trt:             year of treatment = experiment duration at time of measurement. (Not all Nutrient Network sites started their treatments at the same year)
#        trt:                  Control, NP, NPK, K, N, PK, NPK+Fence etc. categorical treatment varialbe
#        block:                block within site (5 blocks in Cedar Creek LTER NutNet site)
#        plot:                 unique plot number within Cdar Creek LTER NutNet site (1-60)
#        subplot:              a random subplot (out of 4) is used for cover measurements (always the same per plot)
#        N:                    0: no N added, 1: N added
#        P:                    0: no P added, 1: P added
#        K:                    0: no K added, 1: K added
#        Exclose:              0: no Fence, 1: Fence
#        Nitrogen.fertilizer:  Additional treatments within the N treatments: variable amount of N 0, 1, 5, 10 gN/m2/y            
#        Family:               Taxonomy
#        Taxon:                Species
#        live                  dead or alive (0/1)
#        local_provenance      locally assigned provenance of taxon
#        local_lifeform        locally assigned lifeform of taxon
#        local_lifespan        locally assinged lifespan of taxon
#        functional_group      = category in biomass data
#        N_fixer               nitrogen fixin yes or no (0/1)     
#        max_cover             maximum observed percent cover of taxon in plot for year indicated

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
L1_dir <- Sys.getenv("L1DIR")
list.files(L0_dir)

# Load data ####
species_list_CDR <- read.csv(paste(L0_dir, "CDR/cc_plant_species.csv", sep = "/"))

burns <- read.csv(paste(L0_dir, "CDR/CCESR_Burn_record_by_Burn_Unit_and_Year.csv", sep = "/"))

climate <- read.csv(paste(L0_dir, "CDR/e080_Daily climate summary.csv", sep="/"))

e001e002_anpp <- read.csv(paste(L0_dir, "CDR/E001 E002 Aboveground Biomass for ML through 2022.csv", sep = "/")) %>%
  mutate(Field = ifelse(Field %in% "a", "A", Field))

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

#no longer including this as it is a subset of e001/e002
#e097_anpp <- 
#  read.table(paste(L0_dir, "CDR/e097_Plant_aboveground_biomass_data.txt", sep = "/"), 
#             sep  = "\t", 
#             skip = 1
#  ) %>%
#  rename("field"       = "V1", 
#         "exp"         = "V2", 
#         "plot"        = "V3",
#         "date"        = "V4",
#         "ntrt"        = "V5", 
#         "nadd"        = "V6", 
#         "ntrtreceived"= "V7", 
#         "species"     = "V8",  
#         "biomass"     = "V9"
#  )


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


e061_anpp <- read.delim(paste(L0_dir, "CDR/e061_Plant_aboveground_biomass_data.txt", sep = "/")) %>%
  rename("field" = "Field.number.letter",
         "experiment" =  "Experiment.number",
         "macroplot" = "Macroplot.Number",
         "plot" = "Plot.number",
         "treatment" = "Bird.treatment..1.BAL.birds.allow..2.BEX.birds.exclude.",
         "fert_trt" = "Nitrogen.Treatment",
         "species" = "Species.Name",
         "mass_g_m_2" = "Species.Biomass..g.m2." ,
         "year" = "Year"
  ) %>%
  mutate(year = as.integer(paste(year)),
         field = as.factor(field),
         macroplot = as.integer(paste(macroplot)), # Note: Warnign NAs introduced by coercion is ok. It replaces empty cells with NA
         plot = as.integer(paste(plot)),
         treatment = as.factor(toupper(treatment)),
         fert_trt = factor(fert_trt),
         species = factor(species),
         species = str_squish(species),
         species = str_to_sentence(species),
         species = as.factor(species),
         field   = "B",
         experiment = "e061") %>%
  # biomass was collected in 1989, 1990 and then consecutively only from 1996-2003. Focus on this part!
  filter(year %in% c(1996:2003))


e247_anpp <- read.delim(paste(L0_dir, "CDR/e247_Aboveground_Standing_Crop_Biomass.txt", sep = "/"))

e247_cover <- read.delim(paste(L0_dir, "CDR/e247_Plant_Species_Composition_percent_cover.txt", sep = "/")) %>%
  rename(
    "species" = "Taxon",
    "abundance" = "max_cover",
    "nitrogen_amount" = "Nitrogen.added.to.plot") %>%
  mutate(
    species = str_squish(species),
    species = str_to_sentence(species),
  )




# Clean data ####
#climate data to be included in metadata for all datasets
climate_growingseason <- climate %>%
  separate(Date, into=c("Month", "Day", "Year"), remove = F) %>%
  mutate(Date = as.POSIXct(Date, format = "%m/%d/%Y"),
         Day = yday(Date)) %>% 
  filter(Day > 114.5 & Day < 228.5) %>%
  group_by(Year) %>%
  mutate(AvgTemp.defF. = (MaxTemp.degF.+MinTemp.degF.)/2) %>% #average the daily max and min to a single daily avg value
  summarize(growtemp = mean(AvgTemp.defF.),
            growprecip = sum(Precip.inches.)) %>%
  mutate(growtemp = (growtemp-32) * (5/9), #convert temp to celcius and precipitation to mm
         growprecip = growprecip * 25.4,
         Year = as.numeric(Year)) %>%
  select(Year, growtemp, growprecip)


climate <- climate %>%
  separate(Date, into=c("Month", "Day", "Year")) %>%
  group_by(Year) %>%
  mutate(AvgTemp.defF. = (MaxTemp.degF.+MinTemp.degF.)/2) %>% #average the daily max and min to a single daily avg value
  summarize(temperature = mean(AvgTemp.defF.),
            precipitation = sum(Precip.inches.)) %>%
  mutate(temperature = (temperature-32) * (5/9), #convert temp to celcius and precipitation to mm
         precipitation = precipitation * 25.4,
         Year = as.numeric(Year)) %>%
  select(Year, temperature, precipitation) %>%
  merge(., climate_growingseason, by = "Year")

#note that 1962 is clearly an error but we don't use that year of data anyway in our datasets right???
#now can add temp (mean temp for a given year) and precip (total precip) to metadata by merging based on years

##e001 and e002####
str(e001e002_anpp)
names(e001e002_anpp)

#adding study information and renaming columns to match master data sheet
e001e002_anpp <- e001e002_anpp %>%
  mutate(site="CDR",
         higher_order_organization = paste("Experiment", Exp, " field", Field), # note I added the field string, so that it is not a random A. Not sure this is needed
         species = as.factor(Species),
         uniqueid = paste(higher_order_organization, " plot", Plot)) 

#adding climate information
e001e002_anpp <- inner_join(e001e002_anpp, climate, by="Year", multiple="all")

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

#new column for disturbance (e001 is undisturbed, e002 is disturbed)
e001e002_anpp <- e001e002_anpp %>%
  mutate(disturbance = case_when(Exp == 1 ~ "undisturbed",
                                 Exp == 2 ~ "disturbed"))
#time since last burn
df <- burns %>%
  select(Year, Field.A..e001.) %>%
  filter(!Year %in% "2000 Fall") %>%
  mutate(Year = ifelse(Year %in% "2000 Spring", "2000", Year),
         Year = as.numeric(paste(Year)),
         your_column = Field.A..e001.,
         time_since_fire = 0)

for (i in 2:nrow(df)) {
  if (is.na(df$your_column[i]) || df$your_column[i] == "") {
    df$time_since_fire[i] <- df$time_since_fire[i - 1] + 1
  } else {
    df$time_since_fire[i] <- 0
  }
}


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

# merge(e245_anpp,
#       species_list_CDR %>% 
#         select(Species, ITISRecognizedName) %>%
#         rename(species = Species),
#       by = "species",
#       all.x = T) %>% 
#   select(species, ITISRecognizedName) %>%
#   unique() %>%
#   arrange(species) %>%
#   filter(ITISRecognizedName %in% c(NA, ""))
# 
# merge(e061_anpp,
#       species_list_CDR %>% 
#         select(Species, ITISRecognizedName) %>%
#         rename(species = Species),
#       by = "species",
#       all.x = T) %>% 
#   select(species, ITISRecognizedName) %>%
#   unique() %>%
#   arrange(species) %>%
#   filter(ITISRecognizedName %in% c(NA, ""))

#merge(e247_cover,
#      species_list_CDR %>% 
#        select(Species, ITISRecognizedName) %>%
#        rename(species = Species),
#      by = "species",
#      all.x = T) %>% 
#  select(species, ITISRecognizedName) %>%
##  unique() %>%
#  arrange(species) %>%
#  filter(ITISRecognizedName %in% c(NA, "")) %>%
#  filter(!species %in% genus_sp_in_biomass & !species %in% non_plant_things_in_biomass & !species %in% maybe_plant_things_in_biomass)

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
                                 "Woody",
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
  merge(., 
        df %>% select(Year, time_since_fire), 
        by.y = "Year",
        all.x = TRUE) %>%
  clean_names(.) %>%
  mutate(diversity_manipulated = "naturally_assembled",
         treatment = ifelse(n_trt %in% 1, "control", "treatment"),
         source = "https://doi.org/10.6073/pasta/2eba7aac6b347d27a92208e03fd3f8ea; https://doi.org/10.6073/pasta/66724d71711b80d520fa33a690f962b2",
         treatment_comment = "") %>%
  select(site, year, plot, uniqueid, higher_order_organization, temperature, precipitation, treatment, growtemp, growprecip,
         nutrients_added, nitrogen_amount, grazing, fire_frequency, time_since_fire,
         disturbance, source, treatment_comment, diversity_manipulated)

#e001 still needs to add temp, precip, and other variables that the master datasheet will have

#combine rows that have same species but different biomass - this would be due to error I assume (they measured biomass of a species and entered it, then had another of the same species and added that entry as well)
e001e002_anpp <- e001e002_anpp %>%
  group_by(Year, site, Plot, uniqueid, higher_order_organization, species) %>% # this removes nitr_add and n_atmn_n_add columns which we don't want for cleaned data
  summarize(Biomass.g.m2.=sum(Biomass.g.m2.))

e001e002_anpp <- e001e002_anpp %>%
  group_by(Year, site, Plot, uniqueid, higher_order_organization) %>%
  mutate(relative_abundance = Biomass.g.m2./sum(Biomass.g.m2.),
         abundance = Biomass.g.m2.,
         original_measurement_unit = "biomass_g/m2") %>%
  clean_names(.)

e001e002_anpp <- e001e002_anpp %>%
  select(year, site, plot, higher_order_organization, uniqueid, species, abundance, relative_abundance, original_measurement_unit)

#e001 and e002 anpp data is now in correct format

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
                          Transect == "Y" ~ 4),
         uniqueid = paste(higher_order_organization," plot", Plot)) 

#adding climate information
e054_anpp <- inner_join(e054_anpp, climate, by="Year", multiple="all")

#make new column that designates if fertilized or not
e054_anpp <- e054_anpp %>%
  mutate(nutrients_added = "no_fertilizer",
         nitrogen_amount = 0)

#make new column designating fence treatment and burn treatment (based on cdr experimental design info from website - see flowchart for e001 and e002)
e054_anpp <- e054_anpp %>%
  mutate(grazing = "ungrazed",
         fire_frequency = case_when(Burnt_started2007 == "Unburned" ~ 0,
                                    Burnt_started2007 == "Burned" ~ 5,),
         fire_frequency = ifelse(Year < 2007, 0, fire_frequency))




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
  mutate(disturbance = "undisturbed",
         diversity_manipulated = "naturally_assembled",
         treatment = "control",
         treatment_comment = "fields differ in time since abandonment and after 2007 in burn/no burn treatment, but we entered all plots as control, since there is no real control vs. treatment contrast.",
         source = "https://doi.org/10.6073/pasta/02d38edbe0860ef0a0555ff3e495ca1a",
         time_since_fire = NA,
         ) %>%
  select(site, year, plot, higher_order_organization, uniqueid, 
         temperature, precipitation, growtemp, growprecip, treatment, disturbance,
         nutrients_added, nitrogen_amount, grazing, fire_frequency, time_since_fire,
         source, treatment_comment, diversity_manipulated)


df <- burns %>%
  select(Year,  X108) %>% #26 unburned
  filter(!Year %in% "2000 Fall") %>%
  mutate(Year = ifelse(Year %in% "2000 Spring", "2000", Year),
         Year = as.numeric(paste(Year)),
         your_column = X108,
         time_since_fire = 0)

for (i in 2:nrow(df)) {
  if (is.na(df$your_column[i]) || df$your_column[i] == "") {
    df$time_since_fire[i] <- df$time_since_fire[i - 1] + 1
  } else {
    df$time_since_fire[i] <- 0
  }
}

df_full <- data.frame()
colnames(df_full) <- c("Year", "time_since_fire", "OldField")

for (j in paste("X", unique(e054_anpp$OldField[e054_anpp$OldField != 26]), sep = "")){
  df <- burns %>%
    select(Year, all_of(j)) %>%
    filter(!Year %in% "2000 Fall") %>%
    mutate(Year = ifelse(Year %in% "2000 Spring", "2000", Year),
           Year = as.numeric(paste(Year)),
           time_since_fire = 0) %>%
    rename("your_column" = j)
  
  
  for (i in 2:nrow(df)) {
    if (is.na(df$your_column[i]) || df$your_column[i] == "") {
      df$time_since_fire[i] <- df$time_since_fire[i - 1] + 1
    } else {
      df$time_since_fire[i] <- 0
    }
  }

  df$OldField <- rep(gsub(j, pattern = "X", replacement = ""))
  
  df_full <- rbind.data.frame(df_full, df[,c("Year", "time_since_fire", "OldField")])
  }

#e054 still needs to time since fire and other variables that the master datasheet will have

#combine rows that have same species but different biomass - this would be due to error I assume (they measured biomass of a species and entered it, then had another of the same species and added that entry as well)
e054_anpp <- e054_anpp %>%
  group_by(Year, site, Plot, higher_order_organization, uniqueid, species) %>% # this removes nitr_add and n_atmn_n_add columns which we don't want for cleaned data
  summarize(Biomass=sum(Biomass))

e054_anpp <- e054_anpp %>%
  group_by(Year, site, Plot, higher_order_organization, uniqueid) %>%
  mutate(relative_abundance = Biomass/sum(Biomass),
         abundance = Biomass,
         original_measurement_unit = "biomass_g/m2") %>%
  filter(abundance != 0) %>%
  clean_names(.)

e054_anpp <- e054_anpp %>%
  select(year, site, plot, higher_order_organization, uniqueid, species, abundance, relative_abundance, original_measurement_unit)


#e054_anpp matches data template - metadata still needs more info - time since last burn, etc..

##e245####

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
         relative_abundance        = mass_g_m_2/sum(mass_g_m_2),
         uniqueid = paste(higher_order_organization," plot", plot)) %>%
  rename("abundance" = mass_g_m_2) %>%
  select(year, site, plot, higher_order_organization, uniqueid, species, abundance, relative_abundance, original_measurement_unit, treatment)

#adding climate information
e245_anpp <- climate %>%
  mutate(year = Year) %>%
  select(year, temperature, precipitation, growtemp, growprecip) %>%
  inner_join(e245_anpp, ., by="year", multiple="all")

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




### create metadata for this experiment

# e245 is in burn unit 108. use a subset of the burns data to calculate years since fire for each year in the burn unit 108
df <- burns %>%
  select(Year, X108) %>%
  filter(!Year %in% "2000 Fall") %>%
  mutate(Year = ifelse(Year %in% "2000 Spring", "2000", Year),
         Year = as.numeric(paste(Year)),
         your_column = X108,
         time_since_fire = 0)

for (i in 2:nrow(df)) {
  if (is.na(df$your_column[i]) || df$your_column[i] == "") {
    df$time_since_fire[i] <- df$time_since_fire[i - 1] + 1
  } else {
    df$time_since_fire[i] <- 0
  }
}


e245_metadata <- e245_anpp %>%
  select(year, site, plot, higher_order_organization, uniqueid, temperature, precipitation, growtemp, growprecip, treatment) %>%
  rename(treatment_comment = treatment) %>%
  unique() %>%
  mutate(
    treatment = ifelse(treatment_comment %in% "Control", "control", "treatment"),
    nutrients_added = "no_fertilizer",
    nitrogen_amount = 0,
    disturbance = "undisturbed",
    grazing = ifelse(treatment_comment %in% "Fenced", "ungrazed", "grazed"),
    fire_frequency = 1/3,
    diversity_manipulated = "naturally_assembled",
    source = "https://doi.org/10.6073/pasta/303607d5f92929a4b20ba127c47d21f0 (Accessed 2024-05-08"
  ) %>%
  merge(., 
        df %>% select(Year, time_since_fire), 
        by.x = "year", 
        by.y = "Year",
        all.x = TRUE) # %>%
# merge(.,
#       meteodata)
# select(year, site, plot, higher_order_organization, temperature, 
#        precipitation, nutrients_added,  treatment, nutrients_added, 
#        nitrogen_amount, disturbance, grazing, fire_frequency, time_since_fire,
#        treatment_comment, diversity_manipulated)


##e061####
e061_anpp  <- e061_anpp %>%
  mutate(site = "CDR",
         #plot = paste("e061_macroplot", macroplot, "_plot", plot, "_", treatment, sep= ""),
         higher_order_organization = paste("e061_fieldB","macroplot", macroplot,sep= ""),
         uniqueid = paste(higher_order_organization," plot", plot),
         abundance = mass_g_m_2,
         original_measurement_unit = "biomass_g/m2",
         species = ifelse(species %in% "Chenopodiium album", "Chenopodium album", paste(species)),
         species = ifelse(species %in% "Polyganum cilinode", "Polygonum cilinode", paste(species)),
         species = ifelse(species %in% "Schizhyrium scoparium", "Schizachyrium scoparium", paste(species)),
         species = ifelse(species %in% "Selaria lutescens (glauca)", "Setaria lutescens (glauca)", paste(species))) %>%
  filter(!species %in% non_plant_things_in_biomass) %>%
  group_by(site, field, experiment, macroplot, plot, treatment, year) %>%
  mutate(relative_abundance        = abundance/sum(abundance )) 

  
# create metadata  

# notes for metadata: 
# e061 is located in field B in the same macroplots than e004 and e060. Field B
# was fenced until 2004. So there is no grazing. The  burn information data 
# set has a column names "Part of field B(e004?)". This field was last burned in
# 1995. So it can be classified as fire frequency of zero (data used here starts 
# 1996). Time since fire is therefore year-1995


e061_metadata  <- e061_anpp %>%
  ungroup() %>%
  select(year, site, plot, higher_order_organization, uniqueid, treatment, fert_trt) %>%
  unique() %>%
  mutate(treatment_comment = ifelse(treatment %in% "BEX", "birds excluded", "control"),
         treatment = ifelse(treatment %in% "BEX", "treatment", "control"),
         nutrients_added = "no_fertilizer",
         nitrogen_amount = 0,
         disturbance = "undisturbed",
         grazing = "ungrazed",
         diversity_manipulated = "naturally_assembled",
         fire_frequency = 0,
         time_since_fire = year-1995,
         source = "e061") # %>%

#adding climate information
e061_metadata <- climate %>%
  mutate(year = Year) %>%
  select(year, temperature, precipitation, growtemp, growprecip) %>%
  inner_join(e061_metadata, ., by="year", multiple="all")

e061_metadata <- e061_metadata %>%
  select(-fert_trt)


e061_anpp <- e061_anpp %>%
  ungroup() %>%
  select(year, site, plot, higher_order_organization,uniqueid, species, abundance, relative_abundance, original_measurement_unit)


## e247 ####
e247_cover <- e247_cover %>%
  mutate(site = "CDR",
         #plot = paste("e247_block_", block, "plot_", plot, sep = ""),
         higher_order_organization = paste("e247_block_", block, sep = ""),
         uniqueid = paste(higher_order_organization," plot", plot),
         # these species names do not exist in the cedar creek species list. One could look them up and match them, so they would be easier to match for other things.
         # species = ifelse(species %in% "Achillea millefolium",             "", species),
         # species = ifelse(species %in% "Agrostis capillaris",              "", species),
         # species = ifelse(species %in% "Agrostis stolonifera",             "", species),
         # species = ifelse(species %in% "Ambrosia artemisiifolia",          "", species),
         # species = ifelse(species %in% "Anaphalis margaritacea",           "", species),
         # species = ifelse(species %in% "Andropogon gerardii",              "", species),
         # species = ifelse(species %in% "Anemone canadensis",               "", species),
         # species = ifelse(species %in% "Conyza canadensis",                "", species),
         # species = ifelse(species %in% "Cyperus grayi",                    "", species),
         # species = ifelse(species %in% "Cyperus lupulinus",                "", species),
         # species = ifelse(species %in% "Echinacea serotina",               "", species),
         # species = ifelse(species %in% "Elymus repens",                    "", species),
         # species = ifelse(species %in% "Fallopia convolvulus",             "", species),
         # species = ifelse(species %in% "Lychnis latifolia ssp. Alba",      "", species),
         # species = ifelse(species %in% "Panicum acuminatum",               "", species),
         # species = ifelse(species %in% "Pennisetum glaucum",               "", species),
         # species = ifelse(species %in% "Rudbeckia hirta",                  "", species),
         # species = ifelse(species %in% "Rudbeckia hirta var. Pulcherrima", "", species),
         # species = ifelse(species %in% "Silene latifolia",                 "", species),
         # species = ifelse(species %in% "Symphyotrichum boreale",           "", species),
         # species = ifelse(species %in% "Toxicodendron radicans",           "", species),
         # species = ifelse(species %in% "Tragopogon dubius",                "", species),
         original_measurement_unit = "%cover") %>%
  group_by(year, plot, uniqueid) %>%
  mutate(relative_abundance        = abundance/sum(abundance ))

e247_metadata <- e247_cover %>%
  select(year, site, plot, higher_order_organization, uniqueid, trt, Exclose, nitrogen_amount) %>%
  unique() %>%
  mutate(source = "e247",
         treatment = ifelse(trt %in% "Control", "control", "treatment"),
         nutrients_added = trt,
         nutrients_added = gsub(nutrients_added, pattern = "+Fence",  replacement = ""),
         nutrients_added = gsub(nutrients_added, pattern = "Fence",   replacement = "no_fertlizier"),
         nutrients_added = gsub(nutrients_added, pattern = "Control", replacement = "no_fertilizer"),
         disturbance = "undisturbed",
         grazing = ifelse(trt %in% c("NPK+Fence", "Fence"), "ungrazed", "grazed"),
         fire_frequency = NA,
         time_since_fire = NA,
         diversity_manipulated = "naturally_assembled",
         measurement_scale_biomass = "0.2m^2",
         measurement_scale_cover = "1m^2",
         source = "https://doi.org/10.6073/pasta/303607d5f92929a4b20ba127c47d21f0",
         treatment_comment = NA) # %>%

#adding climate information
e247_metadata <- climate %>%
  mutate(year = Year) %>%
  select(year, temperature, precipitation, growtemp, growprecip) %>%
  inner_join(e247_metadata, ., by="year", multiple="all")

e247_metadata <- e247_metadata %>%
  select(year, site, source, higher_order_organization, uniqueid, plot, treatment,
         nitrogen_amount, nutrients_added, disturbance, grazing, fire_frequency,
         time_since_fire, diversity_manipulated, temperature, precipitation, growtemp, growprecip, measurement_scale_biomass, measurement_scale_cover,
         source, treatment_comment)


e247_cover <- e247_cover %>%
  select(year, site, plot, higher_order_organization, uniqueid, species, abundance, relative_abundance, original_measurement_unit)

#ANPP and Richness#####
e001e002_metrics <- e001e002_anpp %>%
  group_by(year, site, higher_order_organization, plot,  uniqueid, original_measurement_unit) %>%
  summarize(plot_biomass = sum(abundance),
            richness = n()) %>%
  mutate(measurement_scale_biomass = "0.3m^2",
         measurement_scale_cover = NA)

e054_metrics <- e054_anpp %>%
  group_by(year, site, higher_order_organization, plot,  uniqueid, original_measurement_unit) %>%
  summarize(plot_biomass = sum(abundance),
            richness = n()) %>%
  mutate(measurement_scale_biomass = "0.3m^2",
         measurement_scale_cover = NA)

e245_metrics <- e245_anpp %>%
  group_by(year, site, higher_order_organization, plot,  uniqueid, original_measurement_unit) %>%
  summarize(plot_biomass = sum(abundance),
            richness = n()) %>%
  mutate(measurement_scale_biomass = "0.2m^2",
         measurement_scale_cover = NA)

e061_metrics <- e061_anpp %>%
  group_by(year, site, higher_order_organization, plot,  uniqueid, original_measurement_unit) %>%
  summarize(plot_biomass = sum(abundance),
            richness = n()) %>%
  mutate(measurement_scale_biomass = "0.3m^2",
         measurement_scale_cover = NA)

#for e247 need to calculate richness from cover data and biomass from anpp data then merge dataframes
e247_anpp <- e247_anpp %>%
  group_by(year, block, plot) %>%
  summarize(plot_biomass = sum(Aboveground.biomass)) %>%
  ungroup() %>%
  select(year, plot, plot_biomass) %>%
  mutate(measurement_scale_biomass = "0.2m^2")

e247_richness <- e247_cover %>%
  group_by(year, site, plot, higher_order_organization, uniqueid, original_measurement_unit) %>%
  summarize(richness = n()) %>%
  mutate(measurement_scale_cover = "1m^2")

e247_metrics <- merge(e247_richness, e247_anpp, by=c("year", "plot"))

#Combine CDR plotlevel  datasets####

cdr_data <- e001e002_metrics %>% mutate(dataset = "e001_e002") %>%
  rbind(e054_metrics %>% mutate(dataset = "e054")) %>%
  rbind(e245_metrics %>% mutate(dataset = "e245")) %>%
  rbind(e061_metrics %>% mutate(dataset = "e061")) %>%
  rbind(e247_metrics %>% mutate(dataset = "e247"))

#so this has all the bare minimum variables - no metadata information - need to confirm and work on checking off all variables for each data set - then can consider making a master metadata df too????

#Combine CDR specieslevel datasets ######
cdr_sp_data <- e001e002_anpp %>% mutate(dataset = "e001_e002") %>%
  select(year,         site,         dataset,      plot,         higher_order_organization,
         uniqueid,     species,      abundance,    relative_abundance,
         original_measurement_unit) %>%
  rbind(e054_anpp %>% mutate(dataset = "e054") %>%
          select(year,         site,         dataset,      plot,         higher_order_organization,
                 uniqueid,     species,      abundance,    relative_abundance,
                 original_measurement_unit)) %>%
  rbind(e061_anpp %>% mutate(dataset = "e061") %>%
          select(year,         site,         dataset,      plot,         higher_order_organization,
                 uniqueid,     species,      abundance,    relative_abundance,
                 original_measurement_unit)) %>%
  rbind(e245_anpp %>% ungroup() %>% mutate(dataset = "e245") %>%
          select(year,         site,         dataset,      plot,         higher_order_organization,
                 uniqueid,     species,      abundance,    relative_abundance,
                 original_measurement_unit)) %>%
  rbind(e247_cover %>% mutate(dataset = "e247") %>%
          select(year,         site,         dataset,      plot,         higher_order_organization,
                 uniqueid,     species,      abundance,    relative_abundance,
                 original_measurement_unit))

#Combine CDR specieslevel abundance (biomass based/cover based in case of nutnet)


#Combine CDR metadata #####
cdr_metadata <- e001e002_metadata %>% mutate(dataset = "e001_e002") %>%
  select(year,            site,            dataset,       plot,          higher_order_organization,
         uniqueid,        temperature,     precipitation, growtemp, growprecip, treatment,
         nutrients_added, nitrogen_amount, disturbance,   grazing,
         fire_frequency,  time_since_fire, source,        treatment_comment,
         diversity_manipulated) %>%
  rbind(e054_metadata %>% mutate(dataset = "e054") %>%
          select(year,            site,            dataset,       plot,          higher_order_organization,
                 uniqueid,        temperature,     precipitation, growtemp, growprecip, treatment,
                 nutrients_added, nitrogen_amount, disturbance,   grazing,
                 fire_frequency,  time_since_fire, source,        treatment_comment,
                 diversity_manipulated)) %>%
  rbind(e061_metadata %>% mutate(dataset = "e061") %>%
          select(year,            site,            dataset,       plot,          higher_order_organization,
                 uniqueid,        temperature,     precipitation, growtemp, growprecip, treatment,
                 nutrients_added, nitrogen_amount, disturbance,   grazing,
                 fire_frequency,  time_since_fire, source,        treatment_comment,
                 diversity_manipulated))%>%
  rbind(e245_metadata %>% mutate(dataset = "e245") %>%
          select(year,            site,            dataset,       plot,          higher_order_organization,
                 uniqueid,        temperature,     precipitation, growtemp, growprecip, treatment,
                 nutrients_added, nitrogen_amount, disturbance,   grazing,
                 fire_frequency,  time_since_fire, source,        treatment_comment,
                 diversity_manipulated))%>%
  rbind(e247_metadata %>% mutate(dataset = "e247") %>%
          select(year,            site,            dataset,       plot,          higher_order_organization,
                 uniqueid,        temperature,     precipitation, treatment, growtemp, growprecip,
                 nutrients_added, nitrogen_amount, disturbance,   grazing,
                 fire_frequency,  time_since_fire, source,        treatment_comment,
                 diversity_manipulated))


#SAVE####
write.csv(cdr_data,     paste(L1_dir, "CDR_plotlevel_metrics.csv", sep = "/"))
write.csv(cdr_sp_data,  paste(L1_dir, "CDR_specieslevel_abundance.csv", sep = "/"))
write.csv(cdr_metadata, paste(L1_dir, "CDR_metadata.csv", sep = "/"))
