# TITLE:        LTER Grassland Rock: KBS ANPP biomass and plant composition cleanup
# AUTHORS:      Caitlin Broderick and Ashley Darst
# COLLABORATORS:  
# DATA INPUT:   Data imported as csv files from shared Google drive L0 folder
# DATA OUTPUT:    
# PROJECT:      LTER Grassland Rock
# DATE:         October 2023

# this code is to clean MCSE (and microplot treatment) data, merge them, and link with temp/precip
# 1-17-24 - added GLBRC BCSE and scale-up code. 
# To do: check columns - what do we want to add? Also, we did not yet match up
# "Treatment" and "Replicate" columns, since dif sampling schemes, etc.

# Clear all existing data
rm(list=ls())

# Load packages
library(tidyverse)
library(gtools)
library(lubridate)

# Set working directory 
L0_dir <- Sys.getenv("L0DIR")
L1_dir <- Sys.getenv("L1DIR")
list.files(L0_dir)


# Read in CSV files
# bring in mcse data 
        # note: The leading number in the file name is the table code on the KBS LTER site. Additionally, KBS generally puts 
                # a lot of metadata rows at the top, and then has a weird system where they include units under the column names. 
                # Caitlin removed the metadata at the top and combined these (e.g so biomass_g_m2 was one column name instead of taking two rows).
                # Caitlin replaced these "straight-from KBS-website" files with the slightly modified ones on the shared Google drive. Essentially, 
                # the raw download data from KBS is not usable in R.

#bring in MCSE
mcse <- read.csv(file.path(L0_dir, "KBS/291-biomass+compilation+for+herbaceous+systems+1698767008_L0.csv"), stringsAsFactors = FALSE)
#mcse <- read.csv("291-biomass+compilation+for+herbaceous+systems+1698767008.csv", stringsAsFactors = FALSE)

# bring in microplot data
micro <- read.csv(file.path(L0_dir, "KBS/154-early+successional+microplot+biomass+sorted+to+species++1698756518_L0.csv"), stringsAsFactors = FALSE)
#micro <- read.csv("154-early+successional+microplot+biomass+sorted+to+species++1698756518.csv", stringsAsFactors = FALSE)

#bring in GLBRC 
glbrc <- read.csv(file.path(L0_dir, "KBS/269-above+ground+biomass+of+herbaceous+perennial+crops+treatments+5+6+7+9+10+1702908396_L0.csv"), stringsAsFactors = FALSE)


#bring in GLBRC scale-up
glbrc_scaleup <- read.csv(file.path(L0_dir, "KBS/180-biomass+of+the+glbrc+scale+up+fields+1702908404_L0.csv"), stringsAsFactors = FALSE)

# bring in NutNet (biomass, sp comp)
nutnet_cover <- read.csv(file.path(L0_dir, "KBS/KBS_NutNet_full_cover_data.csv"), stringsAsFactors = FALSE)
nutnet_bio <- read.csv(file.path(L0_dir, "KBS/KBS_NutNet_full_biomass_data.csv"), stringsAsFactors = FALSE)


# bring in temp precip data
weatherdaily <- read.csv(file.path(L0_dir, "KBS/7-lter+weather+station+daily+precip+and+air+temp+1698756534_L0.csv"))
#weatherdaily <- read.csv("7-lter+weather+station+daily+precip+and+air+temp+1698756534.csv", stringsAsFactors = FALSE)

# bring in full weather data
weatherdaily_all <- read.csv(file.path(L0_dir, "KBS/12-lter+weather+station+daily+weather+all+variates+1707331925.csv"))


names(mcse)
names(micro)
names(glbrc)
names(glbrc_scaleup)
names(nutnet_bio)
names(nutnet_cover)
names(weatherdaily)

# remove "x" column from nutnet files
nutnet_bio <- nutnet_bio[,-1]
nutnet_cover <- nutnet_cover[,-1]
names(nutnet_bio)
names(nutnet_cover)

names(mcse) <- tolower(names(mcse))
names(micro) <- tolower(names(micro))
names(glbrc) <- tolower(names(glbrc))
names(glbrc_scaleup) <- tolower(names(glbrc_scaleup))
names(nutnet_bio) <- tolower(names(nutnet_bio))
names(nutnet_cover) <- tolower(names(nutnet_cover))
names(weatherdaily) <- tolower(names(weatherdaily))


############################
############################
###### mcse 
############################
############################
names(mcse)
unique(mcse$treatment)
unique(mcse$campaign)

# make column indicating that this is 291: main MCSE table
mcse$source <- "MCSE (Table 291)"
#mcse$measurement_scale_biomass <- 1
#mcse$measurement_scale_cover <- 1
mcse$area_sampled_bio <- 1
mcse$area_sampled_cover <- 1
mcse$date <- lubridate::mdy(mcse$date)


# filter to only be t7
t7 <- mcse %>% 
  filter(treatment == "T7" & campaign == "Peak Biomass")




# add columns
t7$nutrients_added <- "no_fertilizer"
t7$nitrogen_amount <- NA
t7$disturbance <- "undisturbed"
t7$grazing <- "ungrazed"
t7$fire_frequency <- 1
t7$time_since_fire <- with(t7, ifelse(year == 2007, 2, 1)) # not really true, only a few months since fire # not burned in 2007
t7$experiment <- "mcse"
t7$cover_method <- "pseudo percent cover from mass"


# We'll need to go through this more and see if we need to change or delete anything (like none plant stuff). 
# I wonder if we should make a master list of species between our 3 sites?
unique(t7$species)
# REMOVE: "UnSorted" COULD BE MANY SP, "Standing Dead" and "Surface Litter" FROM LAST YEAR
# OK for sp richness and ANPP calcs:
# Genus sp., for example (""Hieracium sp. (*)")
# "Unknown dicot (*)"
# "Unknown grass" 
# "Dicots"
# "Monocots"
# "Moss"  

# Biomass categories for all datasets
genus_sp_in_biomass <- c("Hieracium sp. (*)",
                         "Desmodium sp.",
                         "Rumex species",
                         "Vitis sp.",
                         "Cornus sp.",
                         "Setaria sp. (*)",
                         "Monarda sp.",
                         "Malus spp. (*)",
                         "Festuca sp. (*)",
                         "Celastrus sp.",
                         "Rosa sp. (*)",
                         "Fragaria sp. (*)",
                         "Rubus sp. (*)",
                         "Trifolium sp. (*)",
                         "Lonicera spp.",
                         "Geum sp.",
                         "Crataegus spp.",
                         "Rhamnus sp.",
                         "Solidago sp. (*)",
                         "Melilotus sp.",
                         "Vicia sp.",
                         "Carex sp.",
                         "Solanum sp. (*)",
                         "Cirsium sp. (*)",
                         "Populus sp. (*)",
                         "Bromus sp. (*)",
                         "Aster sp. (*)",
                         "Plantago sp.",
                         "Acer spp.",
                         "Lactuca sp. (*)",
                         "Veronica sp.",
                         "Prunus spp. (*)",
                         "Cerastium sp.",
                         "Sonchus sp.",
                         "Poa sp. (*)",
                         "Diervilla sp.",
                         "Brassica sp. (*)",
                         "Juncus sp. (*)",
                         "Composite sp. (*)",
                         "Mimosa sp.",
                         "Erigeron sp.",
                         "Agrostis sp.",
                         "Galium sp. (*)",
                         "Lamium sp.",
                         "Lolium sp. (*)",
                         "Morus spp. (*)",
                         "Muhlenbergia species",
                         "Polygonum sp. (*)",
                         "Amaranthus species",
                         "Cenchrus sp.",
                         "Anthemis sp.",
                         "Apocynum sp. (*)",
                         "Crepis sp.",
                         "Digitaria sp",
                         "Epilobium sp.",
                         "Euphorbia sp.",
                         "Lathyrus sp.",
                         "Lepidium sp.",
                         "Lespedeza sp. (*)",
                         "Morrenia Sp.",
                         "Oxalis sp.",
                         "Panicum sp",
                         "Potentilla sp.",
                         "Quercus spp.",
                         "Ribes spp.",
                         "Silene sp.",
                         "Sorghum sp. (*)",
                         "Verbascum sp.",
                         "Viola sp. (*)",
                         "Rhus spp.",
                         "Selaginella sp.",
                         "ACER SP.",
                         "ALLIUM SP.",
                         "HIERACIUM SP.",
                         "JUNCUS SP.",
                         "LONICERA SP.",
                         "MALUS SP.",
                         "MELILOTUS SP.",
                         "PLANTAGO SP.",
                         "PRUNUS SP.",
                         "SETARIA SP.",
                         "TRIFOLIUM",
                         "ULMUS SP.",
                         "VERONICA SP.",
                         "VITIS SP. (kbs.us)"
                         )
non_plant_things_in_biomass <- c("Standing Dead",
                                 "Moss",
                                 "Lichen",
                                 "Surface Litter",
                                 "Fungi",
                                 "Leaf Litter",
                                 "BRYOPHYTE SP.",
                                 "GROUND",
                                 "OTHER ANIMAL DIGGING",
                                 "OTHER ANIMAL DROPPINGS",
                                 "OTHER LITTER",
                                 "OTHER WOODY OVERSTORY"
                                 )

maybe_plant_things_in_biomass <- c("UnSorted",
                                   "Unknown dicot (*)",
                                   "Unknown grass",
                                   "unknown Sedge",
                                   "Unkown Fabaceae",
                                   "Unknown monocot (*)",
                                   "unknown Brassicaceae",
                                   "Woody",
                                   "another unknown dicot",
                                   "unknown Asteraceae",
                                   "Unknown Orchidaceae",
                                   "Aster basal leaves",
                                   "Dicots",
                                   "Monocots",
                                   "Composite Basal Leaves",
                                   "Unknown",
                                   "unknown Elm",
                                   "another unknown woody species",
                                   "Dead Dicot",
                                   "Ferns",
                                   "Unknown Dicot",
                                   "Unknown Monocot",
                                   "Unknown Rosacae",
                                   "Unknown Solanaceae",
                                   "UNKNOWN",
                                   "UNKNOWN ASTERACEAE",
                                   "UNKNOWN BRASSICACEAE",
                                   "UNKNOWN GRASS",
                                   "UNKNOWN PTERIDOPHYTA"
                                   )

t7_nounknown <- t7 %>% 
  filter(!species %in% non_plant_things_in_biomass)


# check to make sure all grouping variables good - number of rows should be 1
t7_nounknown %>% 
  group_by(year, treatment, replicate, station, nutrients_added, nitrogen_amount,
           disturbance, grazing, fire_frequency, time_since_fire, species) %>% 
  count() # if there were >1, more than 1 sp entry for each "group"

# calculate total ANPP and richness: sum up ANPP for each plot, count rows for each plot
anpp_rich_t7 <- t7_nounknown %>% 
  group_by(year, treatment, replicate, station, experiment, nutrients_added, nitrogen_amount,
           disturbance, grazing, fire_frequency, time_since_fire, area_sampled_bio, area_sampled_cover, source) %>% 
  summarise(plot_biomass = sum(biomass_g_m2), # add up biomass
            plot_richness = n()) # count number of rows for richness




############################
###### microplots in t7
############################

# note that microplots had mutliple sampling events per season for some years
  # i will remove this based on sven graph

micro$date <- lubridate::mdy(micro$date)

micro$year <- lubridate::year(micro$date)

micro <- micro %>% 
  rename("area_sampled_bio" = "area_sampled")
micro$area_sampled_cover <- micro$area_sampled_bio


# this one will have month to filter by !!!!
micro$month <- lubridate::month(micro$date)

table(lubridate::month(micro$date))
table( micro$year, lubridate::month(micro$date))

# make column indicating that this is 154: microplot from MCSE
micro$source <- "MCSE Microplot (Table 154)"

# this is already summed across i think subplots (the 0.1 m2 samples) cuz total area sampled is 1
# anpp per treatment replicate disturbance regime....

names(micro)
unique(micro$species)
#REMOVE:
"Unknown" 
"UnSorted"
"Surface Litter"
"Standing Dead"
# leaving in unknown grass/forb, only known to family/genus, etc.

# add columns
micro$nutrients_added <- with(micro, ifelse(fertilized_microplot == "fertilized", "N", "no_fertilizer"))
micro$nitrogen_amount <- with(micro, ifelse(fertilized_microplot == "fertilized", 12.3, NA)) # g N m^-2 as granular ammonium nitrate 
micro$disturbance <- with(micro, ifelse(disturbed_microplot == "disturbed", "disturbed", "undisturbed"))
micro$grazing <- "ungrazed"
micro$fire_frequency <- 1
micro$time_since_fire <- with(micro, ifelse(year == 2007, 2, 1)) # not really true, only a few months since fire # not burned in 2007
micro$experiment <- "microplots"
micro$cover_method <- "pseudo percent cover from mass"


micro_nounknown <- micro %>% 
  filter(!species %in% non_plant_things_in_biomass)

# check to make sure all grouping variables good - number of rows should be 1
micro_nounknown %>% 
  group_by(year, treatment, replicate,  nutrients_added, nitrogen_amount,
           disturbance, grazing, fire_frequency, time_since_fire, species) %>% 
  count()  # there are >1 : there are repeats
  # species found in different sampling months in the same year same plot....

# remove sampling events that occurred at nonpeak times
  # sven email 3-6-2024
micro_norepeatsampling_a <-micro_nounknown %>% 
  filter (year> 1989) # took out 1989 not all trts implemented yet
micro_norepeatsampling_b <- micro_norepeatsampling_a %>% 
  filter(  !(year == 1991 & month == 7)) #remove july sampling 1991 keep august
micro_norepeatsampling_c <- micro_norepeatsampling_b %>% 
  filter(  !(year == 1993 & month == 7))#remove july sampling 1993 keep august
micro_clean <- micro_norepeatsampling_c %>% 
  filter(  !(year == 1996 & month == 6)) #remove july sampling 1996 keep august


# now lets check
micro_clean %>% 
  group_by(year, treatment, replicate,  nutrients_added, nitrogen_amount,
           disturbance, grazing, fire_frequency, time_since_fire, species) %>% 
  count() # just down to the unknowns now. 
  # not sure where to go from here. maybe just leave it?
  # make decisions later?


table( micro$year, micro$month)
table( micro_clean$year, micro_clean$month)
# un/disturbed plots often harvested in different months



anpp_rich_micro <-  micro_clean %>% 
  group_by(year, treatment, replicate, experiment, nutrients_added,
           nitrogen_amount, disturbance, grazing, fire_frequency, time_since_fire, area_sampled_bio, area_sampled_cover, source) %>% 
  summarise (plot_biomass = sum (biomass_g_m2), # get plot biomass
             plot_richness = n()) # get plot richness

head(anpp_rich_micro)






#########################################################################
############ MCSE: combine these datasets into ANPP/rich df and SPCOMP df
########################################################################


#FIRST:
# combine these ANPP datasets for microplot and normal Mcse
names(anpp_rich_t7)
names(anpp_rich_micro)


# merge main MCSE with microplot
anpp_rich_KBS_T7 <- rbind(anpp_rich_t7, anpp_rich_micro)

head(anpp_rich_KBS_T7 )
unique(anpp_rich_KBS_T7$treatment) # dope.

hist(anpp_rich_KBS_T7$plot_biomass)


# NEXT : GET Percent cover dataset (sp comp)

head(t7)
head(micro)


# merge main MCSE with microplot
t7_with_ANPP <- merge(t7_nounknown, anpp_rich_t7, by = c("year", "treatment", "station",
                                          "replicate", "experiment", "nutrients_added",
                                          "nitrogen_amount", "disturbance", "grazing",
                                          "fire_frequency", "time_since_fire", "area_sampled_bio", "area_sampled_cover", "source"))
head(t7_with_ANPP)

# get psesudo percent cover by dividing plant by total for ANPP...
t7_with_ANPP$relative_abundance <- t7_with_ANPP$biomass_g_m2 / t7_with_ANPP$plot_biomass
head(t7_with_ANPP)

# micro
micro_with_ANPP <- merge(micro_clean, anpp_rich_micro, by = c("year", "treatment",
                                                   "replicate", "experiment", "nutrients_added",
                                                   "nitrogen_amount", "disturbance", "grazing",
                                                   "fire_frequency", "time_since_fire", "area_sampled_bio", "area_sampled_cover", "source"))
head(micro_with_ANPP)

micro_with_ANPP$relative_abundance <- micro_with_ANPP$biomass_g_m2 / micro_with_ANPP$plot_biomass

head(micro_with_ANPP)


# bind together
head(t7_with_ANPP)
head(micro_with_ANPP)
allt7_SpComp <- dplyr:: bind_rows (t7_with_ANPP, micro_with_ANPP)

allt7_SpComp



############################
############################
###### glbrc
############################
############################

names(glbrc)
unique(glbrc$treatment)
# https://lter.kbs.msu.edu/research/long-term-experiments/glbrc-intensive-experiment/
# G10 =  restored prairie
# G5 = Switchgrass
# G6 = Miscanthus
# G7 Native Grasses  - a mix of 4 species
# G9 = Early successional
unique(glbrc$campaign)

table(glbrc$area_sampled_m2)




# make column indicating that this is 291: main MCSE table
glbrc$source <- "GLBRC (Table 269)"

glbrc$date <- lubridate::mdy(glbrc$date)

glbrc$nutrients_added <- "no_fertilizer"

#this is a placeholder. check to make sure no fert. 
glbrc <- glbrc %>% 
  rename("area_sampled_bio" = "area_sampled_m2")
glbrc$area_sampled_cover <- glbrc$area_sampled_bio

glbrc <- glbrc %>%  # rename "site" column, need to use that for later.
  rename(  "glbrc_site" ="site")
names(glbrc)

# filter to only be peak biomass, and trts we want, and year
# NOTE: Include native grassses ? 4 planted species. Did not here
glbrc_grassland <- glbrc %>% 
  filter( (treatment == "G10" | # restored prairie
             treatment == "G9") & # Early successional 
            campaign == "peak biomass" &
            year <2018) # 2018 on is not sorted yet
unique(glbrc_grassland$treatment)

glbrc_grassland$nitrogen_amount <- NA # are subplots included???
glbrc_grassland$disturbance <- "undisturbed"
glbrc_grassland$grazing <- "ungrazed"
glbrc_grassland$fire_frequency <- NA
glbrc_grassland$time_since_fire <- NA
glbrc_grassland$experiment <- "glbrc"
glbrc_grassland$cover_method <- "pseudo percent cover from mass"

unique(glbrc_grassland$species)
# remove: "UnSorted",   ,"Unknown" ,  "Standing Dead",  "Surface Litter"   

glbrc_grassland_nounknown <- glbrc_grassland %>% 
  filter(!species %in% non_plant_things_in_biomass)

# Remove observations when subplots were lumped together for harvest
glbrc_grassland_nounknown <- glbrc_grassland_nounknown %>%
  filter(station != "")

# calculate total ANPP: sum up ANPP for ALL species
anpp_rich_glbrc <- glbrc_grassland_nounknown %>% 
  group_by(year, treatment, glbrc_site, replicate, station,experiment, nutrients_added, area_sampled_bio,area_sampled_cover,
           nitrogen_amount, disturbance, grazing, fire_frequency, time_since_fire, source) %>% 
  summarise(plot_biomass = sum(biomass_g_m2),
            plot_richness = n())

anpp_rich_glbrc



############################
###### glbrc scaleup
############################

names(glbrc_scaleup)
unique(glbrc_scaleup$treatment)
#https://lter.kbs.msu.edu/maps/images/glbrc_scaleup_sites_history_and_naming_conventions.pdf
# M2 =  CRP-Prairie
# L3 = AGR-Prairie


# make column indicating what data table this is from KBS website
glbrc_scaleup$source <- "GLBRC Scale-Up (Table 180)"

glbrc_scaleup$date <- lubridate::mdy(glbrc_scaleup$date)

glbrc_scaleup$nutrients_added <- "no_fertilizer"
#this is a placeholder. check to make sure no fert. 

glbrc_scaleup <- glbrc_scaleup %>%  # rename "site" column, need to use that for later.
  rename(  "glbrc_site" ="site")

# Add area sampled information based on https://data.sustainability.glbrc.org/protocols/156
glbrc_scaleup$area_sampled_bio <- 1
glbrc_scaleup$area_sampled_cover <- 1

# filter to only be trts we caare about
# ALSO REMOVE 2009, there are only 2 species (Bromus, )
glbrc_scaleup_grassland <- glbrc_scaleup %>% 
  filter( (treatment == "M2" | # CRP --> Prairie 
             treatment == "L3" ) &  # AGR --> Prairie
             fraction == "whole" & # the first year they were separating grain, crops...
           year > 2009) # actually just remove that first year, looks like it was crops  
unique(glbrc_scaleup_grassland$treatment)
unique(glbrc_scaleup_grassland$year)
unique(glbrc_scaleup_grassland$fraction)

glbrc_scaleup_grassland$nitrogen_amount <- NA # has been fertilized twice, include?
glbrc_scaleup_grassland$disturbance <- "undisturbed"
glbrc_scaleup_grassland$grazing <- "ungrazed"
glbrc_scaleup_grassland$fire_frequency <- NA
glbrc_scaleup_grassland$time_since_fire <- NA
glbrc_scaleup_grassland$experiment <- "glbrc_scaleup"
glbrc_scaleup_grassland$cover_method <- "pseudo percent cover from mass"

unique(glbrc_scaleup_grassland$species)
# remove"UnSorted" ,    "Surface Litter",     "Standing Dead"

glbrc_scaleup_grassland_nounknown <- glbrc_scaleup_grassland %>% 
  filter(!species %in% non_plant_things_in_biomass)


# calculate total ANPP and rich: sum up ANPP for each plot, count rows
names(glbrc_scaleup_grassland_nounknown)
anpp_rich_glbrc_scaleup <- glbrc_scaleup_grassland_nounknown %>% 
  group_by(year, treatment,  station,glbrc_site, experiment, nutrients_added,
           nitrogen_amount, disturbance, grazing, fire_frequency, time_since_fire, source, area_sampled_cover, area_sampled_bio) %>% 
  summarise(plot_biomass = sum(biomass_g_m2),
            plot_richness = n()) # SOMETHING GOING ON IN 2009!

anpp_rich_glbrc_scaleup




####################################################################
########### GLBRC: combine these datasets into ANPP df and SPCOMP df
####################################################################


#FIRST:
# combine these ANPP datasets
names(anpp_rich_glbrc)
names(anpp_rich_glbrc_scaleup)


# merge main BCSE with scale up (ANPP)
anpp_rich_glbrc_KBS <- bind_rows(anpp_rich_glbrc, anpp_rich_glbrc_scaleup)
head(anpp_rich_glbrc_KBS )
unique(anpp_rich_glbrc_KBS$treatment) 

hist(anpp_rich_KBS_T7$plot_biomass)
hist(anpp_rich_glbrc_KBS$plot_biomass)


# NEXT : GET Percent cover dataset (sp comp) for GLBRC

head(glbrc_grassland_nounknown )
head(glbrc_scaleup_grassland_nounknown)


# BCSE - merge main witih scaleup
glbrc_BCSE_with_ANPP <- merge(glbrc_grassland_nounknown, anpp_rich_glbrc, by = c("year", "treatment", "station", "glbrc_site", "area_sampled_bio", "area_sampled_cover",
                                                                  "replicate", "experiment", "nutrients_added",
                                                                  "nitrogen_amount", "disturbance", "grazing",
                                                                  "fire_frequency", "time_since_fire", "source"))
head(glbrc_BCSE_with_ANPP)

# get psesudo percent cover by dividing plant by total for ANPP...
glbrc_BCSE_with_ANPP$relative_abundance <- glbrc_BCSE_with_ANPP$biomass_g_m2 / glbrc_BCSE_with_ANPP$plot_biomass
head(glbrc_BCSE_with_ANPP)
hist(glbrc_BCSE_with_ANPP$relative_abundance)


# Scaleup
glbrc_scaleup_with_ANPP <- merge(glbrc_scaleup_grassland_nounknown, anpp_rich_glbrc_scaleup, by = c("year", "treatment", "station" , "glbrc_site", 
                                                                                     "experiment", "nutrients_added", "nitrogen_amount", "disturbance", "grazing",
                                                                                     "fire_frequency", "time_since_fire", "source"))
names(glbrc_scaleup_with_ANPP)
head(glbrc_scaleup_with_ANPP)

glbrc_scaleup_with_ANPP$relative_abundance <- glbrc_scaleup_with_ANPP$biomass_g_m2 / glbrc_scaleup_with_ANPP$plot_biomass

head(glbrc_scaleup_with_ANPP)
hist(glbrc_scaleup_with_ANPP$relative_abundance)

# bind together the GLBRC datasets, they should include species percent cover
head(glbrc_BCSE_with_ANPP)
head(glbrc_scaleup_with_ANPP)
allGLBRC_SpComp <- dplyr:: bind_rows (glbrc_BCSE_with_ANPP, glbrc_scaleup_with_ANPP)

allGLBRC_SpComp
names(allGLBRC_SpComp)





############################
############################
###### nutnet
############################
############################


# ANPP
names(nutnet_bio)
unique (nutnet_bio$year_trt) # how many years has the trt been applied
unique (nutnet_bio$trt) # rename this to treatment
nutnet_bio$area_sampled_bio <- 0.2
# column management 
nutnet_bio <- nutnet_bio %>% rename (treatment = trt) # rename to match other datasets
nutnet_bio$source <- "NutNet (Not publicly available)" # make column indicating where this came from
  unique (nutnet_bio$treatment)
  unique(micro$fertilized_microplot)
nutnet_bio$nutrients_added <-  dplyr::case_match(nutnet_bio$treatment, "Fence" ~"no_fertilizer", "Control"~"no_fertilizer" ,.default = nutnet_bio$treatment)
  unique(nutnet_bio$nutrients_added)
nutnet_bio$nitrogen_amount <- with(nutnet_bio, 
                                   ifelse(nutrients_added == "NK" | nutrients_added == "N" | nutrients_added == "NP" | nutrients_added == "NPK" | nutrients_added == "NPK+Fence", 10, NA)) # g N m^-2 as granular ammonium nitrate 
  unique(nutnet_bio$nitrogen_amount) # only amount nitrogen # should we add columns for other elements? (each is 10)
nutnet_bio$area_sampled_bio
  
#micro$disturbance <- with(micro, ifelse(disturbed_microplot == "disturbed", "disturbed", "undisturbed"))
nutnet_bio$grazing <- "ungrazed"
nutnet_bio$disturbance <- "NA" # CHECK THIS!!!! what to do with fence !!!!!
nutnet_bio$fire_frequency <- NA
nutnet_bio$time_since_fire <- NA
nutnet_bio$experiment <- "nutnet"


# remove dead litter
unique(nutnet_bio$category)
nutnet_bio_nolitter <- nutnet_bio %>% 
  filter(category != "LITTER"  )
# calculate total ANPP: sum up ANPP for vasc an non vasc plants
anpp_nutnet <- nutnet_bio_nolitter %>% 
  group_by(year, treatment,   block, plot, subplot,experiment, nutrients_added, area_sampled_bio,
           nitrogen_amount, disturbance, grazing, fire_frequency, time_since_fire, source) %>% 
  summarise(plot_biomass = sum(mass))

anpp_nutnet





# richness
names(nutnet_cover)
unique (nutnet_cover$year_trt) # how many years has the trt been applied
unique (nutnet_cover$trt) # rename this to treatment

nutnet_cover$area_sampled_cover <- 1
# column management 
nutnet_cover <- nutnet_cover %>% rename (treatment = trt) # rename to match other datasets
nutnet_cover <- nutnet_cover %>% rename (species = taxon)
nutnet_cover$source <- "NutNet (Not publicly available)" # make column indicating where this came from
unique (nutnet_cover$treatment)
unique(micro$fertilized_microplot)
nutnet_cover$nutrients_added <-  dplyr::case_match(nutnet_cover$treatment, "Fence" ~"no_fertilizer", "Control"~"no_fertilizer" ,.default = nutnet_cover$treatment)
unique(nutnet_cover$nutrients_added)
nutnet_cover$nitrogen_amount <- with(nutnet_cover, 
                                     ifelse(nutrients_added == "NK" | nutrients_added == "N" | nutrients_added == "NP" | nutrients_added == "NPK" | nutrients_added == "NPK+Fence", 10, NA)) # g N m^-2 as granular ammonium nitrate 
unique(nutnet_cover$nitrogen_amount) # only amount nitrogen # should we add columns for other elements? (each is 10)
#micro$disturbance <- with(micro, ifelse(disturbed_microplot == "disturbed", "disturbed", "undisturbed"))
nutnet_cover$grazing <- "ungrazed"
nutnet_cover$disturbance <- "NA" # CHECK THIS!!!! what to do with fence !!!!!
nutnet_cover$fire_frequency <- NA
nutnet_cover$time_since_fire <- NA
nutnet_cover$experiment <- "nutnet"

nutnet_cover$cover_method <- "percent cover to nearest 1 percent"

# remove dead litter and unknown
unique(nutnet_cover$species)
nutnet_cover_nounknown <- nutnet_cover %>% 
  filter(!species %in% non_plant_things_in_biomass)

names(nutnet_cover_nounknown )

nutnet_cover_nounknown <- nutnet_cover_nounknown %>% 
  group_by(year, treatment,   block, plot, subplot, experiment, nutrients_added, area_sampled_cover,
           nitrogen_amount, disturbance, grazing, fire_frequency, time_since_fire, source) %>%
  mutate (relative_abundance = max_cover / sum(max_cover))

# calculate richness
rich_nutnet <- nutnet_cover_nounknown %>% 
  group_by(year, treatment,   block, plot, subplot,experiment, nutrients_added, area_sampled_cover,
           nitrogen_amount, disturbance, grazing, fire_frequency, time_since_fire, source) %>% 
  summarise(plot_richness = n())

names(rich_nutnet)
names(anpp_nutnet)

anpp_rich_nutnet <- merge (rich_nutnet, anpp_nutnet, 
                           by = c("year", "treatment",  "block"  , "plot"  ,         
                            "subplot"  , "experiment" ,  "nutrients_added", "nitrogen_amount", "disturbance" ,   
                             "grazing",  "fire_frequency"  ,"time_since_fire", "source"))

# make a species cover file that includes plot biomass, rich, anpp....
nutnet_cover_with_ANPP <- merge(nutnet_cover_nounknown, anpp_rich_nutnet,
                                by = c("year", "treatment", "block"  , "plot"  ,         
                                      "subplot"  , "experiment" ,  "nutrients_added", "nitrogen_amount", "disturbance" ,   
                                      "grazing",  "fire_frequency"  ,"time_since_fire", "source"))







########################################
# bind all KBS dataset together???
########################################

# species comp
head(allt7_SpComp)
head(allGLBRC_SpComp)
head(nutnet_cover_with_ANPP)

allkbsdata_spcomp <- dplyr:: bind_rows (allt7_SpComp, allGLBRC_SpComp, nutnet_cover_with_ANPP)
head(allkbsdata_spcomp)


# anpp only
head(anpp_rich_KBS_T7)
head(anpp_rich_glbrc_KBS)
head(anpp_rich_nutnet)
allkbsdata_anpp <- dplyr:: bind_rows (anpp_rich_KBS_T7, anpp_rich_glbrc_KBS, anpp_rich_nutnet)

names(allkbsdata_anpp)

###################################
###################################
#then join up with temp and precip
###################################
###################################


head(weatherdaily)
weatherdaily$date <- lubridate::mdy(weatherdaily$date)
head(weatherdaily)

weatheryear <- weatherdaily %>% 
  group_by(year) %>% 
  summarise (meantemp =mean(air_temp_mean, na.rm = TRUE), # na.rm cuz missing obs for temp.
             annualprecip = sum(precipitation))

weatheryear

# calculate growing season precip (85-248)
# need to add this in

weathergrow <- weatherdaily
weathergrow$yday <- yday(weathergrow$date)
weathergrowyear <- weathergrow %>% 
  group_by(year) %>% 
  filter(yday %in% (85:248)) %>%
  summarise (growtemp =mean(air_temp_mean, na.rm = TRUE), # na.rm cuz missing obs for temp.
             growprecip = sum(precipitation))

allweather_year <- merge (weatheryear, weathergrowyear, by = "year")


# Merge weather with sp comp data for both MCSE and GLBRC
allkbsdata_spcomp_tp <- merge (allkbsdata_spcomp, allweather_year, by = "year")

allkbsdata_spcomp_tp

names(allkbsdata_spcomp_tp)

# add site column
# made this LTER_Site for now. because GLBRC hamde 
allkbsdata_spcomp_tp$site <- "KBS"

idcols <- c("site", "experiment","treatment", "station", "replicate" ,"block", "plot", "subplot", "nutrients_added", "disturbance")
plotid <- c("site", "experiment", "treatment", "station", "replicate" ,"block", "plot", "subplot","nutrients_added", "disturbance")
higher_level_organization <-  c("site", "experiment", "treatment", "station", "replicate" ,"block")


allkbsdata_spcomp_tp$unique_id <- apply( allkbsdata_spcomp_tp[ , idcols ] , 1 , paste , collapse = "_" )
allkbsdata_spcomp_tp$plot_id <- apply( allkbsdata_spcomp_tp[ , plotid ] , 1 , paste , collapse = "_" )
allkbsdata_spcomp_tp$higher_level_organization <- apply( allkbsdata_spcomp_tp[ , higher_level_organization ] , 1 , paste , collapse = "_" )
allkbsdata_spcomp_tp$unique_id
allkbsdata_spcomp_tp$plot_id
allkbsdata_spcomp_tp$higher_level_organization
# write a new .csv with the cleaned and merged data and upload to the shared google drive L1 folder

allkbsdata_spcomp_tp$unique_id <- str_remove_all(allkbsdata_spcomp_tp$unique_id, "_NA")
allkbsdata_spcomp_tp$plot_id<- str_remove_all(allkbsdata_spcomp_tp$plot_id, "_NA")
allkbsdata_spcomp_tp$higher_level_organization<- str_remove_all(allkbsdata_spcomp_tp$higher_level_organization, "_NA")

allkbsdata_spcomp_tp$unique_id
allkbsdata_spcomp_tp$plot_id
allkbsdata_spcomp_tp$higher_level_organization

#write.csv(allkbsdata_spcomp_tp, file.path(L1_dir, "./KBS_MCSE_GLBRC_NUTNET_SpComp.csv"), row.names=F)
#Do not run: write.csv(allt7_SpComp_tp, "KBS_MCSE_T7_SpComp.csv")







allkbsdata_anpp_tp <- merge(allkbsdata_anpp , allweather_year , by = "year")

# add site column
# made this LTER_Site for now. because GLBRC hamde 
allkbsdata_anpp_tp$site <- "KBS"


allkbsdata_anpp_tp$unique_id <- apply(allkbsdata_anpp_tp[ , idcols ], 1 , paste , collapse = "_" )
allkbsdata_anpp_tp$plot_id <- apply(allkbsdata_anpp_tp[ , plotid ], 1 , paste , collapse = "_" )
allkbsdata_anpp_tp$higher_level_organization <- apply(allkbsdata_anpp_tp[ , higher_level_organization ], 1 , paste , collapse = "_" )


allkbsdata_anpp_tp$unique_id
allkbsdata_anpp_tp$plot_id 
allkbsdata_anpp_tp$higher_level_organization 

allkbsdata_anpp_tp$unique_id <- str_remove_all(allkbsdata_anpp_tp$unique_id, "_NA")
allkbsdata_anpp_tp$plot_id<- str_remove_all(allkbsdata_anpp_tp$plot_id, "_NA")
allkbsdata_anpp_tp$higher_level_organization <- str_remove_all(allkbsdata_anpp_tp$higher_level_organization , "_NA")

allkbsdata_anpp_tp$unique_id
allkbsdata_anpp_tp$plot_id
allkbsdata_anpp_tp$higher_level_organization 

# write a new .csv with the cleaned and merged data and upload to the shared google drive L1 folder

#write.csv(allkbsdata_anpp_tp, file.path(L1_dir, "./KBS_MCSE_GLBRC_NUTNET_ANPP_RICH.csv"), row.names=F)
#write.csv(allt7_ANPP_tp, "KBS_MCSE_T7_ANPP.csv")

ggplot(allkbsdata_anpp_tp, aes (x = annualprecip, y = plot_biomass)) + 
  geom_point()


ggplot(allkbsdata_anpp_tp, aes (x = annualprecip, y = plot_richness)) + 
  geom_point()


ggplot(allkbsdata_anpp_tp, aes (x = growprecip, y = plot_richness)) + 
  geom_point()




# 6-17-2024:



# first i guess we calculate diversity. 

# calculate shannon diversity   - based on ashleys code
head(allkbsdata_spcomp_tp)
kbs_species <- allkbsdata_spcomp_tp %>%
  mutate(proportions = (relative_abundance/100)*(log(relative_abundance/100)))
kbs_species$proportions[is.na(kbs_species$proportions)] <- 0  # Fix Nas to zeros
kbs_shannon <- kbs_species %>%
  #group_by(year, treatment, station, replicate, disturbance, nutrients_added) %>%
  group_by(unique_id) %>% 
  summarise(shannon = -1*sum(proportions)) # single value for each plot
kbs_shannon

names(kbs_shannon)

#calculate evenness


kbs_plot_level_metrics <- merge ( allkbsdata_anpp_tp, kbs_shannon, by = "unique_id")
kbs_plot_level_metrics$plot_evenness <- kbs_plot_level_metrics$shannon / log (kbs_plot_level_metrics$plot_richness)

names(kbs_plot_level_metrics)


# here, we are doing two things:
  # (1) remove microplots - inconsistent sampling
  # (2) breaking up into main data and metadata
# for each of two datasets:
  # allkbsdata_spcomp_tp
  # allkbsdata_anpp_tp


# make new cols 

# subset cols for c



cols_plot_level_metrics <- c("year", "site","higher_level_organization" , "plot_id","unique_id", "plot_biomass", "plot_richness", "plot_evenness", "area_sampled_bio" , "area_sampled_cover", "shannon", "source")



kbs_plot_level_metrics_sel <- kbs_plot_level_metrics %>% select(all_of (cols_plot_level_metrics ))

head(kbs_plot_level_metrics_sel)

# match columns names with CDR and KNZ
kbs_plot_level_metrics_sel <- rename(kbs_plot_level_metrics_sel, plot = plot_id)
kbs_plot_level_metrics_sel <- rename(kbs_plot_level_metrics_sel, higher_order_organization = higher_level_organization)
kbs_plot_level_metrics_sel <- rename(kbs_plot_level_metrics_sel, uniqueid = unique_id)
kbs_plot_level_metrics_sel <- rename(kbs_plot_level_metrics_sel, richness = plot_richness)
kbs_plot_level_metrics_sel <- rename(kbs_plot_level_metrics_sel, evenness = plot_evenness)
kbs_plot_level_metrics_sel <- rename(kbs_plot_level_metrics_sel, measurement_scale_biomass = area_sampled_bio)
kbs_plot_level_metrics_sel <- rename(kbs_plot_level_metrics_sel, measurement_scale_cover = area_sampled_cover)


write.csv(kbs_plot_level_metrics_sel, file.path(L1_dir, "./KBS_plot_level_metrics.csv"), row.names=F)


cols_metadata <- c("year", "site", "experiment", "plot_id", "higher_level_organization" , "unique_id", 
                             "meantemp", "annualprecip", "growtemp", "growprecip", "treatment", "nutrients_added", 
                             "nitrogen_amount", "disturbance", "grazing", "fire_frequency", "time_since_fire", "source")


kbs_meta <- kbs_plot_level_metrics %>% select(all_of (cols_metadata ))
# match columns names with CDR and KNZ
kbs_meta <- rename(kbs_meta, plot = plot_id)
kbs_meta <- rename(kbs_meta, higher_order_organization = higher_level_organization)
kbs_meta <- rename(kbs_meta, uniqueid = unique_id)

# Convert experiment specific treatment values into control
controls <- c("G9", "G10", "Control", "L3", "M2")
kbs_meta$treatment <- replace(kbs_meta$treatment, kbs_meta$treatment %in% controls, "control")

kbs_meta %>%
  filter(treatment == "T7" & disturbance == "disturbed")

kbs_meta$treatment[kbs_meta$treatment == "T7" & kbs_meta$nutrients_added == "N" & kbs_meta$disturbance == "undisturbed"] <- "N"
kbs_meta$treatment[kbs_meta$treatment == "T7" & kbs_meta$nutrients_added == "N" & kbs_meta$disturbance == "disturbed"] <- "N + tilled"
kbs_meta$treatment[kbs_meta$treatment == "T7" & kbs_meta$nutrients_added != "N" & kbs_meta$disturbance == "disturbed"] <- "tilled"


write.csv(kbs_meta, file.path(L1_dir, "./KBS_metadata.csv"), row.names=F)


names(allkbsdata_spcomp_tp)

cols_species_level_abundance <- c("year", "site","higher_level_organization" , "plot_id","unique_id", "species" , "relative_abundance" , "cover_method", "area_sampled_bio" , "area_sampled_cover", "biomass_g_m2")
kbs_species_level_abundance <- allkbsdata_spcomp_tp %>% select(all_of (cols_species_level_abundance ))

kbs_species_level_abundance
names(kbs_species_level_abundance)
head(kbs_species_level_abundance)
# match columns names with CDR and KNZ
kbs_species_level_abundance <- rename(kbs_species_level_abundance, plot = plot_id)
kbs_species_level_abundance <- rename(kbs_species_level_abundance, higher_order_organization = higher_level_organization)
kbs_species_level_abundance <- rename(kbs_species_level_abundance, uniqueid = unique_id)
kbs_species_level_abundance <- rename(kbs_species_level_abundance, abundance = biomass_g_m2)
kbs_species_level_abundance$original_measurement_unit <- "g_m2"
kbs_species_level_abundance <- kbs_species_level_abundance %>%
  select(-c(area_sampled_bio, area_sampled_cover))

write.csv(kbs_species_level_abundance, file.path(L1_dir, "./KBS_species_level_abundance.csv"), row.names=F)

