# TITLE:        LTER Grassland Rock: KBS ANPP biomass and plant composition cleanup
# AUTHORS:      Caitlin Broderick
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


# bring in temp precip data
weatherdaily <- read.csv(file.path(L0_dir, "KBS/7-lter+weather+station+daily+precip+and+air+temp+1698756534_L0.csv"))
#weatherdaily <- read.csv("7-lter+weather+station+daily+precip+and+air+temp+1698756534.csv", stringsAsFactors = FALSE)

# bring in full weather data
weatherdaily_all <- read.csv(file.path(L0_dir, "KBS/12-lter+weather+station+daily+weather+all+variates+1707331925.csv"))


head(mcse)
head(micro)
head(glbrc)
head(glbrc_scaleup)
head(weatherdaily)


############################
############################
###### mcse 
############################
############################

unique(mcse$Treatment)
unique(mcse$Campaign)

# make column indicating that this is 291: main MCSE table
mcse$Source <- "MCSE (Table 291)"

mcse$Date <- lubridate::mdy(mcse$Date)


# filter to only be t7
t7 <- mcse %>% 
  filter(Treatment == "T7" & Campaign == "Peak Biomass")

# add columns
t7$Nutrients_added <- "no_fertilizer"

# We'll need to go through this more and see if we need to change or delete anything (like none plant stuff). 
# I wonder if we should make a master list of species between our 3 sites?
unique(t7$Species)
# REMOVE: "UnSorted" COULD BE MANY SP, "Standing Dead" and "Surface Litter" FROM LAST YEAR
# OK for sp richness and ANPP calcs:
# Genus sp., for example (""Hieracium sp. (*)")
# "Unknown dicot (*)"
# "Unknown grass" 
# "Dicots"
# "Monocots"
# "Moss"  

t7_nounknown <- t7 %>% 
  filter(Species != "UnSorted" & 
           Species != "Surface Litter" &
           Species != "Standing Dead")
# calculate total ANPP and richness: sum up ANPP for each plot, count rows for each plot
anpp_rich_t7 <- t7_nounknown %>% 
  group_by(Year, Treatment, Replicate, Station, Source, Nutrients_added) %>% 
  summarise(plot_biomass = sum(Biomass_g_m2), # add up biomass
            plot_richness = n()) # count number of rows for richness




############################
###### microplots in t7
############################

micro$Date <- lubridate::mdy(micro$Date)

micro$Year <- lubridate::year(micro$Date)

# make column indicating that this is 154: microplot from MCSE
micro$Source <- "MCSE Microplot (Table 154)"

# this is already summed across i think subplots (the 0.1 m2 samples) cuz total area sampled is 1
# anpp per treatment replicate disturbance regime....

names(micro)
unique(micro$Species)
#REMOVE:
"Unknown" 
"UnSorted"
"Surface Litter"
"Standing Dead"
# leaving in unknown grass/forb, only known to family/genus, etc.

# add columns
micro$Nutrients_added <- "N+"

micro_nounknown <- micro %>% 
  filter(Species != "Unknown" &  # should we take out unknown????
           Species != "UnSorted" &  # take out.
           Species != "Surface Litter" &  # take out
           Species != "Standing Dead") # take out
anpp_rich_micro <-  micro_nounknown %>% 
  group_by(Year, Treatment, Replicate, Disturbed_Microplot, Fertilized_Microplot,Source, Nutrients_added) %>% 
  summarise (plot_biomass = sum (Biomass_g_m2), # get plot biomass
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
unique(anpp_rich_KBS_T7$Treatment) # dope.

hist(anpp_rich_KBS_T7$plot_biomass)


# NEXT : GET Percent cover dataset (sp comp)

head(t7)
head(micro)


# merge main MCSE with microplot
t7_with_ANPP <- merge(t7_nounknown, anpp_rich_t7, by = c("Year", "Treatment", "Station",
                                          "Replicate", "Source", "Nutrients_added"))
head(t7_with_ANPP)

# get psesudo percent cover by dividing plant by total for ANPP...
t7_with_ANPP$Pseudo_PercCover <- t7_with_ANPP$Biomass_g_m2 / t7_with_ANPP$plot_biomass * 100
head(t7_with_ANPP)

# micro
micro_with_ANPP <- merge(micro_nounknown, anpp_rich_micro, by = c("Year", "Treatment", "Disturbed_Microplot", "Fertilized_Microplot",
                                                   "Replicate", "Source", "Nutrients_added"))
head(micro_with_ANPP)

micro_with_ANPP$Pseudo_PercCover <- micro_with_ANPP$Biomass_g_m2 / micro_with_ANPP$plot_biomass * 100

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
unique(glbrc$Treatment)
# https://lter.kbs.msu.edu/research/long-term-experiments/glbrc-intensive-experiment/
# G10 =  restored prairie
# G5 = Switchgrass
# G6 = Miscanthus
# G7 Native Grasses  - a mix of 4 species
# G9 = Early successional
unique(glbrc$Campaign)

# make column indicating that this is 291: main MCSE table
glbrc$Source <- "GLBRC (Table 269)"

glbrc$Date <- lubridate::mdy(glbrc$Date)

glbrc$Nutrients_added <- "no_fertilizer"
#this is a placeholder. check to make sure no fert. 

# filter to only be peak biomass, and trts we want, and year
# NOTE: Include native grassses ? 4 planted species. Did not here
glbrc_grassland <- glbrc %>% 
  filter( (Treatment == "G10" | # restored prairie
             Treatment == "G9") & # Early successional 
            Campaign == "peak biomass" &
            Year <2018) # 2018 on is not sorted yet
unique(glbrc_grassland$Treatment)


unique(glbrc_grassland$Species)
# remove: "UnSorted",   ,"Unknown" ,  "Standing Dead",  "Surface Litter"   

glbrc_grassland_nounknown <- glbrc_grassland %>% 
  filter(Species != "UnSorted" & 
           Species != "Unknown" & 
           Species != "Standing Dead" & 
           Species != "Surface Litter"  )
# calculate total ANPP: sum up ANPP for ALL species
anpp_rich_glbrc <- glbrc_grassland_nounknown %>% 
  group_by(Year, Treatment, Site, Replicate, Station,Source, Nutrients_added) %>% 
  summarise(plot_biomass = sum(Biomass_g_m2),
            plot_richness = n())

anpp_rich_glbrc



############################
###### glbrc scaleup
############################

names(glbrc_scaleup)
unique(glbrc_scaleup$Treatment)
#https://lter.kbs.msu.edu/maps/images/glbrc_scaleup_sites_history_and_naming_conventions.pdf
# M2 =  CRP-Prairie
# L3 = AGR-Prairie


# make column indicating what data table this is from KBS website
glbrc_scaleup$Source <- "GLBRC Scale-Up (Table 180)"

glbrc_scaleup$Date <- lubridate::mdy(glbrc_scaleup$Date)

glbrc_scaleup$Nutrients_added <- "no_fertilizer"
#this is a placeholder. check to make sure no fert. 


# filter to only be trts we caare about
# ALSO REMOVE 2009, there are only 2 species (Bromus, )
glbrc_scaleup_grassland <- glbrc_scaleup %>% 
  filter( (Treatment == "M2" | # CRP --> Prairie 
             Treatment == "L3" ) &  # AGR --> Prairie
             Fraction == "whole" & # the first year they were separating grain, crops...
           Year > 2009) # actually just remove that first year, looks like it was crops  
unique(glbrc_scaleup_grassland$Treatment)
unique(glbrc_scaleup_grassland$Year)
unique(glbrc_scaleup_grassland$Fraction)

unique(glbrc_scaleup_grassland$Species)
# remove"UnSorted" ,    "Surface Litter",     "Standing Dead"

glbrc_scaleup_grassland_nounknown <- glbrc_scaleup_grassland %>% 
  filter(Species != "UnSorted" &
           Species !=  "Surface Litter" &
           Species != "Standing Dead" )

# calculate total ANPP and rich: sum up ANPP for each plot, count rows
anpp_rich_glbrc_scaleup <- glbrc_scaleup_grassland_nounknown %>% 
  group_by(Year, Treatment,  Station,Source, Nutrients_added) %>% 
  summarise(plot_biomass = sum(Biomass_g_m2),
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
anpp_rich_glbrc_KBS <- rbind(anpp_rich_glbrc, anpp_rich_glbrc_scaleup)

head(anpp_rich_glbrc_KBS )
unique(anpp_rich_glbrc_KBS$Treatment) 

hist(anpp_rich_KBS_T7$plot_biomass)
hist(anpp_rich_glbrc_KBS$plot_biomass)


# NEXT : GET Percent cover dataset (sp comp) for GLBRC

head(glbrc_grassland_nounknown )
head(glbrc_scaleup_grassland_nounknown)


# BCSE - merge main witih scaleup
glbrc_BCSE_with_ANPP <- merge(glbrc_grassland_nounknown, anpp_rich_glbrc, by = c("Year", "Treatment", "Station", "Site",
                                                                  "Replicate", "Source", "Nutrients_added"))
head(glbrc_BCSE_with_ANPP)

# get psesudo percent cover by dividing plant by total for ANPP...
glbrc_BCSE_with_ANPP$Pseudo_PercCover <- glbrc_BCSE_with_ANPP$Biomass_g_m2 / glbrc_BCSE_with_ANPP$plot_biomass * 100
head(glbrc_BCSE_with_ANPP)



# Scaleup
glbrc_scaleup_with_ANPP <- merge(glbrc_scaleup_grassland_nounknown, anpp_rich_glbrc_scaleup, by = c("Year", "Treatment", "Station" , 
                                                                                     "Source", "Nutrients_added"))
head(glbrc_scaleup_with_ANPP)

glbrc_scaleup_with_ANPP$Pseudo_PercCover <- glbrc_scaleup_with_ANPP$Biomass_g_m2 / glbrc_scaleup_with_ANPP$plot_biomass * 100

head(glbrc_scaleup_with_ANPP)


# bind together the GLBRC datasets, they should include species percent cover
head(glbrc_BCSE_with_ANPP)
head(glbrc_scaleup_with_ANPP)
allGLBRC_SpComp <- dplyr:: bind_rows (glbrc_BCSE_with_ANPP, glbrc_scaleup_with_ANPP)

allGLBRC_SpComp






#






########################################
# bind all KBS dataset together???
########################################

# species comp
head(allt7_SpComp)
head(allGLBRC_SpComp)

allkbsdata_spcomp <- dplyr:: bind_rows (allt7_SpComp, allGLBRC_SpComp)
head(allkbsdata_spcomp)


# anpp only
head(anpp_rich_KBS_T7)
head(anpp_rich_glbrc_KBS)
allkbsdata_anpp <- dplyr:: bind_rows (anpp_rich_KBS_T7, anpp_rich_glbrc_KBS)



###################################
###################################
#then join up with temp and precip
###################################
###################################


head(weatherdaily)
weatherdaily$Date <- lubridate::mdy(weatherdaily$date)
head(weatherdaily)

weatheryear <- weatherdaily %>% 
  group_by(Year) %>% 
  summarise (meantemp =mean(air_temp_mean, na.rm = TRUE), # na.rm cuz missing obs for temp.
             annualprecip = sum(precipitation))

weatheryear


# Merge weather with sp comp data for both MCSE and GLBRC
allkbsdata_spcomp_tp <- merge (allkbsdata_spcomp, weatheryear, by = "Year")

allkbsdata_spcomp_tp



# add site column
# made this LTER_Site for now. because GLBRC hamde 
allkbsdata_spcomp_tp$lter_site <- "KBS"

names(allkbsdata_spcomp_tp) <- tolower(names(allkbsdata_spcomp_tp))

# write a new .csv with the cleaned and merged data and upload to the shared google drive L1 folder

write.csv(allkbsdata_spcomp_tp, file.path(L1_dir, "./KBS_MCSE_GLBRC_SpComp.csv"), row.names=F)
#Do not run: write.csv(allt7_SpComp_tp, "KBS_MCSE_T7_SpComp.csv")







allkbsdata_anpp_tp <- merge(allkbsdata_anpp , weatheryear , by = "Year")

# add site column
# made this LTER_Site for now. because GLBRC hamde 
allkbsdata_anpp_tp$lter_site <- "KBS"

names(allkbsdata_anpp_tp) <- tolower(names(allkbsdata_anpp_tp))


# write a new .csv with the cleaned and merged data and upload to the shared google drive L1 folder

write.csv(allkbsdata_anpp_tp, file.path(L1_dir, "./KBS_MCSE_GLBRC_ANPP_RICH.csv"), row.names=F)
#write.csv(allt7_ANPP_tp, "KBS_MCSE_T7_ANPP.csv")

ggplot(allkbsdata_anpp_tp, aes (x = annualprecip, y = plot_biomass)) + 
  geom_point()


ggplot(allkbsdata_anpp_tp, aes (x = annualprecip, y = plot_richness)) + 
  geom_point()










# do not think we need this bit at the end anymore.
###################################

###################################
#then join up with temp and precip

head(weatherdaily)
weatherdaily$Date <- lubridate::mdy(weatherdaily$date)
head(weatherdaily)

weatheryear <- weatherdaily %>% 
  group_by(Year) %>% 
  summarise (meantemp =mean(air_temp_mean, na.rm = TRUE), # na.rm cuz missing obs for temp.
             annualprecip = sum(precipitation))

weatheryear


allt7_SpComp_tp <- merge (allt7_SpComp, weatheryear, by = "Year")

allt7_SpComp_tp

# add site column
allt7_SpComp_tp$site <- "KBS"

names(allt7_SpComp_tp) <- tolower(names(allt7_SpComp_tp))

# write a new .csv with the cleaned and merged data and upload to the shared google drive L1 folder
write.csv(allt7_SpComp_tp, file.path(L1_dir, "./KBS_MCSE_T7_SpComp.csv"), row.names=F)
#write.csv(allt7_SpComp_tp, "KBS_MCSE_T7_SpComp.csv")



allt7_ANPP_tp <- merge(anpp_rich_KBS_T7, weatheryear, by = "Year")

allt7_ANPP_tp

# add site column
allt7_ANPP_tp$site <- "KBS"

names(allt7_ANPP_tp) <- tolower(names(allt7_ANPP_tp))

# write a new .csv with the cleaned and merged data and upload to the shared google drive L1 folder
write.csv(allt7_ANPP_tp, file.path(L1_dir, "./KBS_MCSE_T7_ANPP.csv"), row.names=F)
#write.csv(allt7_ANPP_tp, "KBS_MCSE_T7_ANPP.csv")


ggplot(allt7_ANPP_tp, aes (x = annualprecip, y = plot_biomass))
