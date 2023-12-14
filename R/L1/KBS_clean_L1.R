# TITLE:          LTER Grassland Rock: KBS ANPP biomass and plant composition cleanup
# AUTHORS:        Caitlin Broderick
# COLLABORATORS:  
# DATA INPUT:     Data imported as csv files from shared Google drive L0 folder
# DATA OUTPUT:    
# PROJECT:        LTER Grassland Rock
# DATE:  October 2023

# this code is to clean MCSE (and microplot treatment) data, merge them, and link with temp/precip
# TO DO: 

# Clear all existing data
rm(list=ls())

# Load packages
library(tidyverse)
library(gtools)

# Set working directory 
L0_dir <- Sys.getenv("L0DIR")
L1_dir <- Sys.getenv("L1DIR")
list.files(L0_dir)

#setwd("~/Downloads")

# Read in CSV files
# bring in mcse data 
        # note: The leading number in the file name is the table code on the KBS LTER site. Additionally, KBS generally puts 
                # a lot of metadata rows at the top, and then has a weird system where they include units under the column names. 
                # Caitlin removed the metadata at the top and combined these (e.g so biomass_g_m2 was one column name instead of taking two rows).
                # Caitlin replaced these "straight-from KBS-website" files with the slightly modified ones on the shared Google drive. Essentially, 
                # the raw download data from KBS is not usable in R.

mcse <- read.csv(file.path(L0_dir, "KBS/291-biomass+compilation+for+herbaceous+systems+1698767008_L0.csv"), stringsAsFactors = FALSE)
#mcse <- read.csv("291-biomass+compilation+for+herbaceous+systems+1698767008.csv", stringsAsFactors = FALSE)

# bring in microplot data
micro <- read.csv(file.path(L0_dir, "KBS/154-early+successional+microplot+biomass+sorted+to+species++1698756518_L0.csv"), stringsAsFactors = FALSE)
#micro <- read.csv("154-early+successional+microplot+biomass+sorted+to+species++1698756518.csv", stringsAsFactors = FALSE)

# bring in temp precip data
weatherdaily <- read.csv(file.path(L0_dir, "KBS/7-lter+weather+station+daily+precip+and+air+temp+1698756534_L0.csv"))
#weatherdaily <- read.csv("7-lter+weather+station+daily+precip+and+air+temp+1698756534.csv", stringsAsFactors = FALSE)

head(mcse)
head(micro)
head(weatherdaily)


############################
###### mcse 
############################

unique(mcse$Treatment)
unique(mcse$Campaign)

# make column indicating that this is 291: main MCSE table
mcse$Source <- "MCSE (Table 291)"

mcse$Date <- lubridate::mdy(mcse$Date)


# filter to only be t7
t7 <- mcse %>% 
  filter(Treatment == "T7" & Campaign == "Peak Biomass")


# calculate total ANPP: sum up ANPP for each species
anpp_t7 <- t7 %>% 
  group_by(Year, Treatment, Replicate, Station,Source) %>% 
  summarise(anpp = sum(biomass_g_m2))

anpp_t7



############################
###### microplots in t7
############################

micro$Date <- lubridate::mdy(micro$Sample_date)

micro$Year <- lubridate::year(micro$Date)

# make column indicating that this is 154: microplot from MCSE
micro$Source <- "MCSE Microplot (Table 154)"

# this is already summed across i think subplots (the 0.1 m2 samples) cuz total area sampled is 1
# anpp per treatment replicate disturbance regime....

names(micro)
anpp_micro <-  micro %>% 
  group_by(Year, Treatment, Replicate, Disturbed_Microplot, Fertilized_Microplot,Source) %>% 
  summarise (anpp = sum (biomass_g_m2))

head(anpp_micro)




####################################################################
################# combine these datasets into ANPP df and SPCOMP df
####################################################################


#FIRST:
# combine these ANPP datasets
names(anpp_t7)
names(anpp_micro)


# merge main MCSE with microplot
anpp_KBS_T7 <- rbind(anpp_t7, anpp_micro)

head(anpp_KBS_T7 )
unique(anpp_KBS_T7$Treatment) # dope.

hist(anpp_KBS_T7$anpp)


# NEXT : GET Percent cover dataset (sp comp)

head(t7)
head(micro)


# merge main MCSE with microplot
t7_with_ANPP <- merge(t7, anpp_t7, by = c("Year", "Treatment", "Station",
                                          "Replicate", "Source"))
head(t7_with_ANPP)

# get psesudo percent cover by dividing plant by total for ANPP...
t7_with_ANPP$Pseudo_PercCover <- t7_with_ANPP$biomass_g_m2 / t7_with_ANPP$anpp * 100
head(t7_with_ANPP)

# micro
micro_with_ANPP <- merge(micro, anpp_micro, by = c("Year", "Treatment", "Disturbed_Microplot", "Fertilized_Microplot",
                                                   "Replicate", "Source"))
head(micro_with_ANPP)

micro_with_ANPP$Pseudo_PercCover <- micro_with_ANPP$biomass_g_m2 / micro_with_ANPP$anpp * 100

head(micro_with_ANPP)


# bind together
head(t7_with_ANPP)
head(micro_with_ANPP)
allt7_SpComp <- dplyr:: bind_rows (t7_with_ANPP, micro_with_ANPP)

allt7_SpComp




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

# write a new .csv with the cleaned and merged data and upload to the shared google drive L1 folder
write.csv(allt7_SpComp_tp, file.path(L1_dir, "./KBS_MCSE_T7_SpComp.csv"), row.names=F)
#write.csv(allt7_SpComp_tp, "KBS_MCSE_T7_SpComp.csv")






allt7_ANPP_tp <- merge(anpp_KBS_T7, weatheryear, by = "Year")

allt7_ANPP_tp

# write a new .csv with the cleaned and merged data and upload to the shared google drive L1 folder
write.csv(allt7_ANPP_tp, file.path(L1_dir, "./KBS_MCSE_T7_ANPP.csv"), row.names=F)
#write.csv(allt7_ANPP_tp, "KBS_MCSE_T7_ANPP.csv")


ggplot(allt7_ANPP_tp, aes (x = annualprecip, y = ANPP))
