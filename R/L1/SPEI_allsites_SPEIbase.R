# TITLE:        LTER Grassland Rock: SPEI from SPEIbase for all sites
# AUTHORS:      Caitlin Broderick
# COLLABORATORS:  
# DATA INPUT:   Data imported as csv files from shared Google drive L0 folder
# DATA OUTPUT:  (not yet, but generate a list of years /spei/spei categories to match up with plant data)
# PROJECT:      LTER Grassland Rock
# DATE:         April 2024



# SPEIbase: https://spei.csic.es/database.html  <- same as Isbell paper. 

#Previous results suggest that primary productivity
#responds to approximately annual water balances in temperate grasslands

# and isbell paper used different spei values - i am using SPEI 12 as the default
# this was the first they considered - annual water balances
#note i left the lat long in the file names for easy retrievel of different SPEI metric


# Clear all existing data
rm(list=ls())

# Load packaages
library(tidyverse)
library(lubridate)

# Set working directory 
L0_dir <- Sys.getenv("L0DIR")
L1_dir <- Sys.getenv("L1DIR")
list.files(L0_dir)

kbs <- read.csv(file.path(L0_dir, "KBS/KBS_SPEI_-85.250000_42.250000.csv"), stringsAsFactors = FALSE)
knz <- read.csv(file.path(L0_dir, "KNZ/KNZ_SPEI_-96.750000_39.250000.csv"), stringsAsFactors = FALSE)
cdr <- read.csv(file.path(L0_dir, "CDR/CDR_SPEI_-93.250000_45.250000.csv"), stringsAsFactors = FALSE)


# lubridate the date columns
kbs$dates <- lubridate::mdy(as.character(kbs$dates))
knz$dates <- lubridate::mdy(as.character(knz$dates))
cdr$dates <- lubridate::mdy(as.character(cdr$dates))

# isbell first filtered for peak biomass harvest at each site. 
kbs.aug <- kbs %>% filter(month(dates) == 8 & !is.nan(spei12))
knz.aug <- knz %>% filter(month(dates) == 8 & !is.nan(spei12))
cdr.aug <- kbs %>% filter(month(dates) == 8 & !is.nan(spei12))


# mean? just curious. should be about zero
mean(kbs.aug$spei12) 
mean(knz.aug$spei12) 
mean(cdr.aug$spei12)  # looking good. 





# calculate interquartile range - "normal range"
# this is the diffeerence between upper and lower quartile ranges. 
# Calculate 1 in 10 year event 

quantile(kbs.aug$spei12, c(0.1, 0.25, 0.75, 0.9))
quantile(knz.aug$spei12, c(0.1, 0.25, 0.75, 0.9))
quantile(cdr.aug$spei12, c(0.1, 0.25, 0.75, 0.9))


# add column for classifying years based on quartiles !!!!!
kbs.aug <-kbs.aug   %>% mutate(spei_category = case_when(
  kbs.aug$spei12 < -1.4 ~ "Extreme dry", 
  kbs.aug$spei12 >= -1.4 & kbs.aug$spei12 < -0.7 ~ "Moderate dry", 
  kbs.aug$spei12 >= -0.7 & kbs.aug$spei12 <= 0.8 ~ "Normal", 
  kbs.aug$spei12 > 0.8 & kbs.aug$spei12 <= 1.3 ~ "Moderate wet", 
  kbs.aug$spei12 > 1.3 ~ "Extreme wet"
))


# add column for classifying years based on quartiles !!!!!
knz.aug <-knz.aug   %>% mutate(spei_category = case_when(
  knz.aug$spei12 < -1.4 ~ "Extreme dry", 
  knz.aug$spei12 >= -1.4 & knz.aug$spei12 < -0.7 ~ "Moderate dry", 
  knz.aug$spei12 >= -0.7 & knz.aug$spei12 <= 0.7 ~ "Normal", 
  knz.aug$spei12 > 0.7 & knz.aug$spei12 <= 1.2 ~ "Moderate wet", 
  knz.aug$spei12 > 1.2 ~ "Extreme wet"
))


# add column for classifying years based on quartiles !!!!!
cdr.aug <-cdr.aug   %>% mutate(spei_category = case_when(
  cdr.aug$spei12 < -1.4 ~ "Extreme dry", 
  cdr.aug$spei12 >= -1.4 & cdr.aug$spei12 < -0.7 ~ "Moderate dry", 
  cdr.aug$spei12 >= -0.7 & cdr.aug$spei12 <= 0.8 ~ "Normal", 
  cdr.aug$spei12 > 0.8 & cdr.aug$spei12 <= 1.3 ~ "Moderate wet", 
  cdr.aug$spei12 > 1.3 ~ "Extreme wet"
))


# add year column for graphing purposes
kbs.aug$year <- year(kbs.aug$dates)
knz.aug$year <- year(knz.aug$dates)
cdr.aug$year <- year(cdr.aug$dates)

# order factor

kbs.aug$spei_category <- factor(kbs.aug$spei_category, levels = c("Extreme dry", "Moderate dry", "Normal", "Moderate wet", "Extreme wet"))
knz.aug$spei_category <- factor(knz.aug$spei_category, levels = c("Extreme dry", "Moderate dry", "Normal", "Moderate wet", "Extreme wet"))
cdr.aug$spei_category <- factor(cdr.aug$spei_category, levels = c("Extreme dry", "Moderate dry", "Normal", "Moderate wet", "Extreme wet"))

ggplot(kbs.aug , aes (x = year, y = spei12)) +
  geom_line(alpha = 0.5) + 
  geom_point( aes(color = spei_category), size = 3) + theme_classic() +
  scale_color_manual(values = c("#F5191CFF", "#E78200FF", "#EAC728FF", "#81BB95FF", "#3B99B1FF"))



ggplot(knz.aug , aes (x = year, y = spei12)) +
  geom_line(alpha = 0.5) + 
  geom_point( aes(color = spei_category), size = 3) + theme_classic() +
  scale_color_manual(values = c("#F5191CFF", "#E78200FF", "#EAC728FF", "#81BB95FF", "#3B99B1FF"))



ggplot(cdr.aug , aes (x = year, y = spei12)) +
  geom_line(alpha = 0.5) + 
  geom_point( aes(color = spei_category), size = 3) + theme_classic() +
  scale_color_manual(values = c("#F5191CFF", "#E78200FF", "#EAC728FF", "#81BB95FF", "#3B99B1FF"))
