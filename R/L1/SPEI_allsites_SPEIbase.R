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

kbs <- read.csv(file.path(L0_dir, "KBS/-85.250000_42.250000.csv"), stringsAsFactors = FALSE)
knz <- read.csv(file.path(L0_dir, "KNZ/-96.750000_39.250000.csv"), stringsAsFactors = FALSE)
cdr <- read.csv(file.path(L0_dir, "CDR/-93.250000_45.250000.csv"), stringsAsFactors = FALSE)


# separate columns
# have to do it here! no maatter what I do, opening this in excel goofs it all up
cdr <- cdr %>%  separate_wider_delim(dates.spei12, delim = ";", names = c("date", "spei12"))
kbs <- kbs %>%  separate_wider_delim(dates.spei12, delim = ";", names = c("date", "spei12"))
knz <- knz %>%  separate_wider_delim(dates.spei12, delim = ";", names = c("date", "spei12"))

# lubridate the date columns
kbs$date <- lubridate::ymd(as.character(kbs$date))
knz$date <- lubridate::ymd(as.character(knz$date))
cdr$date <- lubridate::ymd(as.character(cdr$date))


kbs$spei12 <- as.numeric (as.character(kbs$spei12))
knz$spei12 <- as.numeric (as.character(knz$spei12))
cdr$spei12 <- as.numeric (as.character(cdr$spei12))


# add year column for graphing purposes
kbs$year <- year(kbs$date)
knz$year <- year(knz$date)
cdr$year <- year(cdr$date)

max(kbs$year)
max(knz$year)
max(cdr$year)
# isbell first filtered for peak biomass harvest at each site. 
  # also, only want 100 years? can add that later 
kbs.aug <- kbs %>% filter(month(date) == 8 & !is.nan(spei12))
knz.aug <- knz %>% filter(month(date) == 8 & !is.nan(spei12))
cdr.aug <- kbs %>% filter(month(date) == 8 & !is.nan(spei12))


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

# order factor

kbs.aug$spei_category <- factor(kbs.aug$spei_category, levels = c("Extreme dry", "Moderate dry", "Normal", "Moderate wet", "Extreme wet"))
knz.aug$spei_category <- factor(knz.aug$spei_category, levels = c("Extreme dry", "Moderate dry", "Normal", "Moderate wet", "Extreme wet"))
cdr.aug$spei_category <- factor(cdr.aug$spei_category, levels = c("Extreme dry", "Moderate dry", "Normal", "Moderate wet", "Extreme wet"))

kbs_spei <- ggplot(kbs.aug , aes (x = year, y = spei12)) +
  geom_line(alpha = 0.5) + 
  geom_point( aes(color = spei_category), size = 1.25) + theme_bw() +
  scale_x_continuous(breaks = seq(1900, 2025, by = 10)) +
  annotate( "text", label = "KBS",
    x = 1910, y = 2, size = 5, colour = "black") + 
  scale_color_manual(values = c("#F5191CFF", "#E78200FF", "#EAC728FF", "#81BB95FF", "#3B99B1FF"))
kbs_spei
ggsave("KBS_SPEI.png",plot =kbs_spei, dpi = 300, width =8, height = 4, units = "in")


knz_spei <- ggplot(knz.aug , aes (x = year, y = spei12)) +
  geom_line(alpha = 0.5) + 
  geom_point( aes(color = spei_category), size = 1.25) + theme_bw() +
  scale_x_continuous(breaks = seq(1900, 2025, by = 10)) +
  annotate( "text", label = "KNZ",
            x = 1910, y = 2, size = 5, colour = "black") + 
  scale_color_manual(values = c("#F5191CFF", "#E78200FF", "#EAC728FF", "#81BB95FF", "#3B99B1FF"))
knz_spei
ggsave("KNZ_SPEI.png",plot =knz_spei, dpi = 300, width =8, height = 4, units = "in")


cdr_spei <- ggplot(cdr.aug , aes (x = year, y = spei12)) +
  geom_line(alpha = 0.5) + 
  geom_point( aes(color = spei_category), size = 1.25) + theme_bw() +
  scale_x_continuous(breaks = seq(1900, 2025, by = 10)) +
  annotate( "text", label = "CDR",
            x = 1910, y = 2, size = 5, colour = "black") + 
  scale_color_manual(values = c("#F5191CFF", "#E78200FF", "#EAC728FF", "#81BB95FF", "#3B99B1FF"))
cdr_spei
ggsave("CDR_SPEI.png",plot =cdr_spei, dpi = 300, width =8, height = 4, units = "in")

