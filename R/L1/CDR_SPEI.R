# TITLE:        LTER Grassland Rock: CDR SPEI work
# AUTHORS:      Max Zaret
# COLLABORATORS:  
# DATA INPUT:   Data imported as csv files from shared Google drive L0 folder
# DATA OUTPUT:    
# PROJECT:      LTER Grassland Rock
# DATE:         March 2024
# got SPEI R package based on this website:
#https://climatedataguide.ucar.edu/climate-data/standardized-precipitation-evapotranspiration-index-spei

# Clear all existing data
rm(list=ls())

# Load packages
library(tidyverse)
library(gtools)
library(SPEI)
library(lubridate)

# Set working directory 
L0_dir <- Sys.getenv("L0DIR")
L1_dir <- Sys.getenv("L1DIR")
list.files(L0_dir)



# Read in CSV files
met <- read.csv(paste(L0_dir, "CDR/e080_Daily climate summary.csv", sep="/"))
names(met)
?spei 
# input variable: series of climatic water balance  (ppt - PET)
# so have to calculate PET first
# also says "The SPEI and the SPI were defined for monthly data"
# so gonna average/sum meteorological vars by month first
# just doing na.rm for now, we can talk about this later!!!


# add month and doy column. 
met$date <- lubridate::mdy(met$Date)
met$Year <- year(met$date)
met$month <- month(met$date)
met$doy <- yday(met$date)

head(met)

met_month <- met %>% 
  group_by(Year, month) %>% 
  summarise(air_temp_min = mean(((MinTemp.degF.-32)/1.8), na.rm = TRUE),  #convert to C
            air_temp_max = mean(((MaxTemp.degF.-32)/1.8), na.rm = TRUE), #Convert to C
            #air_temp_mean = mean(mean(MaxTemp.degF.,MinTemp.degF.), na.rm = TRUE), #AVERAGE of all the temp stuff?
            precip_sum = sum((Precip.inches.*25.4), na.rm = TRUE)) # convert to mm

# calculate. PET
# method 1: hargreaves
met_month$PET <- SPEI::hargreaves(Tmin= met_month$air_temp_min ,Tmax = met_month$air_temp_max, Pre= met_month$precip_sum, lat=45.24, na.rm = TRUE, verbose = TRUE)
# Tmin: a numeric vector, tsvector, matrix, tsmatrix, or 3-d array of monthly mean daily minimum temperatures, ºC.
# Tmax: a numeric vector, tsvector, matrix, tsmatrix, or 3-d array of monthly mean daily maximum temperatures, ºC.
# Pre: optional, a numeric vector, tsvector, matrix, tsmatrix, or 3-d array of monthly total precipitation, mm.
# looked up the latitude for CDR


# calculate climatic water balance
met_month$BAL <- met_month$precip_sum - met_month$PET
# calculate SPEI by month. (SPEI-1) 
spei <- spei(met_month$BAL, 1, na.rm = T)
spei$fitted
met_month$spei <- as.vector (spei$fitted )

# Calculate SPEI-12
spei_12 <- spei(met_month$BAL, 12, na.rm = T)
spei_12$fitted
met_month$spei_12 <- as.vector(spei_12$fitted )


# whole year: 
ggplot(met_month, aes (x = month, y = spei_12)) +
  geom_hline(yintercept=0, color = "blue") + 
  geom_hline(yintercept=-1, color = "orange") +  # slette paper - moderately dry
  geom_hline(yintercept=-1.5, color = "red") +  # slette paper - severely dry
  geom_hline(yintercept=-2, color = "darkred") +  # slette paper - extremely dry
  geom_line(linewidth = 1.2) +
  facet_wrap(vars(Year)) # notice the very dry 2012

met_year <- met_month %>%
  group_by(Year) %>%
  summarize(spei_12 = mean(spei_12)) %>%
  mutate(drought_cat_12 = case_when(
    spei_12 <= -1.5 ~ "severely_dry",
    spei_12 > -1.5 & spei_12 <= -1 ~ "moderately_dry",
    spei_12 > -1 & spei_12 <= -0.5 ~ "slightly_dry",
    spei_12 > -0.5 & spei_12 < 0.5 ~ "normal",
    spei_12 >= 0.5 & spei_12 < 1 ~ "slightly_wet",
    spei_12 >= 1 & spei_12 < 1.5 ~ "moderately_wet",
    spei_12 >= 1.5 ~ "severely_wet"
  ))

ggplot(met_year, aes (x = Year, y = spei_12)) +
  geom_hline(yintercept=1.5, color = "darkblue") + 
  geom_hline(yintercept=1, color = "blue") + 
  geom_hline(yintercept=0.5, color = "lightblue") + 
  geom_hline(yintercept=0, color = "black") + 
  geom_hline(yintercept=-0.5, color = "gold") +  
  geom_hline(yintercept=-1, color = "orange") +  # slette paper - moderately dry
  geom_hline(yintercept=-1.5, color = "red") +  # slette paper - severely dry
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5, aes(color = drought_cat_12)) +
  scale_color_manual(values = c("severely_dry" = "red", "moderately_dry" = "orange", "slightly_dry" = "gold", "normal" = "black", "slightly_wet" = "lightblue", "moderately_wet" = "blue", "severely_wet" = "darkblue")) +
  theme_bw() +
  theme(legend.position="none") +
  labs(y = "SPEI-12", title = "CDR SPEI")


# filter out SPEI for August (harvest date)
met_month_aug <- met_month %>%
  filter(month == 8)

met_month_aug <- met_month_aug %>%
  mutate(drought_cat_12 = case_when(
    spei_12 <= -1.5 ~ "severely_dry",
    spei_12 > -1.5 & spei_12 <= -1 ~ "moderately_dry",
    spei_12 > -1 & spei_12 <= -0.5 ~ "slightly_dry",
    spei_12 > -0.5 & spei_12 < 0.5 ~ "normal",
    spei_12 >= 0.5 & spei_12 < 1 ~ "slightly_wet",
    spei_12 >= 1 & spei_12 < 1.5 ~ "moderately_wet",
    spei_12 >= 1.5 ~ "severely_wet"
  ))

#for filtering by normal years#
normal_years <- met_year %>%
  filter(drought_cat_12 == "normal") %>%
  select(Year)

# SPEI-12 for August (harvest month)
# pdf("spei12.pdf", width = 5, height = 5)
ggplot(met_month_aug, aes (x = Year, y = spei_12)) +
  geom_hline(yintercept=1.5, color = "darkblue") + 
  geom_hline(yintercept=1, color = "blue") + 
  geom_hline(yintercept=0.5, color = "lightblue") + 
  geom_hline(yintercept=0, color = "black") + 
  geom_hline(yintercept=-0.5, color = "gold") +  
  geom_hline(yintercept=-1, color = "orange") +  # slette paper - moderately dry
  geom_hline(yintercept=-1.5, color = "red") +  # slette paper - severely dry
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5, aes(color = drought_cat_12)) +
  scale_color_manual(values = c("severely_dry" = "red", "moderately_dry" = "orange", "slightly_dry" = "gold", "normal" = "black", "slightly_wet" = "lightblue", "moderately_wet" = "blue", "severely_wet" = "darkblue")) +
  theme_bw() +
  theme(legend.position="none") +
  labs(y = "SPEI-12", title = "CDR SPEI")
# dev.off()
