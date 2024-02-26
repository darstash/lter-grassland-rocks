# TITLE:        LTER Grassland Rock: SPEI work
# AUTHORS:      Caitlin Broderick
# COLLABORATORS:  
# DATA INPUT:   Data imported as csv files from shared Google drive L0 folder
# DATA OUTPUT:    
# PROJECT:      LTER Grassland Rock
# DATE:         Feb 2024
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
met <- read.csv(file.path(L0_dir, "KBS/7-lter+weather+station+daily+precip+and+air+temp+1698756534_L0.csv"), stringsAsFactors = FALSE)
names(met)
?spei 
# input variable: series of climatic water balance  (ppt - PET)
# so have to calculate PET first
# also says "The SPEI and the SPI were defined for monthly data"
  # so gonna average/sum meteorological vars by month first
  # just doing na.rm for now, we can talk about this later!!!



met$date <- lubridate::mdy(met$date)
met$month <- month(met$date)

met_month <- met %>% 
  select(-date) %>% 
  group_by(Year, month) %>% 
  summarise(air_temp_min = mean(air_temp_min, na.rm = TRUE),  #AVERAGE of all the temp stuff?
            air_temp_max = mean(air_temp_max, na.rm = TRUE), #AVERAGE of all the temp stuff?
            air_temp_mean = mean(air_temp_mean, na.rm = TRUE), #AVERAGE of all the temp stuff?
            precip_sum = sum(precipitation, na.rm = TRUE)) # but use the SUM of the months precip...



# calculate. PET
# method 1: hargreaves
met_month$PET <- SPEI::hargreaves(Tmin= met_month$air_temp_min ,Tmax = met_month$air_temp_max, Pre= met_month$precip_sum, lat=42.4, na.rm = TRUE, verbose = TRUE)
 # Tmin: a numeric vector, tsvector, matrix, tsmatrix, or 3-d array of monthly mean daily minimum temperatures, ºC.
 # Tmax: a numeric vector, tsvector, matrix, tsmatrix, or 3-d array of monthly mean daily maximum temperatures, ºC.
 # Pre: optional, a numeric vector, tsvector, matrix, tsmatrix, or 3-d array of monthly total precipitation, mm.
 # looked up the latitude for kbs...

# method 2: thornthwaite. this seemed less complicated but perhaps less precise
#met_month$PET <- SPEI::thornthwaite(met_month$air_temp_mean , lat=42.4, na.rm = TRUE)


# calculate climatic water balance
met_month$BAL <- met_month$precip_sum - met_month$PET

# calculate SPEI by month. 
spei <- spei(met_month$BAL, 1, na.rm = T)
spei$fitted

met_month$spei <- as.vector (spei$fitted )


ggplot(met_month, aes (x = month, y = spei)) +
  geom_hline(yintercept=0, color = "red") + 
  geom_line(linewidth = 1.2) +
  facet_wrap(vars(Year)) # notice the very dry 2012


# calculate SPEI (and total ppt) for each year
met_year <- met_month %>%  
  group_by(Year) %>% 
  summarise(meanspei = mean(spei),
            totalrain = sum(precip_sum)) 
print(met_year, n=40) # cool that this already highlights the ones below 0....


# how close does total RAIN track my calculated SPEI????
ggplot(met_year, aes (x = totalrain, y = meanspei)) +
  geom_point() # not too bad!

# which years are super dry?
ggplot(met_year, aes (x = totalrain, y = meanspei, label = Year)) +
  geom_text() +
  theme_classic() # cool. 
# 1999, 2002, 2005, 2012 are lowest by spei
# 1996 and 1998 also dry esp by total precip



