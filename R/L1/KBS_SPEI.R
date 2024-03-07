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


# add month and doy column. 
met$date <- lubridate::mdy(met$date)
met$month <- month(met$date)
met$doy <- yday(met$date)

head(met)

met_month <- met %>% 
  group_by(Year, month) %>% 
  summarise(air_temp_min = mean(air_temp_min, na.rm = TRUE),  #AVERAGE of all the temp stuff?
            air_temp_max = mean(air_temp_max, na.rm = TRUE), #AVERAGE of all the temp stuff?
            air_temp_mean = mean(air_temp_mean, na.rm = TRUE), #AVERAGE of all the temp stuff?
            precip_sum = sum(precipitation, na.rm = TRUE)) # but use the SUM of the months precip...

met %>% 
  filter(doy > 84 & doy < 249) %>% 
  group_by(month ) %>% 
  count() # based on robinson paper cutoffs far fewer obs for march and sept

#two different ways to get gs precip - toggle back and forth
    # doy cutoff vs month cutoff
met_month_gs <- met %>% 
  #filter(doy > 84 & doy < 249)  %>% # robinson paper: GS 85–248
  filter(month > 3 & month < 9) %>% 
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
# same for gs
met_month_gs$PET <- SPEI::hargreaves(Tmin= met_month_gs$air_temp_min ,Tmax = met_month_gs$air_temp_max, Pre= met_month_gs$precip_sum, lat=42.4, na.rm = TRUE, verbose = TRUE)
# method 2: thornthwaite. this seemed less complicated but perhaps less precise
#met_month$PET <- SPEI::thornthwaite(met_month$air_temp_mean , lat=42.4, na.rm = TRUE)


# calculate climatic water balance
met_month$BAL <- met_month$precip_sum - met_month$PET
met_month_gs$BAL <- met_month_gs$precip_sum - met_month_gs$PET
# calculate SPEI by month. (SPEI-1) 
spei <- spei(met_month$BAL, 1, na.rm = T)
spei$fitted
met_month$spei <- as.vector (spei$fitted )

spei_gs <- spei(met_month_gs$BAL, 1, na.rm = T) # This may be wrong. spei assumes we start in January and it's a 12 month time series
spei_gs$fitted
met_month_gs$spei <- as.vector (spei_gs$fitted )

# Calculate SPEI-12
spei_12 <- spei(met_month$BAL, 12, na.rm = T)
spei_12$fitted
met_month$spei_12 <- as.vector(spei_12$fitted )


# whole year: 
ggplot(met_month, aes (x = month, y = spei)) +
  geom_hline(yintercept=0, color = "blue") + 
  geom_hline(yintercept=-1, color = "orange") +  # slette paper - moderately dry
  geom_hline(yintercept=-1.5, color = "red") +  # slette paper - severely dry
  geom_hline(yintercept=-2, color = "darkred") +  # slette paper - extremely dry
  geom_line(linewidth = 1.2) +
  facet_wrap(vars(Year)) # notice the very dry 2012
ggplot(met_month, aes (x = month, y = spei_12)) +
  geom_hline(yintercept=0, color = "blue") + 
  geom_hline(yintercept=-1, color = "orange") +  # slette paper - moderately dry
  geom_hline(yintercept=-1.5, color = "red") +  # slette paper - severely dry
  geom_hline(yintercept=-2, color = "darkred") +  # slette paper - extremely dry
  geom_line(linewidth = 1.2) +
  facet_wrap(vars(Year)) # notice the very dry 2012

#growing season
ggplot(met_month_gs, aes (x = month, y = spei)) +
  geom_hline(yintercept=0, color = "blue") + 
  geom_hline(yintercept=-1, color = "orange") +  # slette paper - moderately dry
  geom_hline(yintercept=-1.5, color = "red") +  # slette paper - severely dry
  geom_hline(yintercept=-2, color = "darkred") +  # slette paper - extremely dry
  geom_line(linewidth = 1.2) +
  facet_wrap(vars(Year)) # notice the very dry 2012

# filter out SPEI for August (harvest date)
met_month_aug <- met_month %>%
  filter(month == 8)

# calculate SPEI (and total ppt) for each year 
# not sure if you can take a mean like this
met_year <- met_month %>%  
  group_by(Year) %>% 
  summarise(meanspei = mean(spei),
            totalrain = sum(precip_sum)) 
print(met_year, n=40) # cool that this already highlights the ones below 0....


met_year_gs <- met_month_gs %>%  
  group_by(Year) %>% 
  summarise(meanspei = mean(spei),
            totalrain = sum(precip_sum)) 
print(met_year_gs, n=40) # cool that this already highlights the ones below 0....

# how close does total RAIN track my calculated SPEI????
ggplot(met_year, aes (x = totalrain, y = meanspei)) +
  geom_point() # not too bad!


ggplot(met_year_gs, aes (x = totalrain, y = meanspei)) +
  geom_point() # only growing season


# which years are super dry?
ggplot(met_year, aes (x = totalrain, y = meanspei, label = Year)) +
  geom_text() +
  theme_classic() # cool. 
# 1999, 2002, 2005, 2012 are lowest by spei
# 1996 and 1998 also dry esp by total precip


ggplot(met_year_gs, aes (x = totalrain, y = meanspei, label = Year)) +
  geom_text() +
  theme_classic() # 


ggplot(met_year, aes (x = meanspei)) +
  geom_density()

  ggplot(met_year_gs, aes (x = meanspei)) +
           geom_density()
  
met_year <- met_year %>%
  mutate(drought_cat = case_when(
    meanspei <= -1 ~ "moderately_dry",
    meanspei > -1 & meanspei <= -0.5 ~ "slightly_dry",
    meanspei > -0.5 & meanspei < 0.5 ~ "normal",
    meanspei > 0.5 ~ "slightly_wet"
  ))

met_year_gs <- met_year_gs %>%
  mutate(drought_cat = case_when(
    meanspei <= -1 ~ "moderately_dry",
    meanspei > -1 & meanspei <= -0.5 ~ "slightly_dry",
    meanspei > -0.5 & meanspei < 0.5 ~ "normal",
    meanspei > 0.5 ~ "slightly_wet"
  ))
  
# Mean yearly SPEI over time
ggplot(met_year, aes (x = Year, y = meanspei)) +
  geom_hline(yintercept=1, color = "blue") + 
  geom_hline(yintercept=0.5, color = "lightblue") + 
  geom_hline(yintercept=0, color = "black") + 
  geom_hline(yintercept=-0.5, color = "gold") +  
  geom_hline(yintercept=-1, color = "orange") +  # slette paper - moderately dry
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5, aes(color = drought_cat)) +
  scale_color_manual(values = c("slightly_dry" = "gold", "normal" = "black", "slightly_wet" = "lightblue")) +
  theme_bw() +
  theme(legend.position="none")
             

# Mean growing season SPEI over time
ggplot(met_year_gs, aes (x = Year, y = meanspei)) +
  geom_hline(yintercept=1, color = "blue") + 
  geom_hline(yintercept=0.5, color = "lightblue") + 
  geom_hline(yintercept=0, color = "black") + 
  geom_hline(yintercept=-0.5, color = "gold") +  
  geom_hline(yintercept=-1, color = "orange") +  # slette paper - moderately dry
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5, aes(color = drought_cat)) +
  scale_color_manual(values = c("slightly_dry" = "gold", "normal" = "black", "slightly_wet" = "lightblue")) +
  theme_bw() +
  theme(legend.position="none")

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
met_month_aug <- met_month_aug %>%
  mutate(drought_cat_1 = case_when(
    spei <= -1.5 ~ "severely_dry",
    spei > -1.5 & spei <= -1 ~ "moderately_dry",
    spei > -1 & spei <= -0.5 ~ "slightly_dry",
    spei > -0.5 & spei < 0.5 ~ "normal",
    spei >= 0.5 & spei < 1 ~ "slightly_wet",
    spei >= 1 & spei < 1.5 ~ "moderately_wet",
    spei >= 1.5 ~ "severely_wet"
  ))

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
  labs(y = "SPEI-12", title = "KBS SPEI")
# dev.off()

# SPEI-1 for August (harvest month)
# I don't think SPEI-1 makes much sense. SPEI-3 seems like a minimum
ggplot(met_month_aug, aes (x = Year, y = spei)) +
  geom_hline(yintercept=1.5, color = "darkblue") + 
  geom_hline(yintercept=1, color = "blue") + 
  geom_hline(yintercept=0.5, color = "lightblue") + 
  geom_hline(yintercept=0, color = "black") + 
  geom_hline(yintercept=-0.5, color = "gold") +  
  geom_hline(yintercept=-1, color = "orange") +  # slette paper - moderately dry
  geom_hline(yintercept=-1.5, color = "red") +  # slette paper - severely dry
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5, aes(color = drought_cat_1)) +
  scale_color_manual(values = c("severely_dry" = "red", "moderately_dry" = "orange", "slightly_dry" = "gold", "normal" = "black", "slightly_wet" = "lightblue", "moderately_wet" = "blue", "severely_wet" = "darkblue")) +
  theme_bw() +
  theme(legend.position="none")
