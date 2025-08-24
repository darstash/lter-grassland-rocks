# LTER Grassland Rock: Calculate resistance and resilience
# AUTHORS: Rosalie Terry
# COLLABORATORS: Matt Nieland, Joshua Ajowele  
# DATA INPUT: .csv file of plot metrics imported from Google Drive L2 folder
# DATA OUTPUT: .csv of updated plot metrics dataset that includes resistance and resilience  
# PROJECT: LTER Grassland Rock
# DATE: October 2024

# Clear all existing data
rm(list=ls())

# Load packages
library(tidyverse)

# Set working directory 
L0_dir <- Sys.getenv("L0DIR")
L1_dir <- Sys.getenv("L1DIR")
L2_dir <- Sys.getenv("L2DIR")
list.files(L2_dir)

# Read in CSV files
plot_metrics <- read.csv(file.path(L2_dir, "plot_metrics_SPEI_diversity_L2.csv"), stringsAsFactors = FALSE)

# SPEI6 ----
# Create dataframe that contains extreme years at each site
extreme_per_site <- plot_metrics %>% 
  filter(spei6_category == "Extreme wet" | spei6_category == "Extreme dry")
  

# Create dataframe that contains average "normal" biomass for each plot
biomass_normal <- plot_metrics %>% 
  filter(spei6_category == "Normal" | spei6_category == "Moderate wet" | spei6_category == "Moderate dry") %>% # not enough just normal years
  group_by(uniqueid) %>% # group by plot
  summarize(average_normal_biomass = mean(plot_biomass, na.rm=T)) #na.rm=T is required to avoid NAs when the plots are missing values for certain years

# Create dataframe of metrics for plots during extreme climate events that includes resistance and resilience
df <- data.frame(ex_year = integer(), site = factor(), uniqueid = character(), resistance = numeric(), resilience = numeric())

for(i in 1:nrow(extreme_per_site)) {
  subs1 <- plot_metrics %>% 
    filter(year %in% extreme_per_site[i, "year"] &
             site %in% extreme_per_site[i,"site"]) %>% 
    rename("biomass_during" = "plot_biomass") %>% 
    rename("ex_year" = "year")
  subs2 <- plot_metrics %>% 
    filter(year %in% (extreme_per_site[i, "year"] + 1) &
             site %in% extreme_per_site[i,"site"]) %>% 
    rename("biomass_after" = "plot_biomass")
  subs <- merge (subs1, subset(subs2, select = -c(year)), by = c("site", "uniqueid"))
  subs <- merge(subs, biomass_normal, by = c("uniqueid"), all.x = T)
  subs$resistance <- subs$average_normal_biomass / abs((subs$biomass_during - subs$average_normal_biomass)) # resistance following Isbell et al. 2015
  subs$resilience <- abs((subs$biomass_during - subs$average_normal_biomass) / (subs$biomass_after - subs$average_normal_biomass)) # resilience following Isbell et al. 2015
  common_names <- intersect(colnames(df), colnames(subs))
  df <- rbind.data.frame(df[, common_names], subs[, common_names])
}

# save dataframe of plot resistance and resilience to ECEs as .csv in Google Drive L2 folder
write.csv(df, file.path(L2_dir, "./ece_resist_resil_L2.csv"), row.names=F)


# SPEI9 ----
# Create dataframe that contains extreme years at each site
extreme_per_site9 <- plot_metrics %>% 
  filter(spei9_category == "Extreme wet" | spei9_category == "Extreme dry")


# Create dataframe that contains average "normal" biomass for each plot
biomass_normal9 <- plot_metrics %>% 
  filter(spei9_category == "Normal" | spei9_category == "Moderate wet" | spei9_category == "Moderate dry") %>% # not enough just normal years
  group_by(uniqueid) %>% # group by plot
  summarize(average_normal_biomass = mean(plot_biomass, na.rm=T)) #na.rm=T is required to avoid NAs when the plots are missing values for certain years

# Create dataframe of metrics for plots during extreme climate events that includes resistance and resilience
df9 <- data.frame(ex_year = integer(), site = factor(), uniqueid = character(), resistance = numeric(), resilience = numeric())

for(i in 1:nrow(extreme_per_site9)) {
  subs1 <- plot_metrics %>% 
    filter(year %in% extreme_per_site9[i, "year"] &
             site %in% extreme_per_site9[i,"site"]) %>% 
    rename("biomass_during" = "plot_biomass") %>% 
    rename("ex_year" = "year")
  subs2 <- plot_metrics %>% 
    filter(year %in% (extreme_per_site9[i, "year"] + 1) &
             site %in% extreme_per_site9[i,"site"]) %>% 
    rename("biomass_after" = "plot_biomass")
  subs <- merge (subs1, subset(subs2, select = -c(year)), by = c("site", "uniqueid"))
  subs <- merge(subs, biomass_normal9, by = c("uniqueid"), all.x = T)
  subs$resistance <- subs$average_normal_biomass / abs((subs$biomass_during - subs$average_normal_biomass)) # resistance following Isbell et al. 2015
  subs$resilience <- abs((subs$biomass_during - subs$average_normal_biomass) / (subs$biomass_after - subs$average_normal_biomass)) # resilience following Isbell et al. 2015
  common_names <- intersect(colnames(df9), colnames(subs))
  df9 <- rbind.data.frame(df9[, common_names], subs[, common_names])
}

# save dataframe of plot resistance and resilience to ECEs as .csv in Google Drive L2 folder
write.csv(df9, file.path(L2_dir, "./ece_resist_resil_spei9_L2.csv"), row.names=F)


## ONLY normal years ----
# Create dataframe that contains extreme years at each site
extreme_per_site9 <- plot_metrics %>% 
  filter(spei9_category == "Extreme wet" | spei9_category == "Extreme dry")


# Create dataframe that contains average "normal" biomass for each plot
biomass_normal9_norm <- plot_metrics %>% 
  filter(spei9_category == "Normal") %>%
  group_by(uniqueid) %>% # group by plot
  summarize(average_normal_biomass = mean(plot_biomass, na.rm=T)) #na.rm=T is required to avoid NAs when the plots are missing values for certain years

# Create dataframe of metrics for plots during extreme climate events that includes resistance and resilience
df9_norm <- data.frame(ex_year = integer(), site = factor(), uniqueid = character(), resistance = numeric(), resilience = numeric())

for(i in 1:nrow(extreme_per_site9)) {
  subs1 <- plot_metrics %>% 
    filter(year %in% extreme_per_site9[i, "year"] &
             site %in% extreme_per_site9[i,"site"]) %>% 
    rename("biomass_during" = "plot_biomass") %>% 
    rename("ex_year" = "year")
  subs2 <- plot_metrics %>% 
    filter(year %in% (extreme_per_site9[i, "year"] + 1) &
             site %in% extreme_per_site9[i,"site"]) %>% 
    rename("biomass_after" = "plot_biomass")
  subs <- merge (subs1, subset(subs2, select = -c(year)), by = c("site", "uniqueid"))
  subs <- merge(subs, biomass_normal9_norm, by = c("uniqueid"), all.x = T)
  subs$resistance <- subs$average_normal_biomass / abs((subs$biomass_during - subs$average_normal_biomass)) # resistance following Isbell et al. 2015
  subs$resilience <- abs((subs$biomass_during - subs$average_normal_biomass) / (subs$biomass_after - subs$average_normal_biomass)) # resilience following Isbell et al. 2015
  common_names <- intersect(colnames(df9_norm), colnames(subs))
  df9_norm <- rbind.data.frame(df9_norm[, common_names], subs[, common_names])
}

# save dataframe of plot resistance and resilience to ECEs as .csv in Google Drive L2 folder
write.csv(df9_norm, file.path(L2_dir, "./ece_resist_resil_spei9_norm_L2.csv"), row.names=F)


### troubleshooting loop and calculations
##### loop for KBS data only for troubleshooting
#KBS_extreme_per_site <- extreme_per_site %>% 
#  filter(site == "KBS") #134 unique plots, 6 extreme events (1993, 2004, 2009, 2017, 2020)
#
#sum(biomass_normal$uniqueid %in% KBS_extreme_per_site$uniqueid)
#
#
#for(i in 1:nrow(KBS_extreme_per_site)) {
#  subs1 <- plot_metrics %>% 
#    filter(year %in% KBS_extreme_per_site[i, "year"] &
#             site %in% KBS_extreme_per_site[i,"site"]) %>% 
#    rename("biomass_during" = "plot_biomass") %>% 
#    rename("ex_year" = "year")
#  subs2 <- plot_metrics %>% 
#    filter(year %in% (KBS_extreme_per_site[i, "year"] + 1) &
#             site %in% KBS_extreme_per_site[i,"site"]) %>% 
#    rename("biomass_after" = "plot_biomass")
#  subs <- merge (subs1, subset(subs2, select = -c(year)), by = c("site", "uniqueid"))
#  subs <- merge(subs, biomass_normal, by = c("uniqueid"), all.x = T)
#  subs$resistance <- subs$average_normal_biomass / abs((subs$biomass_during - subs$average_normal_biomass))
#  subs$resilience <- abs((subs$biomass_during - subs$average_normal_biomass) / (subs$biomass_after - subs$average_normal_biomass))
#  common_names <- intersect(colnames(df), colnames(subs))
#  df <- rbind.data.frame(df[, common_names], subs[, common_names])
#}
#
#
### same but for CDR
#CDR_extreme_per_site <- extreme_per_site %>% 
#  filter(site == "CDR") #415 unique plots, 12 extreme events (1983, 1984, 1986, 1988, 1993, 2002, 2011, 2012, 2014, 2016, 2019, 2021)
#
#
#for(i in 1:nrow(CDR_extreme_per_site)) {
#  subs1 <- plot_metrics %>% 
#    filter(year %in% CDR_extreme_per_site[i, "year"] &
#             site %in% CDR_extreme_per_site[i,"site"]) %>% 
#    rename("biomass_during" = "plot_biomass") %>% 
#    rename("ex_year" = "year")
#  subs2 <- plot_metrics %>% 
#    filter(year %in% (CDR_extreme_per_site[i, "year"] + 1) &
#             site %in% CDR_extreme_per_site[i,"site"]) %>% 
#    rename("biomass_after" = "plot_biomass")
#  subs <- merge (subs1, subset(subs2, select = -c(year)), by = c("site", "uniqueid"))
#  subs <- merge(subs, biomass_normal, by = c("uniqueid"), all.x = T)
#  subs$resistance <- subs$average_normal_biomass / abs((subs$biomass_during - subs$average_normal_biomass))
#  subs$resilience <- abs((subs$biomass_during - subs$average_normal_biomass) / (subs$biomass_after - subs$average_normal_biomass))
#  common_names <- intersect(colnames(df), colnames(subs))
#  df <- rbind.data.frame(df[, common_names], subs[, common_names])
#}
#
#CDR_ex_plots <- unique(CDR_extreme_per_site$uniqueid)
#CDR_df_plots <- unique(df$uniqueid)
#
#CDR_missing <- !(CDR_ex_plots %in% CDR_df_plots)
#CDR_missing <- CDR_ex_plots[CDR_missing]
#
#CDR_sad <- CDR_extreme_per_site[which(CDR_extreme_per_site$uniqueid %in% CDR_missing),] # 10 plots from e247 in CDR only have 1 year of data - issue with data cleaning
#
##### testing individual loop components
#subs_a <- plot_metrics %>% 
#  filter(year %in% extreme_per_site[1, "year"] &
#           site %in% extreme_per_site[1,"site"]) %>% 
#  rename("biomass_during" = "plot_biomass") %>% 
#  rename("ex_year" = "year")
#
#subs_b <- plot_metrics %>% 
#  filter(year %in% (extreme_per_site[1, "year"] + 1) &
#           site %in% extreme_per_site[1,"site"]) %>% 
#  rename("biomass_after" = "plot_biomass")
#
#subs_c <- merge (subs_a, subset(subs_b, select = -c(year)), by = c("site", "uniqueid"))
#subs_c <- merge(subs_c, biomass_normal, by = c("uniqueid"), all.x = T)
#subs_c$resistance <- subs_c$average_normal_biomass / abs((subs_c$biomass_during - subs_c$average_normal_biomass))
#subs_c$resilience <- abs((subs_c$biomass_during - subs_c$average_normal_biomass) / (subs_c$biomass_after - subs_c$average_normal_biomass))
#
#
#common_names <- intersect(colnames(df), colnames(subs_c))
#df2 <- rbind(df[, common_names], subs_c[, common_names])