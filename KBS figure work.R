# TITLE:         LTER Grassland Rock: KBS - exploring diversity-stability 
# AUTHOR:        Caitlin Broderick
# COLLABORATORS: LTER synthesis group
# DATA INPUT:    Data imported as csvs from shared Google Drive L1 folder
# DATA OUTPUT:
# PROJECT:       LTER Grassland Rock synthesis group
# DATE:          2/15/2024

# Set-up ----
# Clear all existing data
rm(list=ls())

#Load packages
library(tidyverse)
#library(janitor)

# Set working directory 
L0_dir <- Sys.getenv("L0DIR")
L1_dir <- Sys.getenv("L1DIR")
list.files(L1_dir)



kbsdata <- read.csv(file.path(L1_dir, "KBS_MCSE_GLBRC_ANPP_RICH.csv"))

kbs_species <- read.csv(file.path(L1_dir, "KBS_MCSE_GLBRC_SpComp.csv"))

head(kbsdata)
head(kbs_species)
table(kbs_species$species)



# calculate shannon diversity   - based on ashleys code
head(kbs_species)
kbs_species <- kbs_species %>%
  mutate(proportions = (pseudo_perccover/100)*(log(pseudo_perccover/100)))
kbs_species$proportions[is.na(kbs_species$proportions)] <- 0  # Fix Nas to zeros
kbs_shannon <- kbs_species %>%
  group_by(year, treatment, station, replicate, disturbance, nutrients_added) %>%
  summarise(shannon = -1*sum(proportions)) # single value for each plot
kbs_shannon

names(kbsdata)
names(kbs_shannon)

# merge this diversity data with ANPP data
kbsdata_anpp_div <- merge (kbsdata, kbs_shannon, by = c("year", "treatment", "station", "replicate" ,"nutrients_added", "disturbance"))

#calculate evenness
kbsdata_anpp_div$evenness <- kbsdata_anpp_div$shannon / log (kbsdata_anpp_div$plot_richness)


idcols <- c("year", "treatment", "station", "replicate" ,"nutrients_added", "disturbance")
plotid <- c( "treatment", "station", "replicate" ,"nutrients_added", "disturbance")

kbsdata_anpp_div$unique_id <- apply( kbsdata_anpp_div[ , idcols ] , 1 , paste , collapse = "_" )
kbsdata_anpp_div$plot_id <- apply( kbsdata_anpp_div[ , plotid ] , 1 , paste , collapse = "_" )
kbsdata_anpp_div$unique_id
kbsdata_anpp_div$plot_id


# questions - 

# across all expts and trts does diversity increase ANPP?

ggplot(kbsdata_anpp_div, aes (x = plot_richness, y = plot_biomass)) +  # richness
  geom_point() +
  geom_smooth() + 
  xlab("Richness") + ylab("Biomass (g/m2)") +
  theme_classic () 

ggplot(kbsdata_anpp_div, aes (x = evenness, y = plot_biomass)) +  # evenness
  geom_point() +
  geom_smooth() + 
  xlab("Evenness") + ylab("Biomass (g/m2)") +
  theme_classic () 


ggplot(kbsdata_anpp_div, aes (x = shannon, y = plot_biomass)) +  # div
  geom_point() +
  geom_smooth() + 
  xlab("Shannon") + ylab("Biomass (g/m2)") +
  theme_classic () 




# look at 2012 dr
kbs_anpp_div_2012dr <- kbsdata_anpp_div %>% 
  filter(year > 2010 & year < 2014)
head(kbs_anpp_div_2012dr)
table(kbs_anpp_div_2012dr$year)

# average biomass per richness for each year? i dont know


kbs_anpp_div_2012dr_wide <- kbs_anpp_div_2012dr %>% 
  select (c (year, plot_id, plot_biomass)) %>% 
  pivot_wider(names_from= year, values_from = plot_biomass) %>% 
  mutate (perc_change_dr = ((`2012`- `2011`) / `2011`) * 100) %>% 
  mutate (perc_change_recov = ((`2013`- `2012`) / `2012`) * 100)


kbs_rich_2011 <- kbs_anpp_div_2012dr %>% 
  filter(year == 2011)

kbs_anpp_div_2012dr_wide_wrich <- merge (kbs_anpp_div_2012dr_wide, kbs_rich_2011 , by = "plot_id")
 
ggplot(kbs_anpp_div_2012dr_wide_wrich, aes (x = plot_richness, y = perc_change_dr )) +  # evenness
  geom_point(aes (color = as.character(year)) )+
  scale_color_manual(values = c("blue", "red","green")) + 
  geom_hline(yintercept=0, col = "red") + 
  geom_smooth() + 
  xlab("rich") + ylab("perc change with drought")
  





