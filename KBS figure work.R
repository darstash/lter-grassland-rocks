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
# ID cols - plot ID sampled at particular time point (so basically row ID)
# plot ID - same space sampled year after year


kbsdata_anpp_div$unique_id <- apply( kbsdata_anpp_div[ , idcols ] , 1 , paste , collapse = "_" )
kbsdata_anpp_div$plot_id <- apply( kbsdata_anpp_div[ , plotid ] , 1 , paste , collapse = "_" )
kbsdata_anpp_div$unique_id
kbsdata_anpp_div$plot_id

# Calculate stability of each plot across ALL years
stability <- kbsdata_anpp_div %>%
  group_by(treatment, station, replicate, nutrients_added, disturbance) %>%
  summarize(mean_anpp = mean(plot_biomass),
            sd_anpp = sd(plot_biomass),
            mean_richness = mean(plot_richness),
            mean_shannon = mean(shannon))
stability <- stability %>%
  mutate(stability = mean_anpp/sd_anpp)

# super basic questions - 

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

# richness over time per plot (richness generally increases initially, hasn't plateaued for newer plots)
# problem for "normalizing" richness?
ggplot(kbsdata_anpp_div, aes(x = year, y = plot_richness)) +
  geom_point() +
  facet_wrap(~plot_id)


# top drought years seem like 1999, 2003, 2005, 2012, and 1996 (based on SPEI-12)
# results below are all over the place

################################
# look at 2012 drought
################################
kbs_anpp_div_2012dr <- kbsdata_anpp_div %>% 
  filter(year > 2010 & year < 2014) # subset 2011, 2012, 2013 data (before, during, after)
head(kbs_anpp_div_2012dr)
table(kbs_anpp_div_2012dr$year)



# convert to wide format 
kbs_anpp_div_2012dr_wide <- kbs_anpp_div_2012dr %>% 
  select (c (year, plot_id, plot_biomass)) %>%  # select columns I will use
  pivot_wider(names_from= year, values_from = plot_biomass) %>%  # pivot wider so each year is own column
  mutate (perc_change_dr = ((`2012`- `2011`) / `2011`) * 100) %>% # add column for perc change during drought
  mutate (perc_change_recov = ((`2013`- `2011`) / `2011`) * 100) %>% # add column for perc recovery back to 2011
  mutate (resilience = ((`2012`- `2011`) / (`2013` - `2011`))) %>% # add column for resilience (as defined by Isbell 2015)
  mutate (resistance = ((`2011`) / (`2012` - `2011`))) # add column for resilience (as defined by Isbell 2015)


kbs_rich_2011 <- kbs_anpp_div_2012dr %>% 
  filter(year == 2011) # just get richness for 2011

# merger 2011 richness for each plot with the wide-format biomass data fror 
kbs_anpp_div_2012dr_wide_wrich <- merge (kbs_anpp_div_2012dr_wide, kbs_rich_2011 , by = "plot_id")
 

# drought resistance
ggplot(kbs_anpp_div_2012dr_wide_wrich, aes (x = plot_richness, y = perc_change_dr )) +  
  geom_point(aes (color = as.character(year)) )+
  scale_color_manual(values = c("blue", "red","green")) + 
  geom_hline(yintercept=0, col = "red") + 
  geom_smooth(method = "lm") + 
  xlab("rich") + ylab("perc change with drought") 
  # can add facet by treatment here! such as nutrients_added
  # does N addition change sensitivity to dr
  
# drought recovery
ggplot(kbs_anpp_div_2012dr_wide_wrich, aes (x = plot_richness, y = perc_change_recov )) + 
  geom_point(aes (color = as.character(year)) )+
  scale_color_manual(values = c("blue", "red","green")) + 
  geom_hline(yintercept=0, col = "red") + 
  geom_smooth(method = "lm") + 
  xlab("rich") + ylab("perc recovery") 

# Resistance
ggplot(kbs_anpp_div_2012dr_wide_wrich, aes (x = plot_richness, y = resistance)) + 
  geom_point()+
  geom_hline(yintercept=0, col = "red") + 
  geom_smooth(method = "lm") + 
  xlab("rich") + ylab("resistance") 

# Resilience
ggplot(kbs_anpp_div_2012dr_wide_wrich, aes (x = plot_richness, y = resilience)) + 
  geom_point()+
  geom_hline(yintercept=0, col = "red") + 
  geom_smooth(method = "lm") + 
  xlab("rich") + ylab("resilience") 


################################
# look at 1999 drought
################################
kbs_anpp_div_1999dr <- kbsdata_anpp_div %>% 
  filter(year > 1997 & year < 2001) # subset 1998, 1999, 2000 data (before, during, after)
head(kbs_anpp_div_1999dr)
table(kbs_anpp_div_1999dr$year) # less for 1998?



# convert to wide format 
kbs_anpp_div_1999dr_wide <- kbs_anpp_div_1999dr %>% 
  select (c (year, plot_id, plot_biomass)) %>%  # select columns I will use
  pivot_wider(names_from= year, values_from = plot_biomass) %>%  # pivot wider so each year is own column
  mutate (perc_change_dr = ((`1999`- `1998`) / `1998`) * 100) %>% # add column for perc change during drought
  mutate (perc_change_recov = ((`2000`- `1998`) / `1998`) * 100) %>% # add column for perc recovery back to 2011
  mutate (resilience = ((`1999`- `1998`) / (`2000` - `1998`))) %>% # add column for resilience (as defined by Isbell 2015)
  mutate (resistance = ((`1998`) / (`1999` - `1998`))) # add column for resilience (as defined by Isbell 2015)


kbs_rich_1998 <- kbs_anpp_div_1999dr %>% 
  filter(year == 1998) # just get richness for 1998

# merger 2011 richness for each plot with the wide-format biomass data fror 
kbs_anpp_div_1999dr_wide_wrich <- merge (kbs_anpp_div_1999dr_wide, kbs_rich_1998 , by = "plot_id")
table(kbs_anpp_div_1999dr_wide_wrich$year) # dropped three plots


# drought resistance
ggplot(kbs_anpp_div_1999dr_wide_wrich, aes (x = plot_richness, y = perc_change_dr )) +  
  geom_point()+
  geom_hline(yintercept=0, col = "red") + 
  geom_smooth(method = "lm") + 
  xlab("rich") + ylab("perc change with drought") 
# can add facet by treatment here! such as nutrients_added
# does N addition change sensitivity to dr

# drought recovery
ggplot(kbs_anpp_div_1999dr_wide_wrich, aes (x = plot_richness, y = perc_change_recov )) + 
  geom_point()+
  geom_hline(yintercept=0, col = "red") + 
  geom_smooth(method = "lm") + 
  xlab("rich") + ylab("perc recovery") 

# Resistance
ggplot(kbs_anpp_div_1999dr_wide_wrich, aes (x = plot_richness, y = resistance)) + 
  geom_point()+
  geom_hline(yintercept=0, col = "red") + 
  geom_smooth(method = "lm") + 
  xlab("rich") + ylab("resistance") 

# Resilience
ggplot(kbs_anpp_div_1999dr_wide_wrich, aes (x = plot_richness, y = resilience)) + 
  geom_point()+
  geom_hline(yintercept=0, col = "red") + 
  geom_smooth(method = "lm") + 
  xlab("rich") + ylab("resilience") 


################################
# look at 2003 drought
################################
kbs_anpp_div_2003dr <- kbsdata_anpp_div %>% 
  filter(year > 2001 & year < 2005) # subset 2002, 2003, 2004 data (before, during, after)
head(kbs_anpp_div_2003dr)
table(kbs_anpp_div_2003dr$year) # all 54 plots



# convert to wide format 
kbs_anpp_div_2003dr_wide <- kbs_anpp_div_2003dr %>% 
  select (c (year, plot_id, plot_biomass)) %>%  # select columns I will use
  pivot_wider(names_from= year, values_from = plot_biomass) %>%  # pivot wider so each year is own column
  mutate (perc_change_dr = ((`2003`- `2002`) / `2002`) * 100) %>% # add column for perc change during drought
  mutate (perc_change_recov = ((`2004`- `2002`) / `2002`) * 100) %>% # add column for perc recovery back to 2011
  mutate (resilience = ((`2003`- `2002`) / (`2004` - `2002`))) %>% # add column for resilience (as defined by Isbell 2015)
  mutate (resistance = ((`2002`) / (`2003` - `2002`))) # add column for resilience (as defined by Isbell 2015)


kbs_rich_2002 <- kbs_anpp_div_2003dr %>% 
  filter(year == 2002) # just get richness for 2002

# merger 2011 richness for each plot with the wide-format biomass data fror 
kbs_anpp_div_2003dr_wide_wrich <- merge (kbs_anpp_div_2003dr_wide, kbs_rich_2002 , by = "plot_id")
table(kbs_anpp_div_2003dr_wide_wrich$year) # all plots kept


# drought resistance
ggplot(kbs_anpp_div_2003dr_wide_wrich, aes (x = plot_richness, y = perc_change_dr )) +  
  geom_point()+
  geom_hline(yintercept=0, col = "red") + 
  geom_smooth(method = "lm") + 
  xlab("rich") + ylab("perc change with drought") 
# can add facet by treatment here! such as nutrients_added
# does N addition change sensitivity to dr

# drought recovery
ggplot(kbs_anpp_div_2003dr_wide_wrich, aes (x = plot_richness, y = perc_change_recov )) + 
  geom_point()+
  geom_hline(yintercept=0, col = "red") + 
  geom_smooth(method = "lm") + 
  xlab("rich") + ylab("perc recovery") 

# Resistance # weird outlier
ggplot(kbs_anpp_div_2003dr_wide_wrich, aes (x = plot_richness, y = resistance)) + 
  geom_point()+
  geom_hline(yintercept=0, col = "red") + 
  geom_smooth(method = "lm") + 
  xlab("rich") + ylab("resistance") 

# Resilience
ggplot(kbs_anpp_div_2003dr_wide_wrich, aes (x = plot_richness, y = resilience)) + 
  geom_point()+
  geom_hline(yintercept=0, col = "red") + 
  geom_smooth(method = "lm") + 
  xlab("rich") + ylab("resilience") 


################################
# look at stability
################################

# All years stability 
# weird, opposite of expected
stability %>%
  ggplot(aes(x = mean_richness, y = stability, col = nutrients_added)) +
  geom_point() +
  geom_smooth(method = "lm")
stability %>%
  ggplot(aes(x = mean_richness, y = stability, col = nutrients_added)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~disturbance)

