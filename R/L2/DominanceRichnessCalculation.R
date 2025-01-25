# LTER Grassland Rock: Calculate richness and dominance
# AUTHORS: ...,Joshua Ajowele
# COLLABORATORS:  
# DATA INPUT: .csv file of plot metrics and species abundance imported from Google Drive L2 folder
# DATA OUTPUT: .csv of updated plot metrics dataset that includes richness and dominance  
# PROJECT: LTER Grassland Rock
# DATE: October 2024: Last updated: Jan 25, 2025

# Clear all existing data
rm(list=ls())

#a community metric package that provides options for various evenness metrics.
install.packages("codyn")
library(codyn)
library(tidyverse)


L2_dir <- Sys.getenv("L2DIR")

list.files(L2_dir)

species_abundance_SPEI <- read.csv(file.path(L2_dir, "species_abundance_SPEI.csv"), stringsAsFactors = F)
plot_metrics_SPEI <- read.csv(file.path(L2_dir, "plot_metrics_SPEI.csv"), stringsAsFactors = F)

species_abundance_SPEI_Metric <- species_abundance_SPEI %>% 
  filter(relative_abundance > 0) %>% 
  group_by(year, site, dataset, plot, higher_order_organization, uniqueid,
           cover_method,
           spei12, spei3, spei6, spei9, spei6_category, spei12_category) %>% 
  summarize(Berger_Parker = max(relative_abundance),
            Richness = n())

ggplot(species_abundance_SPEI_Metric, aes(x = Berger_Parker, y = Richness))+
  geom_point(color = "black", alpha = 1/3)+
  theme_bw()

#calculate evenness using EVar; Evar is independent of species richness and recommended by Smith and Wilson (1996)
species_abundance_SPEI_evar<-species_abundance_SPEI%>%
  filter(relative_abundance > 0)%>%
  select(year, uniqueid, species, relative_abundance)%>%
  distinct()
#calculates evenness and richness
Evenness_richness<-community_structure(species_abundance_SPEI_evar, time.var="year",
                                       replicate.var = "uniqueid",
                                       abundance.var="relative_abundance",metric="Evar")%>%
  select(-richness)#remove richness since it was already calculated
#monoculture produces NAs for evenness

#determine the species with the average max relative abundance within each plot across all years
dominant_species<-species_abundance_SPEI_evar%>%
  group_by(uniqueid, species)%>%
  summarise(rel_abund=mean(relative_abundance, na.rm=T))%>%
  filter(rel_abund==max(rel_abund))%>%
  mutate(dominant_species_code=paste(uniqueid, species,sep="_"))%>%#to be used to select the relative abudnance of the dominant species in each plot
  select(-rel_abund)

#create dataframe of relative abundance of dominant species
rel_abund_dom_species<-species_abundance_SPEI_evar%>%
  left_join(dominant_species, by=c("uniqueid","species"))%>%
  filter(dominant_species_code!="NA")%>%#retains only the relative abundance of the dominant species
  rename(dominant_relative_abund=relative_abundance)%>%
  select(-species)
  

plot_metrics_SPEI_diversity <- plot_metrics_SPEI %>% 
  right_join(., species_abundance_SPEI_Metric, by = c("year", "site",  "higher_order_organization", "plot", "uniqueid",
                                                     "spei12", "spei3", "spei6", "spei9", "spei6_category", "spei12_category"))%>%
  left_join(Evenness_richness, by = c("year","uniqueid"))%>%
  left_join(rel_abund_dom_species, by=c("year","uniqueid"))#need to find out why NA for some dominant species relative abundance


write.csv(plot_metrics_SPEI_diversity, file.path(L2_dir, "./plot_metrics_SPEI_diversity.csv"), row.names=F)

