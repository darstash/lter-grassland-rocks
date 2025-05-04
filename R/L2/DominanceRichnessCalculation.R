# LTER Grassland Rock: Calculate richness and dominance
# AUTHORS: ...,Joshua Ajowele
# COLLABORATORS:  
# DATA INPUT: .csv file of plot metrics and species abundance imported from Google Drive L2 folder
# DATA OUTPUT: .csv of updated plot metrics dataset that includes richness and dominance  
# PROJECT: LTER Grassland Rock
# DATE: October 2024: Last updated: May, 2025

# Clear all existing data
rm(list=ls())

#a community metric package that provides options for various evenness metrics.
#install.packages("codyn")
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
  select(year, uniqueid, species, species_harmonized,relative_abundance)%>%
  distinct()
#calculates evenness and richness
Evenness_richness<-community_structure(species_abundance_SPEI_evar, time.var="year",
                                       replicate.var = "uniqueid",
                                       abundance.var="relative_abundance",metric="Evar")%>%
  select(-richness)#remove richness since it was already calculated
#monoculture produces NAs for evenness


#determine the species with the average max relative abundance within each plot across all years
dominant_species_year<-species_abundance_SPEI_evar%>%
  group_by(uniqueid, year)%>%
  mutate(relative_abundance_max=max(relative_abundance))%>%
  filter(relative_abundance==relative_abundance_max)%>%#filter dominant species in each plot
  group_by(uniqueid, species, species_harmonized)%>%
  mutate(count=n())%>%#frequency a species was the dominant species in each plot
  group_by(uniqueid)%>%
  mutate(count_max=max(count))%>%
  filter(count==count_max)%>%#filter the species with the highest frequency as dominant species in each plot
  mutate(dominant_species_code=paste(uniqueid, species_harmonized,sep="_"))%>%
  select(uniqueid, species,species_harmonized, dominant_species_code)%>%
  distinct()

#create dataframe of relative abundance of dominant species
rel_abund_dom_species_year<-species_abundance_SPEI_evar%>%
  left_join(dominant_species_year, by=c("uniqueid","species","species_harmonized"))%>%
  filter(dominant_species_code!="NA")%>%#retains only the relative abundance of the dominant species
  rename(dominant_relative_abund=relative_abundance)%>%
  select(-species,-species_harmonized)


plot_metrics_SPEI_diversity <- plot_metrics_SPEI %>% 
  right_join(., species_abundance_SPEI_Metric, by = c("year", "site",  "higher_order_organization", "plot", "uniqueid",
                                                     "spei12", "spei3", "spei6", "spei9", "spei6_category", "spei12_category"))%>%
  left_join(Evenness_richness, by = c("year","uniqueid"))%>%
  left_join(rel_abund_dom_species_year, by=c("year","uniqueid"))%>%
  #create a new column of dominant abundance and convert NAs to zero
  #as NAs are primarily due to the absence of the dominant species
  mutate(dominant_relative_abund_zero=dominant_relative_abund)%>%
  mutate(dominant_relative_abund_zero=ifelse(is.na(dominant_relative_abund_zero),0,dominant_relative_abund_zero))%>%
  #create a column for predictors (richness, etc.) prior to the climate events
  group_by(uniqueid)%>%
  arrange(uniqueid, year)%>%
  mutate(prior_year_rich=lag(Richness),
         prior_year_type=lag(spei9_category),
         prior_year_dom_zero=lag(dominant_relative_abund_zero),
         prior_year_dom=lag(dominant_relative_abund),
         prior_year_evar=lag(Evar),
         prior_year_spei9=lag(spei9))


write.csv(plot_metrics_SPEI_diversity, file.path(L2_dir, "./plot_metrics_SPEI_diversity.csv"), row.names=F)

#investigating NAs in dominant species abundance
plot_Na<-plot_metrics_SPEI_diversity%>%
  select(year,uniqueid,site, dominant_relative_abund, dominant_species_code)%>%
  filter(uniqueid=="004b_D_tu_1")
species_na<-species_abundance_SPEI_evar%>%
  select(year,uniqueid, species)%>%
  filter(uniqueid=="004b_D_tu_1")%>%
  filter(species=="Cornus drummondii")
#NA occurred primarily due to the absence of the dominant species; about 1200 out of 21648 observations affected
#potential solutions: leave as NAs or convert NAs to zero
