# TITLE:        LTER Grassland Rock: Data checking (e054 field 26)
# AUTHORS:      Seraina Cappelli
# COLLABORATORS: Max Zaret, Ashley Darst
# DATA INPUT:   e054 and all CDR, all wroking group data created and manipulated in CDR_explore_L1.R and all_sites_dataset.R
# DATA OUTPUT:  N/A
# PROJECT:      LTER Grassland Rock
# DATE:         01 November 2026

# Ashley noted: "CDR, uniqueid Experiment 54  field 26 only has one year in the 
# dataset. Could be a typo in the name @Seraina Cappelli @Max Zaret"

# checking CDR_explore_L1.R ####

e054_anpp <- read.csv(paste(L0_dir, "CDR/e54_biomass_1221_ML.csv", sep = "/")) 

e054_anpp %>%
  filter(OldField %in% 26) %>%
  str()
# 1007 entrys

e054_anpp %>%
  filter(OldField %in% 26) %>%
  select(Exp, Year, OldField) %>%
  unique() %>%
  nrow()
# 35 Years

e054_anpp %>%
  filter(OldField %in% 26) %>%
  select(Exp, OldField, Transect, Plot) %>%
  unique() %>%
  nrow()
# 5 plots (one per transect. Always plot 1, transect G, R, W, Y, y)


e054_anpp %>%
  filter(OldField %in% 26) %>%
  select(Exp, Year, OldField, Transect, Plot) %>%
  unique() %>%
  nrow()
# 140 Year x plot

e054_anpp %>%
  filter(OldField %in% 26) %>%
  select(Exp, Year, OldField, Transect) %>%
  unique %>%
  nrow()
# 140 Year x plot


#------------------------------------------------------------------------------

# in line 700 things get renamed

e054_anpp %>%
  filter(higher_order_organization %in%  "Experiment 54  field 26") %>%
  data.frame() %>%
  str()
# now justified 1006, then 1005 entries just

e054_anpp %>%
  filter(higher_order_organization %in%  "Experiment 54  field 26") %>%
  data.frame() %>%  
  select(year) %>%
  unique() %>%
  nrow()
# still 35 Years

e054_anpp %>%
  filter(higher_order_organization %in%  "Experiment 54  field 26") %>%
  data.frame() %>%  
  select(plot, higher_order_organization, uniqueid) %>%
  unique() %>%
  nrow()
# 5 plots

e054_anpp %>%
  filter(higher_order_organization %in%  "Experiment 54  field 26") %>%
  data.frame() %>%  
  select(year, plot, higher_order_organization, uniqueid) %>%
  unique() %>%
  nrow()
# 140 combinations (shouldn't it be 136)

#-----------------------------------------------------------------------------#

# metrics
e054_metrics %>%
  data.frame() %>%
  filter(higher_order_organization %in%  "Experiment 54  field 26") %>%
  str()
# 140 observations

#-----------------------------------------------------------------------------#

cdr_data %>%
  data.frame()%>%
  filter(higher_order_organization %in%  "Experiment 54  field 26") %>%
  str()
# 140 observations

cdr_sp_data %>%
  data.frame() %>%
  filter(higher_order_organization %in%  "Experiment 54  field 26")%>%
  str()
# 1005 obseverations

#_____________________________________________________________________________#
# checking all_sites_dataset.R ####
CDR_species_abundance %>%
  filter(higher_order_organization %in%  "Experiment 54  field 26")%>%
  str()
# 1005 observations

CDR_plot_metrics%>%
  filter(higher_order_organization %in%  "Experiment 54  field 26")%>%
  str()
# 140 observations

plot_metrics %>%
  filter(higher_order_organization %in%  "Experiment 54  field 26")%>%
  str()
# 140 observations

species_abundance %>%
  filter(higher_order_organization %in%  "Experiment 54  field 26")%>%
  str()
# 1005 observations

species_abundance %>%
  filter(higher_order_organization %in%  "Experiment 54  field 26")%>%
  select(year) %>%
  unique() %>%
  nrow()
# 35 years

plot_metrics %>%
  filter(higher_order_organization %in%  "Experiment 54  field 26")%>%
  select(year) %>%
  unique() %>%
  nrow()
# 35 years

species_abundance %>%
  filter(higher_order_organization %in%  "Experiment 54  field 26")%>%
  summary()

plot_metrics %>%
  filter(higher_order_organization %in%  "Experiment 54  field 26")%>%
  summary()

#_____________________________________________________________________________#
# checking all_sites_dataset_plus_spei.R ####
plot_spei %>%
  filter(higher_order_organization %in%  "Experiment 54  field 26")%>%
  str()
# 140 observations

sp_spei %>%
  filter(higher_order_organization %in%  "Experiment 54  field 26")%>%
  str()
# 1005 observations

sp_spei %>%
  filter(higher_order_organization %in%  "Experiment 54  field 26")%>%
  select(year) %>%
  unique() %>%
  nrow()
# 35 years

plot_spei %>%
  filter(higher_order_organization %in%  "Experiment 54  field 26")%>%
  select(year) %>%
  unique() %>%
  nrow()
# 35 years

sp_spei %>%
  filter(higher_order_organization %in%  "Experiment 54  field 26")%>% 
  summary()

plot_spei %>%
  filter(higher_order_organization %in%  "Experiment 54  field 26")%>%
  summary()

#_____________________________________________________________________________#
# checking DominanceRichnessCalculation.R ####
plot_metrics_SPEI_diversity %>%
  filter(higher_order_organization %in%  "Experiment 54  field 26")%>%
  str()
# 140 entries

plot_metrics_SPEI_diversity %>%
  filter(higher_order_organization %in%  "Experiment 54  field 26")%>%
  select(year) %>%
  unique() %>%
  nrow()
# 35 years

plot_metrics_SPEI_diversity %>%
  filter(higher_order_organization %in%  "Experiment 54  field 26")%>%
  select(uniqueid) %>%
  unique()
# 5 plots


species_abundance_SPEI_Metric%>%
  data.frame() %>%
  filter(higher_order_organization %in%  "Experiment 54  field 26")%>%
  str()
# 140 entries

species_abundance_SPEI_Metric%>%
  data.frame() %>%
  filter(higher_order_organization %in%  "Experiment 54  field 26")%>%
  select(year) %>%
  unique() %>%
  nrow()
# 35 years

species_abundance_SPEI_Metric%>%
  data.frame() %>%
  filter(higher_order_organization %in%  "Experiment 54  field 26")%>%
  select(uniqueid) %>%
  unique()
# 5 plots