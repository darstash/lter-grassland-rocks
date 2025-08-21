# TITLE:        LTER Grassland Rock: All data assembly
# AUTHORS:      Seraina Cappelli
# COLLABORATORS:  
# DATA INPUT:   CDR_initial_data_wrangling_L1.R, 
#               KNZ_initial_data_wrangling_L1.Rmd, 
#               KBS_initial_data_wrangling_L1.R, 
#               all_sites_dataset.R, 
#               all_sites_dataset_plus_spei.R,
#               dominance_diversity_calculation.R, 
#               resistance_resilience_calculation.R
# DATA OUTPUT:    
# PROJECT:      LTER Grassland Rock
# DATE:         July 2024

# This code is to run all relevant data assembly code in one whoosh, so that
# one doesn't have to open each processing step separately and run it manually
# when one of the baseline datasets has changed.

# Clear all existing data
rm(list=ls())


# site data cleaning ####
#-----------------------#
source("R/L1/CDR_initial_data_wrangling_L1.R")
knitr::knit("R/L1/KNZ_initial_data_wrangling_L1.Rmd") # creates some csv and other funky stuff. only run it if you know you need it
source("R/L1/KBS_initial_data_wrangling_L1.R")


# combine data together and add SPEI ####
#---------------------------------------#
source("R/L1/all_sites_dataset.R")
source("R/L1/all_sites_dataset_plus_spei.R")


# calculate dominance, diversity (richness and evenness), resistance and resilience ####
#-------------------------------------------------------------#
source("R/L2/dominance_diversity_calculation.R")
source("R/L2/resistance_resilience_calculation.R") # be patient. the loop to calculate resistance and resilience is not the fastest.