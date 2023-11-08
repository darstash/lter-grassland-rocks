# TITLE:          LTER Grassland Rock: KBS ANPP biomass and plant composition cleanup
# AUTHORS:        
# COLLABORATORS:  
# DATA INPUT:     Data imported as csv files from shared Google drive L0 folder
# DATA OUTPUT:    
# PROJECT:        LTER Grassland Rock
# DATE:           

# Clear all existing data
rm(list=ls())

#Load packages
library(tidyverse)

# Set working directory 
Sys.getenv("L0DIR")
L0_dir <- Sys.getenv("L0DIR")
L1_dir <- Sys.getenv("L1DIR")

# EXAMPLE