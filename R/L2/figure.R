# TITLE:        LTER Grassland Rock: Create core figures 
# AUTHORS:      Ashley Darst
# COLLABORATORS:  
# DATA INPUT:   Data imported as csv files from shared Google drive L2 folder
# DATA OUTPUT:  Core figures
# PROJECT:      LTER Grassland Rock
# DATE:         October 2024

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
plot <- read.csv(file.path(L2_dir, "plot_metrics_SPEI_diversity.csv"))

# Figure 1:

