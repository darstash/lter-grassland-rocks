# TITLE:        LTER Grassland Rock: Combine SPEI with site 
# AUTHORS:      Caitlin Broderick
# COLLABORATORS:  
# DATA INPUT:   Data imported as csv files from shared Google drive L1 (spei) and L2 (metrics, abundance) folders
# DATA OUTPUT:  plot metrics and sp abundance that have SPEI
# PROJECT:      LTER Grassland Rock
# DATE:         AAugust 2024




rm(list=ls())

# Load packaages
library(tidyverse)
library(lubridate)

# Set working directory 
L0_dir <- Sys.getenv("L0DIR")
L1_dir <- Sys.getenv("L1DIR")
L2_dir <- Sys.getenv("L2DIR")
list.files(L0_dir)
list.files(L1_dir)
list.files(L2_dir)

# want to combine with plot metrics, and species abundance

spei <- read.csv(file.path(L1_dir, "SPEI_12_allsites.csv"), stringsAsFactors = FALSE)

spei <- read.csv(file.path(L1_dir, "SPEI_12_allsites.csv"), stringsAsFactors = FALSE)
plot <- read.csv(file.path(L2_dir, "plot_metrics.csv"), stringsAsFactors = FALSE)
sp <- read.csv(file.path(L2_dir, "species_abundance.csv"), stringsAsFactors = FALSE)

range(spei$year)
range(plot$year) # something wrong with KNZ years
plot_spei <- merge(plot, spei, by = c("year", "site")) # # hmm why losing observations....
  # losing observations from 2023 and some KNZ rows with bad years
sp_spei <- merge(sp, spei, by = c("year", "site")) # # hmm why losing observations....

write.csv(plot_spei, file.path(L2_dir, "./plot_metrics_SPEI.csv"), row.names=F)
write.csv(sp_spei, file.path(L2_dir, "./species_abundance_SPEI.csv"), row.names=F)
