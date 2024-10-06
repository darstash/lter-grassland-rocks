rm(list=ls())
# Name: Calculation of Berger-Parker and diversity each year, uniqueid, site
# Author: Maowei Liang
# Date: 10-05-2024

########################################################## 
# Note: Calculation runs a couple of hours if using laptop
##########################################################

# load the packages
library(dplyr)
library(tidyverse)
library(tidyr)
library(ggpmisc)
library(gridExtra)

# read the datasets
spe_abun.data <- read.csv('species_abundance.csv',header=T)
variable.names(spe_abun.data)

#spe_abun_CDR_KBS.data <- subset(spe_abun.data, site %in% c("CDR", "KBS"))

#spe_abun_CDR.data <- subset(spe_abun.data, site == "CDR")
#spe_abun_KBS.data <- subset(spe_abun.data, site == "KBS")
#spe_abun_KNZ.data <- subset(spe_abun.data, site == "KNZ") VERY LONG

dd.data <- spe_abun.data
# Extract unique sites, years, and uniqueid, sorted
site <- sort(unique(dd.data$site))
year <- sort(unique(dd.data$year))
uniqueid <- sort(unique(dd.data$uniqueid))

# Initialize the main results dataframe
di_do <- data.frame(site = character(), year = numeric(), uniqueid = character(),
                    berger_parker = numeric(),
                    richness = numeric(), 
                    invsimpson = numeric(), 
                    evenness = numeric(), 
                    total_abun = numeric())

# Nested loops to calculate diversity indices for each site, year, and uniqueid
for (ss in site) {
  for (yr in year) {
    for (uid in uniqueid) {
      # Filter data for current site, year, and uniqueid
      dd <- dd.data[dd.data$site == ss & dd.data$year == yr & dd.data$uniqueid == uid, ]
      
      # Remove NA values only from the abundance column
      dd <- dd[!is.na(dd$abundance), ]
      
      if (nrow(dd) > 0) {
        abundance <- dd$abundance
        total_abun <- sum(abundance)  # Calculate total abundance
        proportions <- abundance / total_abun
        
        # Berger-Parker index
        berger_parker <- max(abundance) / total_abun
        
        # Diversity indices
        invsimpson <- 1 / sum(proportions^2)
        richness <- length(unique(dd$species))
        evenness <- ifelse(richness > 0, invsimpson / richness, 0)
        
        # Append results for this site, year, and uniqueid
        di_do <- rbind(di_do, data.frame(site = ss, year = yr, uniqueid = uid,
                                         richness = richness, 
                                         invsimpson = invsimpson, 
                                         evenness = evenness,  
                                         total_abun = total_abun,
                                         berger_parker = berger_parker))
      }
    }
  }
}

#write.csv(di_do,"div_domi_10052024.csv")

# Calculate mean values for each variable in di_do grouped by site
mean_di_do <- di_do %>%
  group_by(site) %>%
  summarize(mean_richness = mean(richness, na.rm = TRUE),
            mean_invsimpson = mean(invsimpson, na.rm = TRUE),
            mean_evenness = mean(evenness, na.rm = TRUE),
            mean_total_abun = mean(total_abun, na.rm = TRUE),
            mean_berger_parker = mean(berger_parker, na.rm = TRUE))

# Convert mean_di_do to a long format to access mean values conveniently
mean_di_do_long <- mean_di_do %>%
  pivot_longer(cols = starts_with("mean_"), names_to = "variable", values_to = "mean_value") %>%
  mutate(variable = gsub("mean_", "", variable)) # Remove "mean_" prefix for easier matching

# Create an empty list to store the plots
plot_list <- list()

# Define the list of variables for which to create density plots
variables <- c("richness", "invsimpson", "evenness", "total_abun", "berger_parker")

# Loop through each variable to create density plots
for (var in variables) {
  # Extract the mean value for each site for the current variable
  mean_values <- mean_di_do_long %>% filter(variable == var)
  
  # Create plot data frame for the current variable
  plot_data <- data.frame(value = di_do[[var]], site = di_do$site)
  
  # Create density plot with ggplot
  p <- ggplot(plot_data, aes(x = value, fill = as.factor(site))) +
    geom_histogram(bins = 30, color = "black", alpha = 0.7) +
    scale_fill_manual(values = c('#D55E00','#0072B2','#009E73')) +
    labs(title = paste("", var), x = var, y = "Frequncy", fill = "Site") +
    theme_bw() +
    theme(legend.position = "bottom",
          axis.text = element_text(size = 12, family = "Arial", color = "black"),
          axis.title = element_text(size = 16, family = "Arial", color = "black"),
          axis.ticks.length = unit(1, "mm"),
          axis.line = element_line(color = "black", size = 0.5)) +
    # Add mean lines for each site
    geom_vline(data = mean_values, aes(xintercept = mean_value, color = site), linetype = "dashed", size = 1) +
    scale_color_manual(values = c('#D55E00', '#0072B2','#009E73'))
  
  # Add to the plot list
  plot_list[[var]] <- p
}

# Combine all plots into a single figure
grid.arrange(grobs = plot_list, ncol = 3)



