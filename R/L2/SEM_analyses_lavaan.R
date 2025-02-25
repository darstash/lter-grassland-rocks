# TITLE:        LTER Grassland Rock: Create lavaan based SEM analyses 
# AUTHORS:      Seraina Cappelli
# COLLABORATORS:
# DATA INPUT:   Data imported as csv files from shared Google drive L2 folder
# DATA OUTPUT:  SEM analyses
# PROJECT:      LTER Grassland Rock
# DATE:         February 2025

library(tidyverse)
# dplyr     1.1.4     # readr     2.1.5
# forcats   1.0.0     # stringr   1.5.1
# ggplot2   3.5.1     # tibble    3.2.1
# lubridate 1.9.4     # tidyr     1.3.1
# purrr     1.0.2     
library(janitor) # 2.2.0
library(lavaan) # 0.6-17
library(DHARMa) # 0.4.6


# prep data like in analyses.R ####

# Set working directory 
L0_dir <- Sys.getenv("L0DIR")
L1_dir <- Sys.getenv("L1DIR")
L2_dir <- Sys.getenv("L2DIR")
list.files(L2_dir)

# Read in CSV files
plot <- read.csv(file.path(L2_dir, "plot_metrics_SPEI_diversity.csv"))
ece <- read.csv(file.path(L2_dir, "ece_resist_resil.csv"))
meta <- read.csv(file.path(L2_dir, "metadata.csv"))

# Only keep distinct rows in ece
ece <- distinct(ece)

# Change ex_year to year
ece <- ece %>%
  rename(year = ex_year)

# Merge plot with resistance and resilience
plot_ece <- left_join(plot, ece)

# Make column with categories for high and low dominance
plot_ece <- plot_ece %>%
  mutate(dom_category = case_when(
    Berger_Parker >= 0.5 ~ "high",
    Berger_Parker < 0.5 ~ "low"
  ))

# Standardize column names
plot_ece <- clean_names(plot_ece)

# Convert Inf values to NA for resilience for some KBS 2015 plots
plot_ece$resilience[plot_ece$resilience == Inf] <- NA

# Add experiment column
plot_ece$experiment <- sub("nutnet.*", "nutnet", plot_ece$higher_order_organization)
plot_ece$experiment <- sub("glbrc_scaleup.*", "glbrc_scaleup", plot_ece$experiment)
plot_ece$experiment <- sub("glbrc_G10.*", "glbrc", plot_ece$experiment)
plot_ece$experiment <- sub("glbrc_G9.*", "glbrc", plot_ece$experiment)
plot_ece$experiment <- sub("mcse.*", "mcse", plot_ece$experiment)
plot_ece$experiment <- sub("microplots.*", "microplots", plot_ece$experiment)
plot_ece$experiment <- sub("Experiment 1.*", "Experiment 1", plot_ece$experiment)
plot_ece$experiment <- sub("001d.*", "001d", plot_ece$experiment)
plot_ece$experiment <- sub("004a.*", "004a", plot_ece$experiment)
plot_ece$experiment <- sub("004b.*", "004b", plot_ece$experiment)
plot_ece$experiment <- sub("002d.*", "002d", plot_ece$experiment)
plot_ece$experiment <- sub("Experiment 54.*", "Experiment 54", plot_ece$experiment)
plot_ece$experiment <- sub("KNZ_WAT01.*", "KNZ_WAT01", plot_ece$experiment)
plot_ece$experiment <- sub("002c.*", "002c", plot_ece$experiment)
plot_ece$experiment <- sub("e061.*", "e061", plot_ece$experiment)
plot_ece$experiment <- sub("e247.*", "e247", plot_ece$experiment)
plot_ece$experiment <- sub("e245.*", "e245", plot_ece$experiment)
plot_ece$experiment[plot_ece$experiment == "A"] <- "NGE"
plot_ece$experiment[plot_ece$experiment == "B"] <- "NGE"
plot_ece$experiment[plot_ece$experiment == "C"] <- "NGE"
plot_ece$experiment[plot_ece$experiment == "D"] <- "NGE"
plot_ece$experiment[plot_ece$experiment == "E"] <- "NGE"
plot_ece$experiment[plot_ece$experiment == "F"] <- "NGE"

# Merge with metadata
plot_ece_meta <- left_join(plot_ece, meta)

# Remove NAs for non-extreme years
plot_ece_rm_na <- plot_ece_meta %>%
  drop_na(resistance)

# Make year a factor
plot_ece_rm_na$year <- as.factor(plot_ece_rm_na$year)

# # Subset to only have control plots
# plot_ece_control <- plot_ece_rm_na %>%
#   filter(treatment == "control")
# 
# # Look if measurment scale cover matters for resilience
# plot_ece_control %>%
#   ggplot(aes(x = measurement_scale_cover, y = log(resilience))) +
#   geom_point()
# # Look if measurment scale cover matters for richness
# plot_ece_control %>%
#   ggplot(aes(x = measurement_scale_cover, y = richness)) +
#   geom_point()

# Try to scale richness and BP by experiment for control plots only
plot_ece_rm_na <- plot_ece_rm_na %>%
  group_by(experiment) %>%
  mutate(richness_scaled = c(scale(richness)),
         berger_parker_scaled = c(scale(berger_parker)),
         evar_scaled=c(scale(evar)),
         dominant_relative_abund_scaled=c(scale(dominant_relative_abund)),
         dominant_relative_abund_zero_scaled=c(scale(dominant_relative_abund_zero)))

# Look if measurment scale cover matters for richness after scaling
plot_ece_rm_na %>%
  ggplot(aes(x = measurement_scale_cover, y = richness_scaled)) +
  geom_point()#no longer matters
rich_area<-lm(richness_scaled~measurement_scale_cover, data=plot_ece_rm_na)
summary(rich_area) 
simres <- simulateResiduals(rich_area)
plot(simres)#looks good enough


# build SEM ####

# ok, wait how did you code fertilizer and burn?
sem1.0 <- 'resistance ~ richness_scaled + dominant_relative_abund_zero_scaled
           richness_scaled ~ fertilizer + fire_frequency
           dominant_relative_abund_zero_scaled ~ fertilizer + fire_frequency
           '
