# TITLE:        LTER Grassland Rock: Tabulate extreme events and sequences 
# AUTHORS:      Ashley Darst
# COLLABORATORS:  
# DATA INPUT:   Data imported as csv files from shared Google drive L2 (plot_metrics_SPEI) folder
# DATA OUTPUT:  Table with number of extreme events by type
# PROJECT:      LTER Grassland Rock
# DATE:         September 2024

# Clear all existing data
rm(list=ls())

# Load packages
library(tidyverse)
library(ggforce)
library(lme4)
library(lmerTest)
library(emmeans)
library(sjPlot)

# Set working directory 
L0_dir <- Sys.getenv("L0DIR")
L1_dir <- Sys.getenv("L1DIR")
L2_dir <- Sys.getenv("L2DIR")
list.files(L2_dir)

# Read in CSV files ----
plot_spei <- read.csv(file.path(L2_dir, "plot_metrics_SPEI.csv"))

# Make sure spei category is a factor
str(plot_spei)
plot_spei$spei_category <- factor(plot_spei$spei_category, levels = c("Extreme dry", "Moderate dry", "Normal", "Moderate wet", "Extreme wet"))

# Summarize spei by year
spei_summary <- plot_spei %>%
  group_by(year, spei_category, site) %>%
  summarize(spei12 = mean(spei12))

# Look at SPEI plots for each site in the range the data exist
spei_summary %>%
  filter(site == "CDR") %>%
  ggplot(aes (x = year, y = spei12)) +
  geom_line(alpha = 0.5) + 
  geom_point( aes(color = spei_category), size = 1.25) + theme_bw() +
  scale_x_continuous(breaks = seq(1980, 2025, by = 10)) +
  annotate( "text", label = "CDR",
            x = 2000, y = 2, size = 5, colour = "black") + 
  scale_color_manual(values = c("#F5191CFF", "#E78200FF", "#EAC728FF", "#81BB95FF", "#3B99B1FF"))

spei_summary %>%
  filter(site == "KBS") %>%
  ggplot(aes (x = year, y = spei12)) +
  geom_line(alpha = 0.5) + 
  geom_point( aes(color = spei_category), size = 1.25) + theme_bw() +
  scale_x_continuous(breaks = seq(1990, 2025, by = 10)) +
  annotate( "text", label = "KBS",
            x = 2000, y = 2, size = 5, colour = "black") + 
  scale_color_manual(values = c("#E78200FF", "#EAC728FF", "#81BB95FF", "#3B99B1FF"))


spei_summary %>%
  filter(site == "KNZ") %>%
  ggplot(aes(x = year, y = spei12)) +
  geom_line(alpha = 0.5) + 
  geom_point( aes(color = spei_category), size = 1.25) + theme_bw() +
  scale_x_continuous(breaks = seq(1990, 2025, by = 10)) +
  annotate( "text", label = "KNZ",
            x = 2000, y = 2, size = 5, colour = "black") + 
  scale_color_manual(values = c("#F5191CFF", "#E78200FF", "#EAC728FF", "#81BB95FF", "#3B99B1FF"))

# Calculate the number of extreme wet years
CDR_ew <- length(which(spei_summary$site == "CDR" & spei_summary$spei_category == "Extreme wet"))
KBS_ew <- length(which(spei_summary$site == "KBS" & spei_summary$spei_category == "Extreme wet"))
KNZ_ew <- length(which(spei_summary$site == "KNZ" & spei_summary$spei_category == "Extreme wet"))
extreme_wet <- c(CDR_ew, KBS_ew, KNZ_ew)
sites <- c("CDR", "KBS", "KNZ")

# Calculate the number of moderate wet years
CDR_mw <- length(which(spei_summary$site == "CDR" & spei_summary$spei_category == "Moderate wet"))
KBS_mw <- length(which(spei_summary$site == "KBS" & spei_summary$spei_category == "Moderate wet"))
KNZ_mw <- length(which(spei_summary$site == "KNZ" & spei_summary$spei_category == "Moderate wet"))
moderate_wet <- c(CDR_mw, KBS_mw, KNZ_mw)

# Calculate the number of extreme dry years
CDR_ed <- length(which(spei_summary$site == "CDR" & spei_summary$spei_category == "Extreme dry"))
KBS_ed <- length(which(spei_summary$site == "KBS" & spei_summary$spei_category == "Extreme dry"))
KNZ_ed <- length(which(spei_summary$site == "KNZ" & spei_summary$spei_category == "Extreme dry"))
extreme_dry <- c(CDR_ed, KBS_ed, KNZ_ed)

# Calculate the number of moderate dry years
CDR_md <- length(which(spei_summary$site == "CDR" & spei_summary$spei_category == "Moderate dry"))
KBS_md <- length(which(spei_summary$site == "KBS" & spei_summary$spei_category == "Moderate dry"))
KNZ_md <- length(which(spei_summary$site == "KNZ" & spei_summary$spei_category == "Moderate dry"))
moderate_dry <- c(CDR_md, KBS_md, KNZ_md)

# Calculate the number of extreme wet then extreme dry years
spei_summary_CDR <- spei_summary %>%
  filter(site == "CDR")
CDR_ewd <- sum(spei_summary_CDR$spei_category[-nrow(spei_summary_CDR)] == "Extreme wet" & spei_summary_CDR$spei_category[-1] == "Extreme dry") # This was written by ChatGPT
spei_summary_KBS <- spei_summary %>%
  filter(site == "KBS")
KBS_ewd <- sum(spei_summary_KBS$spei_category[-nrow(spei_summary_KBS)] == "Extreme wet" & spei_summary_KBS$spei_category[-1] == "Extreme dry")
spei_summary_KNZ <- spei_summary %>%
  filter(site == "KNZ")
KNZ_ewd <- sum(spei_summary_KNZ$spei_category[-nrow(spei_summary_KNZ)] == "Extreme wet" & spei_summary_KNZ$spei_category[-1] == "Extreme dry")

extreme_wet_dry <- c(CDR_ewd, KBS_ewd, KNZ_ewd)

# Calculate the number of extreme dry then extreme wet years
CDR_edw <- sum(spei_summary_CDR$spei_category[-nrow(spei_summary_CDR)] == "Extreme dry" & spei_summary_CDR$spei_category[-1] == "Extreme wet")
KBS_edw <- sum(spei_summary_KBS$spei_category[-nrow(spei_summary_KBS)] == "Extreme dry" & spei_summary_KBS$spei_category[-1] == "Extreme wet")
KNZ_edw <- sum(spei_summary_KNZ$spei_category[-nrow(spei_summary_KNZ)] == "Extreme dry" & spei_summary_KNZ$spei_category[-1] == "Extreme wet")

extreme_dry_wet <- c(CDR_edw, KBS_edw, KNZ_edw)

# Calculate the number of extreme wet then extreme wet years
CDR_eww <- sum(spei_summary_CDR$spei_category[-nrow(spei_summary_CDR)] == "Extreme wet" & spei_summary_CDR$spei_category[-1] == "Extreme wet")
KBS_eww <- sum(spei_summary_KBS$spei_category[-nrow(spei_summary_KBS)] == "Extreme wet" & spei_summary_KBS$spei_category[-1] == "Extreme wet")
KNZ_eww <- sum(spei_summary_KNZ$spei_category[-nrow(spei_summary_KNZ)] == "Extreme wet" & spei_summary_KNZ$spei_category[-1] == "Extreme wet")

extreme_wet_wet <- c(CDR_eww, KBS_eww, KNZ_eww)

# Calculate the number of extreme dry then extreme dry years
CDR_edd <- sum(spei_summary_CDR$spei_category[-nrow(spei_summary_CDR)] == "Extreme dry" & spei_summary_CDR$spei_category[-1] == "Extreme dry")
KBS_edd <- sum(spei_summary_KBS$spei_category[-nrow(spei_summary_KBS)] == "Extreme dry" & spei_summary_KBS$spei_category[-1] == "Extreme dry")
KNZ_edd <- sum(spei_summary_KNZ$spei_category[-nrow(spei_summary_KNZ)] == "Extreme dry" & spei_summary_KNZ$spei_category[-1] == "Extreme dry")

extreme_dry_dry <- c(CDR_edd, KBS_edd, KNZ_edd)

spei_table <- data.frame(sites, extreme_wet, moderate_wet, extreme_dry, moderate_dry, extreme_wet_dry, extreme_dry_wet, extreme_wet_wet, extreme_dry_dry)

# Not enough extreme events. Let's lump all moderate and extreme events into just wet and dry categories
spei_lump <- spei_summary %>%
  mutate(spei_category = ifelse(as.character(spei_category) == "Extreme dry", "Dry", as.character(spei_category))) %>%
  mutate(spei_category = ifelse(as.character(spei_category) == "Moderate dry", "Dry", as.character(spei_category))) %>%
  mutate(spei_category = ifelse(as.character(spei_category) == "Extreme wet", "Wet", as.character(spei_category))) %>%
  mutate(spei_category = ifelse(as.character(spei_category) == "Moderate wet", "Wet", as.character(spei_category)))

# Calculate the number of moderate/extreme wet years
CDR_w <- length(which(spei_lump$site == "CDR" & spei_lump$spei_category == "Wet"))
KBS_w <- length(which(spei_lump$site == "KBS" & spei_lump$spei_category == "Wet"))
KNZ_w <- length(which(spei_lump$site == "KNZ" & spei_lump$spei_category == "Wet"))
wet <- c(CDR_w, KBS_w, KNZ_w)

# Calculate the number of moderate/extreme dry years
CDR_d <- length(which(spei_lump$site == "CDR" & spei_lump$spei_category == "Dry"))
KBS_d <- length(which(spei_lump$site == "KBS" & spei_lump$spei_category == "Dry"))
KNZ_d <- length(which(spei_lump$site == "KNZ" & spei_lump$spei_category == "Dry"))
dry <- c(CDR_d, KBS_d, KNZ_d)

# Calculate the number of wet then dry years
spei_lump_CDR <- spei_lump %>%
  filter(site == "CDR")
CDR_wd <- sum(spei_lump_CDR$spei_category[-nrow(spei_lump_CDR)] == "Wet" & spei_lump_CDR$spei_category[-1] == "Dry")
spei_lump_KBS <- spei_lump %>%
  filter(site == "KBS")
KBS_wd <- sum(spei_lump_KBS$spei_category[-nrow(spei_lump_KBS)] == "Wet" & spei_lump_KBS$spei_category[-1] == "Dry")
spei_lump_KNZ <- spei_lump %>%
  filter(site == "KNZ")
KNZ_wd <- sum(spei_lump_KNZ$spei_category[-nrow(spei_lump_KNZ)] == "Wet" & spei_lump_KNZ$spei_category[-1] == "Dry")

wet_dry <- c(CDR_wd, KBS_wd, KNZ_wd)

# Calculate the number of wet then dry years
spei_lump_CDR <- spei_lump %>%
  filter(site == "CDR")
CDR_wd <- sum(spei_lump_CDR$spei_category[-nrow(spei_lump_CDR)] == "Wet" & spei_lump_CDR$spei_category[-1] == "Dry")
spei_lump_KBS <- spei_lump %>%
  filter(site == "KBS")
KBS_wd <- sum(spei_lump_KBS$spei_category[-nrow(spei_lump_KBS)] == "Wet" & spei_lump_KBS$spei_category[-1] == "Dry")
spei_lump_KNZ <- spei_lump %>%
  filter(site == "KNZ")
KNZ_wd <- sum(spei_lump_KNZ$spei_category[-nrow(spei_lump_KNZ)] == "Wet" & spei_lump_KNZ$spei_category[-1] == "Dry")

wet_dry <- c(CDR_wd, KBS_wd, KNZ_wd)

# Calculate the number of dry then wet years
CDR_dw <- sum(spei_lump_CDR$spei_category[-nrow(spei_lump_CDR)] == "Dry" & spei_lump_CDR$spei_category[-1] == "Wet")
KBS_dw <- sum(spei_lump_KBS$spei_category[-nrow(spei_lump_KBS)] == "Dry" & spei_lump_KBS$spei_category[-1] == "Wet")
KNZ_dw <- sum(spei_lump_KNZ$spei_category[-nrow(spei_lump_KNZ)] == "Dry" & spei_lump_KNZ$spei_category[-1] == "Wet")

dry_wet <- c(CDR_dw, KBS_dw, KNZ_dw)

# Calculate the number of wet then wet years
CDR_ww <- sum(spei_lump_CDR$spei_category[-nrow(spei_lump_CDR)] == "Wet" & spei_lump_CDR$spei_category[-1] == "Wet")
KBS_ww <- sum(spei_lump_KBS$spei_category[-nrow(spei_lump_KBS)] == "Wet" & spei_lump_KBS$spei_category[-1] == "Wet")
KNZ_ww <- sum(spei_lump_KNZ$spei_category[-nrow(spei_lump_KNZ)] == "Wet" & spei_lump_KNZ$spei_category[-1] == "Wet")

wet_wet <- c(CDR_ww, KBS_ww, KNZ_ww)

# Calculate the number of dry then dry years
CDR_dd <- sum(spei_lump_CDR$spei_category[-nrow(spei_lump_CDR)] == "Dry" & spei_lump_CDR$spei_category[-1] == "Dry")
KBS_dd <- sum(spei_lump_KBS$spei_category[-nrow(spei_lump_KBS)] == "Dry" & spei_lump_KBS$spei_category[-1] == "Dry")
KNZ_dd <- sum(spei_lump_KNZ$spei_category[-nrow(spei_lump_KNZ)] == "Dry" & spei_lump_KNZ$spei_category[-1] == "Dry")

dry_dry <- c(CDR_dd, KBS_dd, KNZ_dd)

spei_table_lump <- data.frame(sites, wet, dry, wet_dry, dry_wet, wet_wet, dry_dry)

# Look at plant biomass averages by event category
plot_spei %>%
  ggplot(aes(x = spei_category, y = plot_biomass)) + 
  geom_sina() +
  stat_summary(fun.data = mean_cl_boot, color = "red")

# Are there any significant differences?
biomass.lm <- lmer(plot_biomass ~ spei_category + (1|site) + (1|site:plot), data = plot_spei)
anova(biomass.lm)
emmeans(biomass.lm, list(pairwise ~ spei_category), adjust = "tukey")
plot_model(
  biomass.lm,
  type = "pred",
  terms = c("spei_category")
)

# Look at plant biomass vs spei
plot_spei %>%
  ggplot(aes(x = spei12, y = plot_biomass)) + 
  geom_point()

# Is there a quadratic relationship?
biomass_spei12.lm <- lmer(plot_biomass ~ spei12 + I(spei12^2) + (1|site) + (1|site:plot), data = plot_spei)
summary(biomass_spei12.lm) # quadratic term significant

plot_model(
  biomass_spei12.lm,
  type = "pred",
  terms="spei12[all]",
  show.data = TRUE
)

# Why 3663 missing values for plot biomass
sum(is.na(plot_spei$plot_biomass))

  