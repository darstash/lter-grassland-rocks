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
library(bbmle)
library(MuMIn)
library(DHARMa)
library(car)

# Set working directory 
L0_dir <- Sys.getenv("L0DIR")
L1_dir <- Sys.getenv("L1DIR")
L2_dir <- Sys.getenv("L2DIR")
list.files(L2_dir)

# Read in CSV files ----
plot_spei <- read.csv(file.path(L2_dir, "plot_metrics_SPEI.csv"))
metadata <- read.csv(file.path(L2_dir, "metadata.csv"))

# Make sure spei category is a factor
str(plot_spei)
plot_spei$spei12_category <- factor(plot_spei$spei12_category, levels = c("Extreme dry", "Moderate dry", "Normal", "Moderate wet", "Extreme wet"))

# Summarize spei by year
spei_summary <- plot_spei %>%
  group_by(year, spei12_category, site) %>%
  summarize(spei12 = mean(spei12),
            spei9 = mean(spei9),
            spei6 = mean(spei6),
            spei3 = mean(spei3))

# Merge metadata with spei to get treatment column
metadata_filter <- metadata %>%
  drop_na(uniqueid)
plot_filter <- plot_spei %>%
  drop_na(uniqueid) %>%
  filter(higher_order_organization != "RaMPs_Block1",
         higher_order_organization != "RaMPs_Block2",
         higher_order_organization != "RaMPs_Block3")
x <- plot_filter %>% group_by(site, year, uniqueid, higher_order_organization) %>% summarize(count = n())
plot_trt <- left_join(plot_filter, metadata_filter)

# Look at SPEI plots for each site in the range the data exist
spei_summary %>%
  filter(site == "CDR") %>%
  ggplot(aes (x = year, y = spei12)) +
  geom_line(alpha = 0.5) + 
  geom_point( aes(color = spei12_category), size = 1.25) + theme_bw() +
  scale_x_continuous(breaks = seq(1980, 2025, by = 10)) +
  annotate( "text", label = "CDR",
            x = 2000, y = 2, size = 5, colour = "black") + 
  scale_color_manual(values = c("#F5191CFF", "#E78200FF", "#EAC728FF", "#81BB95FF", "#3B99B1FF"))

spei_summary %>%
  filter(site == "KBS") %>%
  ggplot(aes (x = year, y = spei12)) +
  geom_line(alpha = 0.5) + 
  geom_point( aes(color = spei12_category), size = 1.25) + theme_bw() +
  scale_x_continuous(breaks = seq(1990, 2025, by = 10)) +
  annotate( "text", label = "KBS",
            x = 2000, y = 2, size = 5, colour = "black") + 
  scale_color_manual(values = c("#E78200FF", "#EAC728FF", "#81BB95FF", "#3B99B1FF"))


spei_summary %>%
  filter(site == "KNZ") %>%
  ggplot(aes(x = year, y = spei12)) +
  geom_line(alpha = 0.5) + 
  geom_point( aes(color = spei12_category), size = 1.25) + theme_bw() +
  scale_x_continuous(breaks = seq(1990, 2025, by = 10)) +
  annotate( "text", label = "KNZ",
            x = 2000, y = 2, size = 5, colour = "black") + 
  scale_color_manual(values = c("#F5191CFF", "#E78200FF", "#EAC728FF", "#81BB95FF", "#3B99B1FF"))

# Calculate the number of extreme wet years
CDR_ew <- length(which(spei_summary$site == "CDR" & spei_summary$spei12_category == "Extreme wet"))
KBS_ew <- length(which(spei_summary$site == "KBS" & spei_summary$spei12_category == "Extreme wet"))
KNZ_ew <- length(which(spei_summary$site == "KNZ" & spei_summary$spei12_category == "Extreme wet"))
extreme_wet <- c(CDR_ew, KBS_ew, KNZ_ew)
sites <- c("CDR", "KBS", "KNZ")

# Calculate the number of moderate wet years
CDR_mw <- length(which(spei_summary$site == "CDR" & spei_summary$spei12_category == "Moderate wet"))
KBS_mw <- length(which(spei_summary$site == "KBS" & spei_summary$spei12_category == "Moderate wet"))
KNZ_mw <- length(which(spei_summary$site == "KNZ" & spei_summary$spei12_category == "Moderate wet"))
moderate_wet <- c(CDR_mw, KBS_mw, KNZ_mw)

# Calculate the number of extreme dry years
CDR_ed <- length(which(spei_summary$site == "CDR" & spei_summary$spei12_category == "Extreme dry"))
KBS_ed <- length(which(spei_summary$site == "KBS" & spei_summary$spei12_category == "Extreme dry"))
KNZ_ed <- length(which(spei_summary$site == "KNZ" & spei_summary$spei12_category == "Extreme dry"))
extreme_dry <- c(CDR_ed, KBS_ed, KNZ_ed)

# Calculate the number of moderate dry years
CDR_md <- length(which(spei_summary$site == "CDR" & spei_summary$spei12_category == "Moderate dry"))
KBS_md <- length(which(spei_summary$site == "KBS" & spei_summary$spei12_category == "Moderate dry"))
KNZ_md <- length(which(spei_summary$site == "KNZ" & spei_summary$spei12_category == "Moderate dry"))
moderate_dry <- c(CDR_md, KBS_md, KNZ_md)

# Calculate the number of extreme wet then extreme dry years
spei_summary_CDR <- spei_summary %>%
  filter(site == "CDR")
CDR_ewd <- sum(spei_summary_CDR$spei12_category[-nrow(spei_summary_CDR)] == "Extreme wet" & spei_summary_CDR$spei12_category[-1] == "Extreme dry") # This was written by ChatGPT
spei_summary_KBS <- spei_summary %>%
  filter(site == "KBS")
KBS_ewd <- sum(spei_summary_KBS$spei12_category[-nrow(spei_summary_KBS)] == "Extreme wet" & spei_summary_KBS$spei12_category[-1] == "Extreme dry")
spei_summary_KNZ <- spei_summary %>%
  filter(site == "KNZ")
KNZ_ewd <- sum(spei_summary_KNZ$spei12_category[-nrow(spei_summary_KNZ)] == "Extreme wet" & spei_summary_KNZ$spei12_category[-1] == "Extreme dry")

extreme_wet_dry <- c(CDR_ewd, KBS_ewd, KNZ_ewd)

# Calculate the number of extreme dry then extreme wet years
CDR_edw <- sum(spei_summary_CDR$spei12_category[-nrow(spei_summary_CDR)] == "Extreme dry" & spei_summary_CDR$spei12_category[-1] == "Extreme wet")
KBS_edw <- sum(spei_summary_KBS$spei12_category[-nrow(spei_summary_KBS)] == "Extreme dry" & spei_summary_KBS$spei12_category[-1] == "Extreme wet")
KNZ_edw <- sum(spei_summary_KNZ$spei12_category[-nrow(spei_summary_KNZ)] == "Extreme dry" & spei_summary_KNZ$spei12_category[-1] == "Extreme wet")

extreme_dry_wet <- c(CDR_edw, KBS_edw, KNZ_edw)

# Calculate the number of extreme wet then extreme wet years
CDR_eww <- sum(spei_summary_CDR$spei12_category[-nrow(spei_summary_CDR)] == "Extreme wet" & spei_summary_CDR$spei12_category[-1] == "Extreme wet")
KBS_eww <- sum(spei_summary_KBS$spei12_category[-nrow(spei_summary_KBS)] == "Extreme wet" & spei_summary_KBS$spei12_category[-1] == "Extreme wet")
KNZ_eww <- sum(spei_summary_KNZ$spei12_category[-nrow(spei_summary_KNZ)] == "Extreme wet" & spei_summary_KNZ$spei12_category[-1] == "Extreme wet")

extreme_wet_wet <- c(CDR_eww, KBS_eww, KNZ_eww)

# Calculate the number of extreme dry then extreme dry years
CDR_edd <- sum(spei_summary_CDR$spei12_category[-nrow(spei_summary_CDR)] == "Extreme dry" & spei_summary_CDR$spei12_category[-1] == "Extreme dry")
KBS_edd <- sum(spei_summary_KBS$spei12_category[-nrow(spei_summary_KBS)] == "Extreme dry" & spei_summary_KBS$spei12_category[-1] == "Extreme dry")
KNZ_edd <- sum(spei_summary_KNZ$spei12_category[-nrow(spei_summary_KNZ)] == "Extreme dry" & spei_summary_KNZ$spei12_category[-1] == "Extreme dry")

extreme_dry_dry <- c(CDR_edd, KBS_edd, KNZ_edd)

spei_table <- data.frame(sites, extreme_wet, moderate_wet, extreme_dry, moderate_dry, extreme_wet_dry, extreme_dry_wet, extreme_wet_wet, extreme_dry_dry)

# Not enough extreme events. Let's lump all moderate and extreme events into just wet and dry categories
spei_lump <- spei_summary %>%
  mutate(spei12_category = ifelse(as.character(spei12_category) == "Extreme dry", "Dry", as.character(spei12_category))) %>%
  mutate(spei12_category = ifelse(as.character(spei12_category) == "Moderate dry", "Dry", as.character(spei12_category))) %>%
  mutate(spei12_category = ifelse(as.character(spei12_category) == "Extreme wet", "Wet", as.character(spei12_category))) %>%
  mutate(spei12_category = ifelse(as.character(spei12_category) == "Moderate wet", "Wet", as.character(spei12_category)))

# Calculate the number of moderate/extreme wet years
CDR_w <- length(which(spei_lump$site == "CDR" & spei_lump$spei12_category == "Wet"))
KBS_w <- length(which(spei_lump$site == "KBS" & spei_lump$spei12_category == "Wet"))
KNZ_w <- length(which(spei_lump$site == "KNZ" & spei_lump$spei12_category == "Wet"))
wet <- c(CDR_w, KBS_w, KNZ_w)

# Calculate the number of moderate/extreme dry years
CDR_d <- length(which(spei_lump$site == "CDR" & spei_lump$spei12_category == "Dry"))
KBS_d <- length(which(spei_lump$site == "KBS" & spei_lump$spei12_category == "Dry"))
KNZ_d <- length(which(spei_lump$site == "KNZ" & spei_lump$spei12_category == "Dry"))
dry <- c(CDR_d, KBS_d, KNZ_d)

# Calculate the number of wet then dry years
spei_lump_CDR <- spei_lump %>%
  filter(site == "CDR")
CDR_wd <- sum(spei_lump_CDR$spei12_category[-nrow(spei_lump_CDR)] == "Wet" & spei_lump_CDR$spei12_category[-1] == "Dry")
spei_lump_KBS <- spei_lump %>%
  filter(site == "KBS")
KBS_wd <- sum(spei_lump_KBS$spei12_category[-nrow(spei_lump_KBS)] == "Wet" & spei_lump_KBS$spei12_category[-1] == "Dry")
spei_lump_KNZ <- spei_lump %>%
  filter(site == "KNZ")
KNZ_wd <- sum(spei_lump_KNZ$spei12_category[-nrow(spei_lump_KNZ)] == "Wet" & spei_lump_KNZ$spei12_category[-1] == "Dry")

wet_dry <- c(CDR_wd, KBS_wd, KNZ_wd)

# Calculate the number of wet then dry years
spei_lump_CDR <- spei_lump %>%
  filter(site == "CDR")
CDR_wd <- sum(spei_lump_CDR$spei12_category[-nrow(spei_lump_CDR)] == "Wet" & spei_lump_CDR$spei12_category[-1] == "Dry")
spei_lump_KBS <- spei_lump %>%
  filter(site == "KBS")
KBS_wd <- sum(spei_lump_KBS$spei12_category[-nrow(spei_lump_KBS)] == "Wet" & spei_lump_KBS$spei12_category[-1] == "Dry")
spei_lump_KNZ <- spei_lump %>%
  filter(site == "KNZ")
KNZ_wd <- sum(spei_lump_KNZ$spei12_category[-nrow(spei_lump_KNZ)] == "Wet" & spei_lump_KNZ$spei12_category[-1] == "Dry")

wet_dry <- c(CDR_wd, KBS_wd, KNZ_wd)

# Calculate the number of dry then wet years
CDR_dw <- sum(spei_lump_CDR$spei12_category[-nrow(spei_lump_CDR)] == "Dry" & spei_lump_CDR$spei12_category[-1] == "Wet")
KBS_dw <- sum(spei_lump_KBS$spei12_category[-nrow(spei_lump_KBS)] == "Dry" & spei_lump_KBS$spei12_category[-1] == "Wet")
KNZ_dw <- sum(spei_lump_KNZ$spei12_category[-nrow(spei_lump_KNZ)] == "Dry" & spei_lump_KNZ$spei12_category[-1] == "Wet")

dry_wet <- c(CDR_dw, KBS_dw, KNZ_dw)

# Calculate the number of wet then wet years
CDR_ww <- sum(spei_lump_CDR$spei12_category[-nrow(spei_lump_CDR)] == "Wet" & spei_lump_CDR$spei12_category[-1] == "Wet")
KBS_ww <- sum(spei_lump_KBS$spei12_category[-nrow(spei_lump_KBS)] == "Wet" & spei_lump_KBS$spei12_category[-1] == "Wet")
KNZ_ww <- sum(spei_lump_KNZ$spei12_category[-nrow(spei_lump_KNZ)] == "Wet" & spei_lump_KNZ$spei12_category[-1] == "Wet")

wet_wet <- c(CDR_ww, KBS_ww, KNZ_ww)

# Calculate the number of dry then dry years
CDR_dd <- sum(spei_lump_CDR$spei12_category[-nrow(spei_lump_CDR)] == "Dry" & spei_lump_CDR$spei12_category[-1] == "Dry")
KBS_dd <- sum(spei_lump_KBS$spei12_category[-nrow(spei_lump_KBS)] == "Dry" & spei_lump_KBS$spei12_category[-1] == "Dry")
KNZ_dd <- sum(spei_lump_KNZ$spei12_category[-nrow(spei_lump_KNZ)] == "Dry" & spei_lump_KNZ$spei12_category[-1] == "Dry")

dry_dry <- c(CDR_dd, KBS_dd, KNZ_dd)

spei_table_lump <- data.frame(sites, wet, dry, wet_dry, dry_wet, wet_wet, dry_dry)

# Look at plant biomass averages by event category
plot_spei %>%
  ggplot(aes(x = spei12_category, y = plot_biomass)) + 
  geom_sina() +
  stat_summary(fun.data = mean_cl_boot, color = "red")

# Are there any significant differences?
biomass.lm <- lmer(plot_biomass ~ spei12_category + (1|site) + (1|site:plot), data = plot_spei)
anova(biomass.lm)
emmeans(biomass.lm, list(pairwise ~ spei12_category), adjust = "tukey")
plot_model(
  biomass.lm,
  type = "pred",
  terms = c("spei12_category")
)

# Look at plant biomass vs spei12
plot_spei %>%
  ggplot(aes(x = spei12, y = plot_biomass)) + 
  geom_point()

plot_spei %>%
  ggplot(aes(x = spei12, y = plot_biomass)) + 
  geom_point() +
  facet_wrap(~site)

# Is there a quadratic relationship?
biomass_spei12.lm <- lmer(plot_biomass ~ spei12 + I(spei12^2) + (1|site) + (1|site:plot), data = plot_spei)
summary(biomass_spei12.lm) # quadratic term significant
biomass_spei12.lm2 <- lmer(plot_biomass ~ spei12 + (1|site) + (1|site:plot), data = plot_spei)


plot_model(
  biomass_spei12.lm,
  type = "pred",
  terms="spei12[all]",
  show.data = TRUE
)

# Why 1 missing values for plot biomass
sum(is.na(plot_spei$plot_biomass))


# Compare spei3, 6, 9, 12 models
# They all appear to have potential quadratic terms # Quadratic makes sense conceptually
plot_spei %>%
  ggplot(aes(x = spei3, y = plot_biomass)) + 
  geom_point()

plot_spei %>%
  ggplot(aes(x = spei6, y = plot_biomass)) + 
  geom_point()

plot_spei %>%
  ggplot(aes(x = spei9, y = plot_biomass)) + 
  geom_point()

biomass_spei9.lm <- lmer(plot_biomass ~ spei9 + I(spei9^2) + (1|site) + (1|site:plot), data = plot_spei)
summary(biomass_spei9.lm) # Quadratic term significant
biomass_spei9.lm2 <- lmer(plot_biomass ~ spei9 + (1|site) + (1|site:plot), data = plot_spei)


biomass_spei6.lm <- lmer(plot_biomass ~ spei6 + I(spei6^2) + (1|site) + (1|site:plot), data = plot_spei)
summary(biomass_spei6.lm)  # Quadratic term significant
biomass_spei6.lm2 <- lmer(plot_biomass ~ spei6 + (1|site) + (1|site:plot), data = plot_spei)


biomass_spei3.lm <- lmer(plot_biomass ~ spei3 + I(spei3^2) + (1|site) + (1|site:plot), data = plot_spei)
summary(biomass_spei3.lm) # Quadratic term significant (is this just because of sample size?)
biomass_spei3.lm2 <- lmer(plot_biomass ~ spei3 + (1|site) + (1|site:plot), data = plot_spei)


# Compare models using AICc like Robinson paper
AICctab(biomass_spei3.lm, biomass_spei3.lm2, biomass_spei6.lm, biomass_spei6.lm2, biomass_spei9.lm, biomass_spei9.lm2, biomass_spei12.lm, biomass_spei12.lm2) # spei9 lowest

# Compare model R2 like Robinson
r.squaredGLMM(biomass_spei3.lm)
r.squaredGLMM(biomass_spei3.lm2)
r.squaredGLMM(biomass_spei6.lm)
r.squaredGLMM(biomass_spei6.lm2)
r.squaredGLMM(biomass_spei9.lm)
r.squaredGLMM(biomass_spei9.lm2)
r.squaredGLMM(biomass_spei12.lm)
r.squaredGLMM(biomass_spei12.lm2)

# Plot spei 9 moel
plot_model(
  biomass_spei9.lm,
  type = "pred",
  terms="spei9[all]",
  show.data = TRUE
)

# Look at spei 9 model diagnostic plots
simulationOutput <- simulateResiduals(fittedModel = biomass_spei9.lm, plot = F)
plot(simulationOutput) # something going on here, but could be due to sample size

# Comparing spei by site
# CDR
biomass_cdr_spei12.lm <- lmer(plot_biomass ~ spei12 + I(spei12^2) + (1|plot), data = plot_spei, subset = site == "CDR")
biomass_cdr_spei9.lm <- lmer(plot_biomass ~ spei9 + I(spei9^2) + (1|plot), data = plot_spei, subset = site == "CDR")
biomass_cdr_spei6.lm <- lmer(plot_biomass ~ spei6 + I(spei6^2) + (1|plot), data = plot_spei, subset = site == "CDR")
biomass_cdr_spei3.lm <- lmer(plot_biomass ~ spei3 + I(spei3^2) + (1|plot), data = plot_spei, subset = site == "CDR")

AICctab(biomass_cdr_spei12.lm, biomass_cdr_spei9.lm, biomass_cdr_spei6.lm, biomass_cdr_spei3.lm) # spei6 slightly better than 9

# KBS
biomass_kbs_spei12.lm <- lmer(plot_biomass ~ spei12 + I(spei12^2) + (1|plot), data = plot_spei, subset = site == "KBS")
biomass_kbs_spei9.lm <- lmer(plot_biomass ~ spei9 + I(spei9^2) + (1|plot), data = plot_spei, subset = site == "KBS")
biomass_kbs_spei6.lm <- lmer(plot_biomass ~ spei6 + I(spei6^2) + (1|plot), data = plot_spei, subset = site == "KBS")
biomass_kbs_spei3.lm <- lmer(plot_biomass ~ spei3 + I(spei3^2) + (1|plot), data = plot_spei, subset = site == "KBS")

AICctab(biomass_kbs_spei12.lm, biomass_kbs_spei9.lm, biomass_kbs_spei6.lm, biomass_kbs_spei3.lm) # spei9 better

# KNZ
biomass_knz_spei12.lm <- lmer(plot_biomass ~ spei12 + I(spei12^2) + (1|plot), data = plot_spei, subset = site == "KNZ")
biomass_knz_spei9.lm <- lmer(plot_biomass ~ spei9 + I(spei9^2) + (1|plot), data = plot_spei, subset = site == "KNZ")
biomass_knz_spei6.lm <- lmer(plot_biomass ~ spei6 + I(spei6^2) + (1|plot), data = plot_spei, subset = site == "KNZ")
biomass_knz_spei3.lm <- lmer(plot_biomass ~ spei3 + I(spei3^2) + (1|plot), data = plot_spei, subset = site == "KNZ")

AICctab(biomass_knz_spei12.lm, biomass_knz_spei9.lm, biomass_knz_spei6.lm, biomass_knz_spei3.lm) # spei9 better



# Only look at control biomass ----
plot_control <- plot_trt %>%
  filter(treatment == "control")

# Look at plant biomass averages by event category
plot_control %>%
  ggplot(aes(x = spei12_category, y = plot_biomass)) + 
  geom_sina(alpha = 0.2) +
  stat_summary(fun.data = mean_cl_boot, color = "red")

# Look at plant biomass vs spei12
plot_control %>%
  ggplot(aes(x = spei12, y = plot_biomass)) + 
  geom_point(alpha=0.2) +
  geom_smooth() +
  ggpubr::stat_cor()

plot_control %>%
  ggplot(aes(x = spei12, y = plot_biomass)) + 
  geom_point(alpha  = 0.2) +
  facet_wrap(~site)

# Model comparison with spei12, 9, 6, 3 for ONLY control plots
control_spei3.lm2 <- lmer(plot_biomass ~ spei3 + I(spei3^2) + (1|site) + (1|site:plot), data = plot_control)
summary(control_spei3.lm2) # Quadratic term significant
control_spei3.lm <- lmer(plot_biomass ~ spei3 + (1|site) + (1|site:plot), data = plot_control)
simres <- simulateResiduals(control_spei3.lm2)
plot(simres)

control_spei6.lm2 <- lmer(plot_biomass ~ spei6 + I(spei6^2) + (1|site) + (1|site:plot), data = plot_control)
summary(control_spei6.lm2) # Quadratic term significant
control_spei6.lm <- lmer(plot_biomass ~ spei6 + (1|site) + (1|site:plot), data = plot_control)

control_spei9.lm2 <- lmer(plot_biomass ~ spei9 + I(spei9^2) + (1|site) + (1|site:plot), data = plot_control)
summary(control_spei9.lm2) # Quadratic term significant
control_spei9.lm <- lmer(plot_biomass ~ spei9 + (1|site) + (1|site:plot), data = plot_control)

control_spei12.lm2 <- lmer(plot_biomass ~ spei12 + I(spei12^2) + (1|site) + (1|site:plot), data = plot_control)
summary(control_spei12.lm2) # Quadratic term significant
control_spei12.lm <- lmer(plot_biomass ~ spei12 + (1|site) + (1|site:plot), data = plot_control)

AICctab(control_spei3.lm2, control_spei3.lm, control_spei6.lm2, control_spei6.lm, control_spei9.lm2, control_spei9.lm, control_spei12.lm2, control_spei12.lm)

# SPEI 6 with a quadratic term best model when all sites in one model

# Compare model R2 like Robinson
r.squaredGLMM(control_spei3.lm2)
r.squaredGLMM(control_spei3.lm)
r.squaredGLMM(control_spei6.lm2) # best R2m
r.squaredGLMM(control_spei6.lm)
r.squaredGLMM(control_spei9.lm2)
r.squaredGLMM(control_spei9.lm)
r.squaredGLMM(control_spei12.lm2)
r.squaredGLMM(control_spei12.lm)

# Plot the quadratic SPEI6 model
plot_model(
  control_spei6.lm2,
  type = "pred",
  terms="spei6[all]",
  show.data = TRUE
)

# Add experiment column
plot_control$experiment <- sub("nutnet.*", "nutnet", plot_control$higher_order_organization)
plot_control$experiment <- sub("glbrc_scaleup.*", "glbrc_scaleup", plot_control$experiment)
plot_control$experiment <- sub("glbrc_G10.*", "glbrc", plot_control$experiment)
plot_control$experiment <- sub("glbrc_G9.*", "glbrc", plot_control$experiment)
plot_control$experiment <- sub("mcse.*", "mcse", plot_control$experiment)
plot_control$experiment <- sub("microplots.*", "microplots", plot_control$experiment)
plot_control$experiment <- sub("Experiment 1.*", "Experiment 1", plot_control$experiment)
plot_control$experiment <- sub("001d_A_fl", "001d_fl", plot_control$experiment)
plot_control$experiment <- sub("001d_B_fl", "001d_fl", plot_control$experiment)
plot_control$experiment <- sub("001d_C_fl", "001d_fl", plot_control$experiment)
plot_control$experiment <- sub("001d_D_fl", "001d_fl", plot_control$experiment)
plot_control$experiment <- sub("001d_A_tu", "001d_tu", plot_control$experiment)
plot_control$experiment <- sub("001d_B_tu", "001d_tu", plot_control$experiment)
plot_control$experiment <- sub("001d_C_tu", "001d_tu", plot_control$experiment)
plot_control$experiment <- sub("001d_D_tu", "001d_tu", plot_control$experiment)
plot_control$experiment <- sub("004a_A_fl", "004a_fl", plot_control$experiment)
plot_control$experiment <- sub("004a_B_fl", "004a_fl", plot_control$experiment)
plot_control$experiment <- sub("004a_C_fl", "004a_fl", plot_control$experiment)
plot_control$experiment <- sub("004a_D_fl", "004a_fl", plot_control$experiment)
plot_control$experiment <- sub("004a_A_tu", "004a_tu", plot_control$experiment)
plot_control$experiment <- sub("004a_B_tu", "004a_tu", plot_control$experiment)
plot_control$experiment <- sub("004a_C_tu", "004a_tu", plot_control$experiment)
plot_control$experiment <- sub("004a_D_tu", "004a_tu", plot_control$experiment)
plot_control$experiment <- sub("004b_A_fl", "004b_fl", plot_control$experiment)
plot_control$experiment <- sub("004b_B_fl", "004b_fl", plot_control$experiment)
plot_control$experiment <- sub("004b_C_fl", "004b_fl", plot_control$experiment)
plot_control$experiment <- sub("004b_D_fl", "004b_fl", plot_control$experiment)
plot_control$experiment <- sub("004b_A_tu", "004b_tu", plot_control$experiment)
plot_control$experiment <- sub("004b_B_tu", "004b_tu", plot_control$experiment)
plot_control$experiment <- sub("004b_C_tu", "004b_tu", plot_control$experiment)
plot_control$experiment <- sub("004b_D_tu", "004b_tu", plot_control$experiment)
plot_control$experiment <- sub("002d.*", "002d", plot_control$experiment)
plot_control$experiment <- sub("Experiment 54.*", "Experiment 54", plot_control$experiment)
plot_control$experiment <- sub("KNZ_WAT01.*", "KNZ_WAT01", plot_control$experiment)
plot_control$experiment <- sub("002c.*", "002c", plot_control$experiment)
plot_control$experiment <- sub("e061.*", "e061", plot_control$experiment)
plot_control$experiment <- sub("e247.*", "e247", plot_control$experiment)
plot_control$experiment <- sub("e245.*", "e245", plot_control$experiment)
plot_control$experiment[plot_control$experiment == "A"] <- "NGE"
plot_control$experiment[plot_control$experiment == "B"] <- "NGE"
plot_control$experiment[plot_control$experiment == "C"] <- "NGE"
plot_control$experiment[plot_control$experiment == "D"] <- "NGE"
plot_control$experiment[plot_control$experiment == "E"] <- "NGE"
plot_control$experiment[plot_control$experiment == "F"] <- "NGE"

# Make year factor for random effect
plot_control$year <- as.factor(plot_control$year)

# Try logging biomass and adding more random effects
# Model comparison with spei12, 9, 6, 3 for ONLY control plots
control_spei3.lm2.ln <- lmer(log1p(plot_biomass) ~ spei3 + I(spei3^2) + (1|site/experiment/uniqueid) + (1|year), data = plot_control)
summary(control_spei3.lm2.ln) # Quadratic term significant
control_spei3.lm.ln <- lmer(log1p(plot_biomass) ~ spei3 + (1|site/experiment/uniqueid) + (1|year), data = plot_control)
simres <- simulateResiduals(control_spei3.lm2.ln)
plot(simres) # Better 

control_spei6.lm2.ln <- lmer(log1p(plot_biomass) ~ spei6 + I(spei6^2) + (1|site/experiment/uniqueid)
                             + (1|year), data = plot_control)
summary(control_spei6.lm2.ln) # Quadratic term not significant
control_spei6.lm.ln <- lmer(log1p(plot_biomass) ~ spei6 + (1|site/experiment/uniqueid)
                            + (1|year), data = plot_control)
simres <- simulateResiduals(control_spei6.lm2.ln)
plot(simres) # Better 

control_spei9.lm2.ln <- lmer(log1p(plot_biomass) ~ spei9 + I(spei9^2) + (1|site/experiment/uniqueid) + (1|year), data = plot_control)
summary(control_spei9.lm2.ln) # Quadratic term marginally significant
control_spei9.lm.ln <- lmer(log1p(plot_biomass) ~ spei9 + (1|site/experiment/uniqueid) + (1|year), data = plot_control)
summary(control_spei9.lm.ln)
simres <- simulateResiduals(control_spei9.lm2.ln)
plot(simres) # Better 
simres <- simulateResiduals(control_spei9.lm.ln)
plot(simres) # Okay 

control_spei12.lm2.ln <- lmer(log1p(plot_biomass) ~ spei12 + I(spei12^2) + (1|site/experiment/uniqueid) + (1|year), data = plot_control, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) # Failed to converge with defaults
summary(control_spei12.lm2.ln)
control_spei12.lm.ln <- lmer(log1p(plot_biomass) ~ spei12 + (1|site/experiment/uniqueid) + (1|year), data = plot_control, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) # Failed to converge with defaults
simres <- simulateResiduals(control_spei12.lm2.ln)
plot(simres) # Better 

AICctab(control_spei3.lm2.ln, control_spei3.lm.ln, control_spei6.lm2.ln, control_spei6.lm.ln, control_spei9.lm2.ln, control_spei9.lm.ln, control_spei12.lm2.ln, control_spei12.lm.ln)

# Compare model R2 like Robinson
r.squaredGLMM(control_spei3.lm2.ln)
r.squaredGLMM(control_spei3.lm.ln)
r.squaredGLMM(control_spei6.lm2.ln)
r.squaredGLMM(control_spei6.lm.ln)
r.squaredGLMM(control_spei9.lm2.ln) # best R2m
r.squaredGLMM(control_spei9.lm.ln)  # second best R2m
r.squaredGLMM(control_spei12.lm2.ln)
r.squaredGLMM(control_spei12.lm.ln)

# Plot the SPEI9 model
plot_model(
  control_spei9.lm.ln,
  type = "pred",
  terms="spei9",
  show.data = TRUE
)

# Log response of SPEI ----
# Using SPEI3
plot_control_lrr3 <- plot_control %>%
  group_by(uniqueid) %>%
  arrange(uniqueid, year) %>%
  mutate(prior_year_biomass = lag(plot_biomass),
         prior_year_type = lag(spei3_category)) %>%
  ungroup()
normal_biomass3 <- plot_control_lrr3  %>%
  group_by(uniqueid) %>%
  filter(spei3_category == "Normal" | spei3_category == "Moderate wet" | spei3_category == "Moderate dry") %>%
  summarize(average_normal_biomass = mean(plot_biomass, na.rm=T))
plot_control_lrr3 <- left_join(plot_control_lrr3, normal_biomass3)

lrr3 <- plot_control_lrr3 %>%
  filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
  filter(spei3_category == "Extreme wet" | spei3_category == "Extreme dry") %>%
  mutate(LRR = log(plot_biomass / average_normal_biomass))  # Calculate the log response ratio

# Take absolute value of SPEI3
lrr3 <- lrr3 %>%
  mutate(abs_spei3 = abs(spei3))

lrr.lm3 <- lmer(LRR ~ abs_spei3*spei3_category + (1|site/experiment/uniqueid) + (1|year), data = lrr3)
summary(lrr.lm3)
simres <- simulateResiduals(lrr.lm3)
plot(simres) # Weird

# Plot the LRR model
plot_model(
  lrr.lm3,
  type = "pred",
  terms= c("abs_spei3", "spei3_category"),
  show.data = TRUE
)

# Using SPEI6
plot_control_lrr <- plot_control %>%
  group_by(uniqueid) %>%
  arrange(uniqueid, year) %>%
  mutate(prior_year_biomass = lag(plot_biomass),
         prior_year_type = lag(spei6_category)) %>%
  ungroup()
normal_biomass <- plot_control_lrr  %>%
  group_by(uniqueid) %>%
  filter(spei6_category == "Normal" | spei6_category == "Moderate wet" | spei6_category == "Moderate dry") %>%
  summarize(average_normal_biomass = mean(plot_biomass, na.rm=T))
plot_control_lrr <- left_join(plot_control_lrr, normal_biomass)

lrr <- plot_control_lrr %>%
  filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
  filter(spei6_category == "Extreme wet" | spei6_category == "Extreme dry") %>%
  mutate(LRR = log(plot_biomass / average_normal_biomass))  # Calculate the log response ratio

# Take absolute value of SPEI6
lrr <- lrr %>%
  mutate(abs_spei6 = abs(spei6))

lrr.lm <- lmer(LRR ~ abs_spei6*spei6_category + (1|site/experiment/uniqueid) + (1|year), data = lrr)
summary(lrr.lm)
simres <- simulateResiduals(lrr.lm)
plot(simres) # Weird

# Plot the LRR model
plot_model(
  lrr.lm,
  type = "pred",
  terms= c("abs_spei6", "spei6_category"),
  show.data = TRUE
)

# Using SPEI9
plot_control_lrr9 <- plot_control %>%
  group_by(uniqueid) %>%
  arrange(uniqueid, year) %>%
  mutate(prior_year_biomass = lag(plot_biomass),
         prior_year_type = lag(spei9_category)) %>%
  ungroup()
normal_biomass9 <- plot_control_lrr9  %>%
  group_by(uniqueid) %>%
  filter(spei9_category == "Normal" | spei9_category == "Moderate wet" | spei9_category == "Moderate dry") %>%
  summarize(average_normal_biomass = mean(plot_biomass, na.rm=T))
plot_control_lrr9 <- left_join(plot_control_lrr9, normal_biomass9)

lrr9 <- plot_control_lrr9 %>%
  filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
  filter(spei9_category == "Extreme wet" | spei9_category == "Extreme dry") %>%
  mutate(LRR = log(plot_biomass / average_normal_biomass))  # Calculate the log response ratio

# Take absolute value of SPEI9
lrr9 <- lrr9 %>%
  mutate(abs_spei9 = abs(spei9))

lrr.lm9 <- lmer(LRR ~ abs_spei9*spei9_category + (1|site/experiment/uniqueid) + (1|year), data = lrr9)
summary(lrr.lm9)
simres <- simulateResiduals(lrr.lm9)
plot(simres)

# Plot the LRR model
plot_model(
  lrr.lm9,
  type = "pred",
  terms= c("abs_spei9", "spei9_category"),
  show.data = TRUE
)

# Using SPEI12
plot_control_lrr12 <- plot_control %>%
  group_by(uniqueid) %>%
  arrange(uniqueid, year) %>%
  mutate(prior_year_biomass = lag(plot_biomass),
         prior_year_type = lag(spei12_category)) %>%
  ungroup()
normal_biomass12 <- plot_control_lrr12  %>%
  group_by(uniqueid) %>%
  filter(spei12_category == "Normal" | spei12_category == "Moderate wet" | spei12_category == "Moderate dry") %>%
  summarize(average_normal_biomass = mean(plot_biomass, na.rm=T))
plot_control_lrr12 <- left_join(plot_control_lrr12, normal_biomass12)

lrr12 <- plot_control_lrr12 %>%
  filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
  filter(spei12_category == "Extreme wet" | spei12_category == "Extreme dry") %>%
  mutate(LRR = log(plot_biomass / average_normal_biomass))  # Calculate the log response ratio

# Take absolute value of SPEI6
lrr12 <- lrr12 %>%
  mutate(abs_spei12 = abs(spei12))

lrr.lm12 <- lmer(LRR ~ abs_spei12*spei12_category + (1|site/experiment/uniqueid) + (1|year), data = lrr12)
summary(lrr.lm12)
simres <- simulateResiduals(lrr.lm12)
plot(simres) # Weird

# Plot the LRR model
plot_model(
  lrr.lm12,
  type = "pred",
  terms= c("abs_spei12", "spei12_category"),
  show.data = TRUE
)

# LRR Nitrigen plots
# Add experiment column (sloppy)
# Add experiment column
plot_trt$experiment <- sub("nutnet.*", "nutnet", plot_trt$higher_order_organization)
plot_trt$experiment <- sub("glbrc_scaleup.*", "glbrc_scaleup", plot_trt$experiment)
plot_trt$experiment <- sub("glbrc_G10.*", "glbrc", plot_trt$experiment)
plot_trt$experiment <- sub("glbrc_G9.*", "glbrc", plot_trt$experiment)
plot_trt$experiment <- sub("mcse.*", "mcse", plot_trt$experiment)
plot_trt$experiment <- sub("microplots.*", "microplots", plot_trt$experiment)
plot_trt$experiment <- sub("Experiment 1.*", "Experiment 1", plot_trt$experiment)
plot_trt$experiment <- sub("001d_A_fl", "001d_fl", plot_trt$experiment)
plot_trt$experiment <- sub("001d_B_fl", "001d_fl", plot_trt$experiment)
plot_trt$experiment <- sub("001d_C_fl", "001d_fl", plot_trt$experiment)
plot_trt$experiment <- sub("001d_D_fl", "001d_fl", plot_trt$experiment)
plot_trt$experiment <- sub("001d_A_tu", "001d_tu", plot_trt$experiment)
plot_trt$experiment <- sub("001d_B_tu", "001d_tu", plot_trt$experiment)
plot_trt$experiment <- sub("001d_C_tu", "001d_tu", plot_trt$experiment)
plot_trt$experiment <- sub("001d_D_tu", "001d_tu", plot_trt$experiment)
plot_trt$experiment <- sub("004a_A_fl", "004a_fl", plot_trt$experiment)
plot_trt$experiment <- sub("004a_B_fl", "004a_fl", plot_trt$experiment)
plot_trt$experiment <- sub("004a_C_fl", "004a_fl", plot_trt$experiment)
plot_trt$experiment <- sub("004a_D_fl", "004a_fl", plot_trt$experiment)
plot_trt$experiment <- sub("004a_A_tu", "004a_tu", plot_trt$experiment)
plot_trt$experiment <- sub("004a_B_tu", "004a_tu", plot_trt$experiment)
plot_trt$experiment <- sub("004a_C_tu", "004a_tu", plot_trt$experiment)
plot_trt$experiment <- sub("004a_D_tu", "004a_tu", plot_trt$experiment)
plot_trt$experiment <- sub("004b_A_fl", "004b_fl", plot_trt$experiment)
plot_trt$experiment <- sub("004b_B_fl", "004b_fl", plot_trt$experiment)
plot_trt$experiment <- sub("004b_C_fl", "004b_fl", plot_trt$experiment)
plot_trt$experiment <- sub("004b_D_fl", "004b_fl", plot_trt$experiment)
plot_trt$experiment <- sub("004b_A_tu", "004b_tu", plot_trt$experiment)
plot_trt$experiment <- sub("004b_B_tu", "004b_tu", plot_trt$experiment)
plot_trt$experiment <- sub("004b_C_tu", "004b_tu", plot_trt$experiment)
plot_trt$experiment <- sub("004b_D_tu", "004b_tu", plot_trt$experiment)
plot_trt$experiment <- sub("002d.*", "002d", plot_trt$experiment)
plot_trt$experiment <- sub("Experiment 54.*", "Experiment 54", plot_trt$experiment)
plot_trt$experiment <- sub("KNZ_WAT01.*", "KNZ_WAT01", plot_trt$experiment)
plot_trt$experiment <- sub("002c.*", "002c", plot_trt$experiment)
plot_trt$experiment <- sub("e061.*", "e061", plot_trt$experiment)
plot_trt$experiment <- sub("e247.*", "e247", plot_trt$experiment)
plot_trt$experiment <- sub("e245.*", "e245", plot_trt$experiment)
plot_trt$experiment[plot_trt$experiment == "A"] <- "NGE"
plot_trt$experiment[plot_trt$experiment == "B"] <- "NGE"
plot_trt$experiment[plot_trt$experiment == "C"] <- "NGE"
plot_trt$experiment[plot_trt$experiment == "D"] <- "NGE"
plot_trt$experiment[plot_trt$experiment == "E"] <- "NGE"
plot_trt$experiment[plot_trt$experiment == "F"] <- "NGE"

# Subset nitrogen only
plot_n <- plot_trt %>%
  filter(nutrients_added == "NPK+" | nutrients_added == "N" | nutrients_added == "NP" | nutrients_added == "NPK" | nutrients_added == "NK" | nutrients_added == "NPK+Fence")

# Using SPEI3
plot_n_lrr3 <- plot_n %>%
  group_by(uniqueid) %>%
  arrange(uniqueid, year) %>%
  mutate(prior_year_biomass = lag(plot_biomass),
         prior_year_type = lag(spei3_category)) %>%
  ungroup()
normal_biomass3n <- plot_n_lrr3  %>%
  group_by(uniqueid) %>%
  filter(spei3_category == "Normal" | spei3_category == "Moderate wet" | spei3_category == "Moderate dry") %>%
  summarize(average_normal_biomass = mean(plot_biomass, na.rm=T))
plot_n_lrr3 <- left_join(plot_n_lrr3, normal_biomass3n)

lrr3n <- plot_n_lrr3 %>%
  filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
  filter(spei3_category == "Extreme wet" | spei3_category == "Extreme dry") %>%
  mutate(LRR = log(plot_biomass / average_normal_biomass))  # Calculate the log response ratio

# Take absolute value of SPEI3
lrr3n <- lrr3n %>%
  mutate(abs_spei3 = abs(spei3))

lrr.lm3n <- lmer(LRR ~ abs_spei3*spei3_category + (1|site/experiment/uniqueid) + (1|year), data = lrr3n)
summary(lrr.lm3n)
simres <- simulateResiduals(lrr.lm3n)
plot(simres) # Weird

# Plot the LRR model
plot_model(
  lrr.lm3n,
  type = "pred",
  terms= c("abs_spei3", "spei3_category"),
  show.data = TRUE
)

# Using SPEI6
plot_n_lrr6 <- plot_n %>%
  group_by(uniqueid) %>%
  arrange(uniqueid, year) %>%
  mutate(prior_year_biomass = lag(plot_biomass),
         prior_year_type = lag(spei6_category)) %>%
  ungroup()
normal_biomass6n <- plot_n_lrr6  %>%
  group_by(uniqueid) %>%
  filter(spei6_category == "Normal" | spei6_category == "Moderate wet" | spei6_category == "Moderate dry") %>%
  summarize(average_normal_biomass = mean(plot_biomass, na.rm=T))
plot_n_lrr6 <- left_join(plot_n_lrr6, normal_biomass6n)

lrr6n <- plot_n_lrr6 %>%
  filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
  filter(spei6_category == "Extreme wet" | spei6_category == "Extreme dry") %>%
  mutate(LRR = log(plot_biomass / average_normal_biomass))  # Calculate the log response ratio

# Take absolute value of SPEI6
lrr6n <- lrr6n %>%
  mutate(abs_spei6 = abs(spei6))

# Remove zero value in LRR
lrr6n <- lrr6n %>% filter(plot_biomass != 0)

lrr.lm6n <- lmer(LRR ~ abs_spei6*spei6_category + (1|site/experiment/uniqueid) + (1|year), data = lrr6n)
summary(lrr.lm6n)
simres <- simulateResiduals(lrr.lm6n)
plot(simres) # Weird

# Plot the LRR model
plot_model(
  lrr.lm6n,
  type = "pred",
  terms= c("abs_spei6", "spei6_category"),
  show.data = TRUE
)

# Using SPEI9
plot_n_lrr9 <- plot_n %>%
  group_by(uniqueid) %>%
  arrange(uniqueid, year) %>%
  mutate(prior_year_biomass = lag(plot_biomass),
         prior_year_type = lag(spei9_category)) %>%
  ungroup()
normal_biomass9n <- plot_n_lrr9  %>%
  group_by(uniqueid) %>%
  filter(spei9_category == "Normal" | spei9_category == "Moderate wet" | spei9_category == "Moderate dry") %>%
  summarize(average_normal_biomass = mean(plot_biomass, na.rm=T))
plot_n_lrr9 <- left_join(plot_n_lrr9, normal_biomass9n)

lrr9n <- plot_n_lrr9 %>%
  filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
  filter(spei9_category == "Extreme wet" | spei9_category == "Extreme dry") %>%
  mutate(LRR = log(plot_biomass / average_normal_biomass))  # Calculate the log response ratio

# Take absolute value of SPEI9
lrr9n <- lrr9n %>%
  mutate(abs_spei9 = abs(spei9))

# Remove zero value in LRR
lrr9n <- lrr9n %>% filter(plot_biomass != 0)

lrr.lm9n <- lmer(LRR ~ abs_spei9*spei9_category + (1|site/experiment/uniqueid) + (1|year), data = lrr9n)
summary(lrr.lm9n)
simres <- simulateResiduals(lrr.lm9n)
plot(simres) # Weird

# Plot the LRR model
plot_model(
  lrr.lm9n,
  type = "pred",
  terms= c("abs_spei9", "spei9_category"),
  show.data = TRUE
)

# Using SPEI12
plot_n_lrr12 <- plot_n %>%
  group_by(uniqueid) %>%
  arrange(uniqueid, year) %>%
  mutate(prior_year_biomass = lag(plot_biomass),
         prior_year_type = lag(spei12_category)) %>%
  ungroup()
normal_biomass12n <- plot_n_lrr12  %>%
  group_by(uniqueid) %>%
  filter(spei12_category == "Normal" | spei12_category == "Moderate wet" | spei12_category == "Moderate dry") %>%
  summarize(average_normal_biomass = mean(plot_biomass, na.rm=T))
plot_n_lrr12 <- left_join(plot_n_lrr12, normal_biomass12n)

lrr12n <- plot_n_lrr12 %>%
  filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
  filter(spei12_category == "Extreme wet" | spei12_category == "Extreme dry") %>%
  mutate(LRR = log(plot_biomass / average_normal_biomass))  # Calculate the log response ratio

# Take absolute value of SPEI12
lrr12n <- lrr12n %>%
  mutate(abs_spei12 = abs(spei12))

# Remove zero value in LRR
lrr12n <- lrr12n %>% filter(plot_biomass != 0)

lrr.lm12n <- lmer(LRR ~ abs_spei12*spei12_category + (1|site/experiment/uniqueid) + (1|year), data = lrr12n, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) # Failed to converge with defaults
summary(lrr.lm12n)
simres <- simulateResiduals(lrr.lm12n)
plot(simres) # Weird

# Plot the LRR model
plot_model(
  lrr.lm12n,
  type = "pred",
  terms= c("abs_spei12", "spei12_category"),
  show.data = TRUE
)


# LRR with nitrogen yes/no as a predictor
plot_nc <- plot_trt %>%
  filter(nutrients_added == "NPK+" | nutrients_added == "N" | nutrients_added == "NP" | nutrients_added == "NPK" | nutrients_added == "NK" | nutrients_added == "NPK+Fence" | nutrients_added == "none" | nutrients_added == "no_fertilizer")

# Add column with nitrogen yes/no
plot_nc <- plot_nc %>%
  mutate(nitrogen = recode(nutrients_added, "NPK+" = "N", "NP" = "N", "NPK" = "N", "NK" = "N", "NPK+Fence" = "N", "none" = "no_fertilizer"))

# Using SPEI3
plot_nc_lrr3 <- plot_nc %>%
  group_by(uniqueid) %>%
  arrange(uniqueid, year) %>%
  mutate(prior_year_biomass = lag(plot_biomass),
         prior_year_type = lag(spei3_category)) %>%
  ungroup()
normal_biomass3nc <- plot_nc_lrr3  %>%
  group_by(uniqueid) %>%
  filter(spei3_category == "Normal" | spei3_category == "Moderate wet" | spei3_category == "Moderate dry") %>%
  summarize(average_normal_biomass = mean(plot_biomass, na.rm=T))
plot_nc_lrr3 <- left_join(plot_nc_lrr3, normal_biomass3nc)

lrr3nc <- plot_nc_lrr3 %>%
  filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
  filter(spei3_category == "Extreme wet" | spei3_category == "Extreme dry") %>%
  mutate(LRR = log(plot_biomass / average_normal_biomass))  # Calculate the log response ratio

# Take absolute value of SPEI3
lrr3nc <- lrr3nc %>%
  mutate(abs_spei3 = abs(spei3))

lrr.lm3nc <- lmer(LRR ~ abs_spei3*spei3_category*nitrogen + (1|site/experiment/uniqueid) + (1|year), data = lrr3nc)
summary(lrr.lm3nc)
Anova(lrr.lm3nc, type = "III")
simres <- simulateResiduals(lrr.lm3nc)
plot(simres)

# Plot the LRR model
plot_model(
  lrr.lm3nc,
  type = "pred",
  terms= c("abs_spei3", "spei3_category", "nitrogen"),
  show.data = TRUE
)

# Using SPEI6
plot_nc_lrr6 <- plot_nc %>%
  group_by(uniqueid) %>%
  arrange(uniqueid, year) %>%
  mutate(prior_year_biomass = lag(plot_biomass),
         prior_year_type = lag(spei6_category)) %>%
  ungroup()
normal_biomass6nc <- plot_nc_lrr6  %>%
  group_by(uniqueid) %>%
  filter(spei6_category == "Normal" | spei6_category == "Moderate wet" | spei6_category == "Moderate dry") %>%
  summarize(average_normal_biomass = mean(plot_biomass, na.rm=T))
plot_nc_lrr6 <- left_join(plot_nc_lrr6, normal_biomass6nc)

lrr6nc <- plot_nc_lrr6 %>%
  filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
  filter(spei6_category == "Extreme wet" | spei6_category == "Extreme dry") %>%
  mutate(LRR = log(plot_biomass / average_normal_biomass))  # Calculate the log response ratio

# Take absolute value of SPEI6
lrr6nc <- lrr6nc %>%
  mutate(abs_spei6 = abs(spei6))

# Remove zero value in LRR
lrr6nc <- lrr6nc %>% filter(plot_biomass != 0)

lrr.lm6nc <- lmer(LRR ~ abs_spei6*spei6_category*nitrogen + (1|site/experiment/uniqueid) + (1|year), data = lrr6nc)
summary(lrr.lm6nc)
Anova(lrr.lm6nc, type = "III")
simres <- simulateResiduals(lrr.lm6nc)
plot(simres)

# Plot the LRR model
plot_model(
  lrr.lm6nc,
  type = "pred",
  terms= c("abs_spei6", "spei6_category", "nitrogen"),
  show.data = TRUE
)

# Using SPEI9
plot_nc_lrr9 <- plot_nc %>%
  group_by(uniqueid) %>%
  arrange(uniqueid, year) %>%
  mutate(prior_year_biomass = lag(plot_biomass),
         prior_year_type = lag(spei9_category)) %>%
  ungroup()
normal_biomass9nc <- plot_nc_lrr9  %>%
  group_by(uniqueid) %>%
  filter(spei9_category == "Normal" | spei9_category == "Moderate wet" | spei9_category == "Moderate dry") %>%
  summarize(average_normal_biomass = mean(plot_biomass, na.rm=T))
plot_nc_lrr9 <- left_join(plot_nc_lrr9, normal_biomass9nc)

lrr9nc <- plot_nc_lrr9 %>%
  filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
  filter(spei9_category == "Extreme wet" | spei9_category == "Extreme dry") %>%
  mutate(LRR = log(plot_biomass / average_normal_biomass))  # Calculate the log response ratio

# Take absolute value of SPEI9
lrr9nc <- lrr9nc %>%
  mutate(abs_spei9 = abs(spei9))

# Remove zero value in LRR
lrr9nc <- lrr9nc %>% filter(plot_biomass != 0)

lrr.lm9nc <- lmer(LRR ~ abs_spei9*spei9_category*nitrogen + (1|site/experiment/uniqueid) + (1|year), data = lrr9nc)
summary(lrr.lm9nc)
Anova(lrr.lm9nc, type = "III")
simres <- simulateResiduals(lrr.lm9nc)
plot(simres)

# Plot the LRR model
plot_model(
  lrr.lm9nc,
  type = "pred",
  terms= c("abs_spei9", "spei9_category", "nitrogen"),
  show.data = TRUE
)

# Using SPEI12
plot_nc_lrr12 <- plot_nc %>%
  group_by(uniqueid) %>%
  arrange(uniqueid, year) %>%
  mutate(prior_year_biomass = lag(plot_biomass),
         prior_year_type = lag(spei12_category)) %>%
  ungroup()
normal_biomass12nc <- plot_nc_lrr12  %>%
  group_by(uniqueid) %>%
  filter(spei12_category == "Normal" | spei12_category == "Moderate wet" | spei12_category == "Moderate dry") %>%
  summarize(average_normal_biomass = mean(plot_biomass, na.rm=T))
plot_nc_lrr12 <- left_join(plot_nc_lrr12, normal_biomass12nc)

lrr12nc <- plot_nc_lrr12 %>%
  filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
  filter(spei12_category == "Extreme wet" | spei12_category == "Extreme dry") %>%
  mutate(LRR = log(plot_biomass / average_normal_biomass))  # Calculate the log response ratio

# Take absolute value of SPEI12
lrr12nc <- lrr12nc %>%
  mutate(abs_spei12 = abs(spei12))

# Remove zero value in LRR
lrr12nc <- lrr12nc %>% filter(plot_biomass != 0)

lrr.lm12nc <- lmer(LRR ~ abs_spei12*spei12_category*nitrogen + (1|site/experiment/uniqueid) + (1|year), data = lrr12nc)
summary(lrr.lm12nc)
Anova(lrr.lm12nc, type = "III")
simres <- simulateResiduals(lrr.lm12nc)
plot(simres)

# Plot the LRR model
plot_model(
  lrr.lm12nc,
  type = "pred",
  terms= c("abs_spei12", "spei12_category", "nitrogen"),
  show.data = TRUE
)

