# TITLE:        LTER Grassland Rock: Create core analyses 
# AUTHORS:      Ashley Darst
# COLLABORATORS:  
# DATA INPUT:   Data imported as csv files from shared Google drive L2 folder
# DATA OUTPUT:  Core analyses
# PROJECT:      LTER Grassland Rock
# DATE:         October 2024

# Clear all existing data
rm(list=ls())

# Load packages
library(tidyverse)
library(lme4)
library(lmerTest)
library(janitor)
library(DHARMa)
library(bbmle)
library(ggeffects)
library(MuMIn)

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

# Merge with metadata
plot_ece_meta <- left_join(plot_ece, meta)

# Remove NAs for non-extreme years
plot_ece_rm_na <- plot_ece_meta %>%
  drop_na(resistance)

# Subset to only have control plots
plot_ece_control <- plot_ece_rm_na %>%
  filter(treatment == "control")

# Analysis 1: resistance vs richness ----
## Control plot only ----

# Model with all variables and maximal random effects (need to add experiment)
resist.control.full <- lmer(log10(resistance) ~ scale(richness)*scale(berger_parker)*site*spei6_category + year + (1|site/uniqueid), data = plot_ece_control, na.action = na.fail) # failed to converge
summary(resist.control.full)
plot(resist.control.full)
simres <- simulateResiduals(resist.control.full)
plot(simres) # Very bad!

# dredge(resist.control.full) # best model only has spei_category lol

# Additive only model
resist.control.add <- lmer(log10(resistance) ~ richness + berger_parker + site + spei6_category + year + (1|site/uniqueid), data = plot_ece_control) # failed to converge

# Three way interaction, drop site, scale richness and dominance
resist.control.int <- lmer(log10(resistance) ~ scale(richness)*scale(berger_parker)*spei6_category + year + (1|site/uniqueid), data = plot_ece_control)
summary(resist.control.int)
simres <- simulateResiduals(resist.control.int)
plot(simres) # Very bad!

# Drop non significant interactions (need to compare AIC)
resist.control <- lmer(log10(resistance) ~ scale(richness)*spei6_category + scale(berger_parker) + year + (1|site/uniqueid), data = plot_ece_control)
summary(resist.control)
simres <- simulateResiduals(resist.control)
plot(simres) # Very bad!

# Make a model based on hypotheses
resist.control.hyp <- lmer(log10(resistance) ~ scale(richness)*spei6_category + scale(berger_parker)*scale(richness) + scale(berger_parker)*spei6_category + year + (1|site/uniqueid), data = plot_ece_control)
summary(resist.control.hyp)
simres <- simulateResiduals(resist.control.hyp)
plot(simres) # Very bad!

AICctab(resist.control.int, resist.control, resist.control.hyp)

# Plot best model
ggpredict(model = resist.control, terms = c("richness", "spei6_category"), back_transform = F) %>%
  plot(show_data = TRUE)




## All plots ----
# Model 1
resist.lm <- lmer(resistance ~ richness + (1|site) + (1|site:plot) + (1|year), data = plot_ece)
summary(resist.lm) # boundary singular fit, site explains 0 variance

# Check residuals and QQplot
simres <- simulateResiduals(resist.lm)
plot(simres) # Very bad!

# Model 2 - try logging resistance
resist.lm.log <- lmer(log10(resistance) ~ richness + (1|site) + (1|site:plot) + (1|year), data = plot_ece)
summary(resist.lm.log)

# Check residuals and QQplot
simres <- simulateResiduals(resist.lm.log)
plot(simres) # Better, but not perfect

# Model 3 - try logging richness
resist.lm.log2 <- lmer(log10(resistance) ~ log10(richness) + (1|site) + (1|site:plot) + (1|year), data = plot_ece)
summary(resist.lm.log2)

# Check residuals and QQplot
simres <- simulateResiduals(resist.lm.log2)
plot(simres) # Didn't help

# Model 4 - add wet vs dry event to model 2
resist.lm.wd <- lmer(log10(resistance) ~ richness + spei6_category + (1|site) + (1|site:plot) + (1|year), data = plot_ece)
summary(resist.lm.wd)
resist.lm.wd.log <- lmer(log10(resistance) ~ log10(richness) + spei6_category + (1|site) + (1|site:plot) + (1|year), data = plot_ece)

# Check residuals and QQplot
simres <- simulateResiduals(resist.lm.wd)
plot(simres) # Didn't help

# Model 5 - try interaction with category
resist.lm.wd.int <- lmer(log10(resistance) ~ richness*spei6_category + (1|site) + (1|site:plot) + (1|year), data = plot_ece)
summary(resist.lm.wd.int)
resist.lm.wd.int.log <- lmer(log10(resistance) ~ log10(richness)*spei6_category + (1|site) + (1|site:plot) + (1|year), data = plot_ece)

# Check residuals and QQplot
simres <- simulateResiduals(resist.lm.wd.int)
plot(simres) # Didn't help

# Model 6 - add dominance
resist.lm.dom <- lmer(log10(resistance) ~ richness*spei6_category + berger_parker + (1|site) + (1|site:plot) + (1|year), data = plot_ece)
summary(resist.lm.dom)
resist.lm.dom.log <- lmer(log10(resistance) ~ log10(richness)*spei6_category + berger_parker + (1|site) + (1|site:plot) + (1|year), data = plot_ece)

# Check residuals and QQplot
simres <- simulateResiduals(resist.lm.dom)
plot(simres) # Didn't help

# Model 7 - add dominance interaction with richness
resist.lm.dom.int <- lmer(log10(resistance) ~ richness*spei6_category + berger_parker*richness + (1|site) + (1|site:plot) + (1|year), data = plot_ece)
summary(resist.lm.dom.int)
resist.lm.dom.int.log <- lmer(log10(resistance) ~ log10(richness)*spei6_category + berger_parker*log10(richness) + (1|site) + (1|site:plot) + (1|year), data = plot_ece)

# Check residuals and QQplot
simres <- simulateResiduals(resist.lm.dom.int)
plot(simres) # Didn't help

# Compare models # Depends on how/if you log richness
AICctab(resist.lm.log, resist.lm.log2, resist.lm.wd, resist.lm.wd.log, resist.lm.wd.int, resist.lm.wd.int.log, resist.lm.dom, resist.lm.dom.log, resist.lm.dom.int, resist.lm.dom.int.log)

# Plot best model
ggpredict(model = resist.lm.log2, terms = c("richness"), back_transform = F) %>%
  plot(show_data = TRUE)

# Random effects like Isbell # Not sure if this is right
resist.lm.r <- lmer(log10(resistance) ~ richness*spei6_category + (1|site:richness:year) + (1|site:richness) + (1|site:year) + (1|site:plot) + (1|site), data = plot_ece, na.action=na.omit)
summary(resist.lm.r)

resist.lm.r2 <- lmer(log10(resistance) ~ richness + spei6_category + (1|site:richness:year) + (1|site:richness) + (1|site:year) + (1|site:plot) + (1|site), data = plot_ece, na.action=na.omit)
resist.lm.r2.log <- lmer(log10(resistance) ~ log10(richness) + spei6_category + (1|site:richness:year) + (1|site:richness) + (1|site:year) + (1|site:plot) + (1|site), data = plot_ece, na.action=na.omit)
resist.lm.r2.log2 <- lmer(log2(resistance) ~ log2(richness) + spei6_category + (1|site:richness:year) + (1|site:richness) + (1|site:year) + (1|site:plot) + (1|site), data = plot_ece, na.action=na.omit)
resist.lm.r3 <- lmer(log10(resistance) ~ richness + spei6_category + berger_parker + (1|site:richness:year) + (1|site:richness) + (1|site:year) + (1|site:plot) + (1|site), data = plot_ece, na.action=na.omit)
resist.lm.r4 <- lmer(log10(resistance) ~ richness*spei6_category + berger_parker + (1|site:richness:year) + (1|site:richness) + (1|site:year) + (1|site:plot) + (1|site), data = plot_ece, na.action=na.omit)
resist.lm.r5 <- lmer(log10(resistance) ~ richness*spei6_category + berger_parker*spei6_category + (1|site:richness:year) + (1|site:richness) + (1|site:year) + (1|site:plot) + (1|site), data = plot_ece, na.action=na.omit)

# Compare models with Isbell random effects
AICctab(resist.lm.r, resist.lm.r2, resist.lm.r2.log, resist.lm.r2.log2, resist.lm.r3, resist.lm.r4, resist.lm.r5)

# Look at model stats
summary(resist.lm.r2.log)
simres <- simulateResiduals(resist.lm.r2.log)
plot(simres) # Not great

# Plot best model
ggpredict(model = resist.lm.r2.log, terms = c("richness", "spei6_category"), back_transform = F) %>%
  plot(show_data = TRUE)


## Analysis 2 - resilience vs richness ----

# Model 1
resil.lm <- lmer(resilience ~ richness + (1|site) + (1|site:plot), data = plot_ece)
summary(resil.lm) # boundary singular fit, site and plot explain 0 variance

# Check residuals and QQplot
simres <- simulateResiduals(resil.lm)
plot(simres) # Very bad!

# Model 2 - try logging resistance
resil.lm.log <- lmer(log10(resilience) ~ richness + (1|site) + (1|site:plot), data = plot_ece)
summary(resil.lm.log) # boundary singular fit

# Check residuals and QQplot
simres <- simulateResiduals(resil.lm.log)
plot(simres) # Better, but not perfect

# Model 3 - try logging richness
resil.lm.log2 <- lmer(log10(resilience) ~ log10(richness) + (1|site) + (1|site:plot), data = plot_ece)
summary(resil.lm.log2)

# Check residuals and QQplot
simres <- simulateResiduals(resil.lm.log2)
plot(simres) # Didn't help

# Model 4 - add wet vs dry event to model 2
resil.lm.wd <- lmer(log10(resilience) ~ richness + spei6_category + (1|site) + (1|site:plot), data = plot_ece)
summary(resil.lm.wd)

# Model 5 - try interaction with category
resil.lm.wd.int <- lmer(log10(resilience) ~ richness*spei6_category + (1|site) + (1|site:plot), data = plot_ece)
summary(resil.lm.wd.int)

# Model 6 - add dominance
resil.lm.dom <- lmer(log10(resilience) ~ richness*spei6_category + berger_parker + (1|site) + (1|site:plot), data = plot_ece)
summary(resil.lm.dom)

# Model 7 - add dominance interaction with richness
resil.lm.dom.int <- lmer(log10(resilience) ~ richness*spei6_category + berger_parker*richness + (1|site) + (1|site:plot), data = plot_ece)
summary(resil.lm.dom.int)

# Compare models 
AICctab(resil.lm.log, resil.lm.wd, resil.lm.wd.int, resil.lm.dom, resil.lm.dom.int)

# Plot best model
ggpredict(model = resil.lm.log, terms = c("richness"), back_transform = F) %>%
  plot(show_data = TRUE)

