# TITLE:        LTER Grassland Rock: Core analyses then MGMT factors
# AUTHORS:      Ashley Darst, then Caitlin Broderick
# COLLABORATORS:  
# DATA INPUT:   Data imported as csv files from shared Google drive L2 folder
# DATA OUTPUT:  Core analyses
# PROJECT:      LTER Grassland Rock
# DATE:         January 2025

# all cleanup code, and initial models, are Ashleys!!!!
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
library(emmeans)

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

# Add experiment column
plot_ece$experiment <- sub("nutnet.*", "nutnet", plot_ece$higher_order_organization)
plot_ece$experiment <- sub("glbrc_scaleup.*", "glbrc_scaleup", plot_ece$experiment)
plot_ece$experiment <- sub("glbrc_G10.*", "glbrc", plot_ece$experiment)
plot_ece$experiment <- sub("glbrc_G9.*", "glbrc", plot_ece$experiment)
plot_ece$experiment <- sub("mcse.*", "mcse", plot_ece$experiment)
plot_ece$experiment <- sub("microplots.*", "microplots", plot_ece$experiment)
plot_ece$experiment <- sub("Experiment 1.*", "Experiment 1", plot_ece$experiment)
plot_ece$experiment <- sub("001d.*", "001d", plot_ece$experiment)
plot_ece$experiment <- sub("004b.*", "004b", plot_ece$experiment)
plot_ece$experiment <- sub("002d.*", "002d", plot_ece$experiment)
plot_ece$experiment <- sub("Experiment 54.*", "Experiment 54", plot_ece$experiment)
plot_ece$experiment <- sub("KNZ_WAT01.*", "KNZ_WAT01", plot_ece$experiment)
plot_ece$experiment <- sub("002c.*", "002c", plot_ece$experiment)
plot_ece$experiment <- sub("e061.*", "e061", plot_ece$experiment)
plot_ece$experiment <- sub("e247.*", "e247", plot_ece$experiment)
plot_ece$experiment <- sub("e245.*", "e245", plot_ece$experiment)

# Merge with metadata
plot_ece_meta <- left_join(plot_ece, meta)

plot_ece_meta  %>% filter(is.na(resistance)) %>% nrow() # 16k NA
plot_ece_meta  %>% filter(!is.na(resistance)) %>% nrow() #4k are not 

# Remove GRAZED PLOTS!  new addition after Konza talk.
table(plot_ece_meta$grazing)
plot_ece_nograze <- plot_ece_meta %>%
  filter(grazing =="ungrazed")


# Remove NAs for non-extreme years
plot_ece_rm_na <- plot_ece_nograze %>%
  drop_na(resistance)

# Make year a factor
str(plot_ece_rm_na)
plot_ece_rm_na$year <- as.factor(plot_ece_rm_na$year)


# make "summary" environmental predictor categories 
plot_ece_rm_na$fire_frequency_cat = dplyr::case_when(plot_ece_rm_na$fire_frequency==0 ~ "unburned",
                                                     is.na(plot_ece_rm_na$fire_frequency )== TRUE ~ "unburned" ,  
                                                     plot_ece_rm_na$fire_frequency>0 & plot_ece_rm_na$fire_frequency < 4 ~"frequent burn",
                                                     plot_ece_rm_na$fire_frequency >= 4 ~ "infrequent burn")

plot_ece_rm_na$nutrients_simple = dplyr::case_when(plot_ece_rm_na$nutrients_added=="no_fertilizer" ~ "no fertilizer",
                                                   plot_ece_rm_na$nutrients_added!="no fertilizer" ~ "fertilized")



# how closely related are different measures of SPEI
plot(plot$spei6, plot$spei12) # very
cor.test(plot$spei6, plot$spei12) # pretty correlated.....

# how correlated are our main predictors, rich and dominance?
plot(plot$Berger_Parker, plot$Richness) # kinda 
cor.test(plot$Berger_Parker, plot$Richness) # 0.6 , so like yes but not enough that we cant have both in model.

cor.test(log(plot$Berger_Parker), plot$Richness)


###################################
# Analysis 1: resistance ----
###################################
## Control plot only ----

# Note: took out all the models that did not converge. 
# Note: Year should be a random effect in my opinion, otherwise uninterpretable. It has to be, right?

# Big model with three way interactions... Three way interaction, drop site, scale richness and dominance - SINGULAR FIT
resist.control.int <- lmer(log10(resistance) ~ scale(richness)*scale(berger_parker)*spei6_category + (1|year) + (1|site/uniqueid), data = plot_ece_control)
summary(resist.control.int)
simres <- simulateResiduals(resist.control.int)
plot(simres) # Not great


# Drop dominance (but it's in our hypothesis so we might want to include it) SINGULAR FIT. best model. 
head(plot_ece_control$uniqueid)
resist.control.rich.ny <- lmer(log10(resistance) ~ scale(richness)*spei6_category + (1|year) + (1|site/uniqueid), data = plot_ece_control)
summary(resist.control.rich.ny)
simres <- simulateResiduals(resist.control.rich.ny)
plot(simres) # Very bad!


# Make a model based on hypotheses - year random! MY personal favorite. 
resist.control.hyp.ny <- lmer(log10(resistance) ~ scale(richness)*spei6_category + scale(berger_parker)*scale(richness) + scale(berger_parker)*spei6_category + (1|year) + (1|site/uniqueid), data = plot_ece_control)
summary(resist.control.hyp.ny)
simres <- simulateResiduals(resist.control.hyp.ny)
plot(simres) # Very bad!

resist.control.hyp.ny.emmeans = emmeans(resist.control.hyp.ny, ~scale(richness)|spei6_category)
summary(resist.control.hyp.ny.emmeans)


# Plot best model
ggpredict(model = resist.control.rich.ny, terms = c("richness", "spei6_category"), back_transform = F) %>%
  plot(show_data = TRUE)

ggpredict(model = resist.control.hyp.ny, terms = c("richness", "spei6_category"), back_transform = F) %>%
  plot(show_data = TRUE) # plot hypothesis model. 


## All plots ----

# no environmental yet!
# Model Hypothesis : Rich*cat , dom*cat 
resist.lm.hyp <- lmer(log10(resistance) ~ scale(richness)*spei6_category + scale(berger_parker)*spei6_category + (1|site) + (1|site:plot) + (1|year), data = plot_ece_rm_na)
summary(resist.lm.hyp)
resist.lm.hyp.log <- lmer(log10(resistance) ~ log10(richness)*spei6_category + berger_parker*spei6_category + (1|site) + (1|site:plot) + (1|year), data = plot_ece_rm_na)
summary(resist.lm.hyp.log)
resist.lm.hyp.log2 <- lmer(log10(resistance) ~ log10(richness)*spei6_category + log10(berger_parker)*spei6_category + (1|site) + (1|site:plot) + (1|year), data = plot_ece_rm_na)
summary(resist.lm.hyp.log2)

# Check residuals and QQplot
simres <- simulateResiduals(resist.lm.hyp)
plot(simres) 
simres <- simulateResiduals(resist.lm.hyp.log)
plot(simres) # Didn't help

resist.lm.dom.emmeans = emmeans(resist.lm.dom, ~scale(richness)|spei6_category)
summary(resist.lm.dom.emmeans)

resist.lm.hyp.log.emmeans = emmeans(resist.lm.hyp.log, ~scale(richness)|spei6_category)
summary(resist.lm.hyp.log.emmeans)

# Compare models # Depends on how/if you log richness
AICctab(resist.lm.hyp,resist.lm.hyp.log,resist.lm.hyp.log2) # interesting, logging richness is good.

# Plot  model
ggpredict(model = resist.lm.hyp, terms = c("richness", "spei6_category"), back_transform = F) %>%
  plot(show_data = TRUE)

ggpredict(model = resist.lm.hyp.log, terms = c("richness", "spei6_category"), back_transform = F) %>%
  plot(show_data = TRUE) # wack. 


# TRY ADDING IN OTHER PREDICTORS !!!!!

# first off 
# additive environment model for RICHNESS - are environment shaping our predictors?
richess.env.add = lmer(log10(richness) ~ nutrients_simple + fire_frequency_cat+ disturbance+ (1|site) + (1|site:plot) + (1|year), data = plot_ece_rm_na)
summary(richess.env.add ) # yes they do esp nutrients. 

#interaction env model for RICHNESS
richess.env = lmer(richness ~ nutrients_simple* fire_frequency_cat* disturbance+ (1|site) + (1|site:plot) + (1|year), data = plot_ece_rm_na)
summary(richess.env ) # yes they do esp nutrients. 

# additive environment model for DOMINANCE - are environemnt shaping our predictors?
dom.env.add = lmer(log10(berger_parker) ~ nutrients_simple + fire_frequency_cat+ disturbance+ (1|site) + (1|site:plot) + (1|year), data = plot_ece_rm_na)
summary(dom.env.add ) # everything

#interaction env model for DOMINANCE
dom.env = lmer(berger_parker ~ nutrients_simple* fire_frequency_cat* disturbance+ (1|site) + (1|site:plot) + (1|year), data = plot_ece_rm_na)
summary(dom.env ) # yes everything. 


######### looking at effects on resistance .
#first off: do other predictors influence dominance or resistance?
# environmental model for RESISTANCE! 
resist.env = lm(log10(resistance) ~ nutrients_simple*fire_frequency_cat , data = plot_ece_rm_na) # removed disturbance too few dissturbed 
summary(resist.env)
plot_ece_rm_na_wresidresist  = plot_ece_rm_na
plot_ece_rm_na_wresidresist$residresist = resid(resist.env)
hist(plot_ece_rm_na_wresidresist$residresist)

# run model on resistance residuals <- from environment 
resist.lm.hyp.resid <- lmer(residresist~ scale(richness)*spei6_category + scale(berger_parker)*spei6_category + (1|site) + (1|site:plot) + (1|year), data = plot_ece_rm_na_wresidresist)
summary(resist.lm.hyp.resid)


plot_ece_rm_na_wresidresist %>% ggplot(aes( x= richness, y=residresist, color=spei6_category)) +
  geom_point() + geom_smooth(method = "lm")


# incorporate dry wet
resist.env.dw = lmer(log10(resistance) ~ nutrients_added*spei6_category + fire_frequency_cat*spei6_category+ disturbance*spei6_category+ (1|site) + (1|site:plot) + (1|year), data = plot_ece_rm_na)
summary(resist.env.dw) # just overall. nutrients and disturbance matter....

emmeans(resist.env.dw  ,~fire_frequency_cat|spei6_category)


plot_ece_rm_na %>% ggplot(aes( x= richness, y=log10(resistance), color=nutrients_simple)) +
  geom_point() + geom_smooth(method = "lm")


resist.env.dw_richdom = lmer(log10(resistance) ~ scale(richness)*spei6_category + scale(berger_parker)*spei6_category+ nutrients_simple*spei6_category + fire_frequency_cat*spei6_category+ disturbance*spei6_category+ (1|site) + (1|site:plot) + (1|year), data = plot_ece_rm_na)
summary(resist.env.dw_richdom)


table(is.na(plot_ece_rm_na$resistance))
table(is.na(plot_ece_rm_na$richness))
table(is.na(plot_ece_rm_na$berger_parker))
table(is.na(plot_ece_rm_na$nutrients_simple))
table(is.na(plot_ece_rm_na$fire_frequency_cat))
table(is.na(plot_ece_rm_na$spei6_category))
table(is.na(plot_ece_rm_na$disturbance))


resist.env.dw_full = lmer(log10(resistance) ~ scale(richness)*spei6_category*nutrients_simple*fire_frequency_cat*spei6_category +
                            scale(berger_parker)*spei6_category*nutrients_simple*fire_frequency_cat*spei6_category +(1|site:plot) + (1|year), data = plot_ece_rm_na, na.action = "na.fail")
summary(resist.env.dw_full)

#dredged = dredge(resist.env.dw_full )
dredged # this does not include iether richness or berger parker !!!
# fire freq, nutrient, sp6. fire*nutrient, fire*sp6
# delta difference is 2.41 . 
# note: this did not include disturbance!!!


# ONE WAY TO DO THIS BETTER: SPLIT UP BY SPEI6 Category!!!
  # way easier than having huge messy interaction !!!
plot_ece_rm_na_dry = plot_ece_rm_na %>% filter(spei6_category=="Extreme dry")
plot_ece_rm_na_wet = plot_ece_rm_na %>% filter(spei6_category=="Extreme wet")



# dry years
resist.env.dw_dry = lmer(log10(resistance) ~ scale(richness)*nutrients_simple*fire_frequency_cat +
                            scale(berger_parker)*nutrients_simple*fire_frequency_cat +(1|site:plot) + (1|year), data = plot_ece_rm_na_dry, na.action = "na.fail")
summary(resist.env.dw_dry)
# do emtrends here?

plot_ece_rm_na_dry %>%  ggplot(aes (x=richness, y = log10(resistance), color=nutrients_simple) ) +
  geom_point() + geom_smooth(method = "lm") + 
  facet_grid(~fire_frequency_cat)

plot_ece_rm_na_dry %>%  ggplot(aes (x=berger_parker, y = log10(resistance), color=nutrients_simple) ) +
  geom_point() + geom_smooth(method = "lm") + 
  facet_grid(~fire_frequency_cat)

ggpredict(model = resist.env.dw_dry , terms = c("richness", "nutrients_simple", "fire_frequency_cat"), back_transform = F) %>%
  plot(show_data = TRUE) # plot model - different than data

ggpredict(model = resist.env.dw_dry , terms = c("berger_parker", "nutrients_simple", "fire_frequency_cat"), back_transform = F) %>%
  plot(show_data = TRUE) # plot model-  different than data 


# wet years
resist.env.dw_wet = lmer(log10(resistance) ~ scale(richness)*nutrients_simple*fire_frequency_cat +
                           scale(berger_parker)*nutrients_simple*fire_frequency_cat +(1|site:plot) + (1|year), data = plot_ece_rm_na_wet, na.action = "na.fail")
summary(resist.env.dw_wet)




# NOW ACTUALLY DO PSEM - this is not working
  # to do - meet iwth Maatt av
library(piecewiseSEM)
plot_ece_rm_na_wet
str(scale(plot_ece_rm_na$berger_parker) )
scale(plot_ece_rm_na$richness) 
modelList_lm <- psem(
  lm(log10(resistance) ~ scale(richness)*spei6_category + scale(berger_parker)*spei6_category + nutrients_simple + fire_frequency_cat , data=plot_ece_rm_na),
  lm(richness ~ nutrients_simple*fire_frequency_cat, data=plot_ece_rm_na),
  lm(berger_parker ~ nutrients_simple*fire_frequency_cat, data=plot_ece_rm_na),
  data=plot_ece_rm_na
)


modelList_lm
basisSet(modelList_lm)
dSep(modelList_lm, .progressBar = FALSE)
summary(modelList_lm )


# model for wet 
names( plot_ece_rm_na_wet)
plot_ece_rm_na_wet_sel = plot_ece_rm_na_wet %>% select(c(richness,berger_parker,nutrients_simple,fire_frequency_cat,resistance))
dim(plot_ece_rm_na_wet_sel)
table(is.na(plot_ece_rm_na_wet_sel))

modelList_lm_wet <- psem(
  lm(log10(resistance) ~ scale(richness) + scale(berger_parker) + nutrients_simple + fire_frequency_cat , data=plot_ece_rm_na_wet_sel),
  lm(richness ~ nutrients_simple*fire_frequency_cat, data=plot_ece_rm_na_wet_sel),
  lm(berger_parker ~ nutrients_simple*fire_frequency_cat, data=plot_ece_rm_na_wet_sel),
  data=plot_ece_rm_na_wet_sel
)

modelList_lm_wet
summary(modelList_lm_wet) # ruh roh
plot(modelList_lm_wet)



# model for dry
names( plot_ece_rm_na_dry)
plot_ece_rm_na_dry_sel = plot_ece_rm_na_dry %>% select(c(richness,berger_parker,nutrients_simple,fire_frequency_cat,resistance))
dim(plot_ece_rm_na_dry_sel)
table(is.na(plot_ece_rm_na_dry_sel))

modelList_lm_dry <- psem(
  lm(log10(resistance) ~ scale(richness) + scale(berger_parker) + nutrients_simple + fire_frequency_cat , data=plot_ece_rm_na_dry_sel),
  lm(richness ~ nutrients_simple*fire_frequency_cat, data=plot_ece_rm_na_dry_sel),
  lm(berger_parker ~ nutrients_simple*fire_frequency_cat, data=plot_ece_rm_na_dry_sel),
  data=plot_ece_rm_na_dry_sel
)

modelList_lm_dry
summary(modelList_lm_dry) # ruh roh
plot(modelList_lm_dry)




#####################
# Analysis 2 - resilience ----
#####################
## Control plots only ----

# note: failed to converge with multiple-way interaction. 
dim(plot_ece_control)


# Move year to random effects
resil.control.rich.ny <- lmer(log10(resilience) ~ scale(richness)*spei6_category + (1|year) + (1|site/uniqueid), data = plot_ece_control) # singular fit
summary(resil.control.rich.ny)
simres <- simulateResiduals(resil.control.rich.ny)
plot(simres) # Very bad!


# Three way interaction (year as random effect)
resil.control.int.ny <- lmer(log10(resilience) ~ scale(richness)*scale(berger_parker)*spei6_category + (1|year) + (1|site/uniqueid), data = plot_ece_control)
summary(resil.control.int)

# Drop non significant interactions (need to compare AIC)
resil.control <- lmer(log10(resilience) ~ scale(richness)*spei6_category + scale(berger_parker)*spei6_category + year + (1|site/uniqueid), data = plot_ece_control)
summary(resil.control)
simres <- simulateResiduals(resil.control)
plot(simres) # Very bad!

# Drop non significant interactions (year as random effect)
resil.control.ny <- lmer(log10(resilience) ~ scale(richness)*spei6_category + scale(berger_parker)*spei6_category + (1|year) + (1|site/uniqueid), data = plot_ece_control)

# Add site
resil.control.rich.site <- lmer(log10(resilience) ~ scale(richness)*spei6_category + site + (1|year) + (1|site/uniqueid), data = plot_ece_control) # failed to converge

# Drop interaction, keep site
resil.control.rich.site.add <- lmer(log10(resilience) ~ scale(richness) + spei6_category + site + (1|year) + (1|site/uniqueid), data = plot_ece_control) # singular fit
summary(resil.control.rich.site.add)
simres <- simulateResiduals(resil.control.rich.site.add)
plot(simres) # Very bad!

# Site and spei interaction
resil.control.site.int <- lmer(log10(resilience) ~ scale(richness) + spei6_category*site + (1|year) + (1|site/uniqueid), data = plot_ece_control) # failed to converge

# Drop richness
resil.control.site.int.nr <- lmer(log10(resilience) ~ spei6_category*site + (1|year) + (1|site/uniqueid), data = plot_ece_control) # failed to converge


AICctab(resil.control.rich, resil.control.rich.ny, resil.control.rich.site.add, resil.control.add.ny, resil.control, resil.control.int, resil.control.int.ny, resil.control.add, resil.control.ny)

# Plot best model
ggpredict(model = resil.control.rich.ny, terms = c("richness", "spei6_category"), back_transform = F) %>%
  plot(show_data = TRUE)



## All plots ----


# environmental model for RESISTANCE!
resil.env = lm(log10(resilience) ~ nutrients_simple*fire_frequency_cat , data = plot_ece_rm_na)
summary(resil.env )
plot_ece_rm_na_wresidresil  = plot_ece_rm_na
plot_ece_rm_na_wresidresil$residresil = resid(resil.env)
hist(plot_ece_rm_na_wresidresil$residresil)

# run model on resistance residuals <- from environment 
resil.lm.hyp.resid <- lmer(residresil~ scale(richness)*spei6_category + scale(berger_parker)*spei6_category + (1|site) + (1|site:plot) + (1|year), data = plot_ece_rm_na_wresidresil)
summary(resil.lm.hyp.resid)


plot_ece_rm_na_wresidresil %>% ggplot(aes( x= richness, y=residresil, color=spei6_category)) +
  geom_point() + geom_smooth(method = "lm")



# ONE WAY TO DO THIS BETTER: SPLIT UP BY SPEI6 Category!!!
# way easier than having huge messy interaction !!!


# dry yearss 
resil.env.dw_dry = lmer(log10(resilience) ~ scale(richness)*nutrients_simple*fire_frequency_cat +
                           scale(berger_parker)*nutrients_simple*fire_frequency_cat +(1|site:plot) + (1|year), data = plot_ece_rm_na_dry, na.action = "na.fail")
summary(resil.env.dw_dry)
# do emtrends here?

plot_ece_rm_na_dry %>%  ggplot(aes (x=richness, y = log10(resilience), color=nutrients_simple) ) +
  geom_point() + geom_smooth(method = "lm") + 
  facet_grid(~fire_frequency_cat) # this is not sig

plot_ece_rm_na_dry %>%  ggplot(aes (x=berger_parker, y = log10(resilience), color=nutrients_simple) ) +
  geom_point() + geom_smooth(method = "lm")  # this is marginal but looks exactly the same




# wet years
resil.env.dw_wet = lmer(log10(resilience) ~ scale(richness)*nutrients_simple*fire_frequency_cat +
                           scale(berger_parker)*nutrients_simple*fire_frequency_cat +(1|site:plot) + (1|year), data = plot_ece_rm_na_wet, na.action = "na.fail")
summary(resil.env.dw_wet) # nothing is sig.








# Model 1
dim(plot_ece)
dim(plot_ece_rm_na)
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



