# TITLE:        LTER Grassland Rock: Create core analyses 
# AUTHORS:      Ashley Darst, Joshua Ajowele
# COLLABORATORS:  
# DATA INPUT:   Data imported as csv files from shared Google drive L2 folder
# DATA OUTPUT:  Core analyses
# PROJECT:      LTER Grassland Rock
# DATE:         October 2024 , last updated: April 2025

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
library(car)
library(GGally)
library(performance)
library(sjPlot)
library(kableExtra)

# Set working directory 
L0_dir <- Sys.getenv("L0DIR")
L1_dir <- Sys.getenv("L1DIR")
L2_dir <- Sys.getenv("L2DIR")
list.files(L2_dir)

# Read in CSV files
plot <- read.csv(file.path(L2_dir, "plot_metrics_SPEI_diversity_L2.csv"))
#ece <- read.csv(file.path(L2_dir, "ece_resist_resil_L2.csv"))#based on SPEI6
meta <- read.csv(file.path(L1_dir, "metadata_L1.csv"))
#ece_9<-read.csv(file.path(L2_dir, "ece_resist_resil_spei9_L2.csv"))#calculated based on SPEI9
ece_9_norm<-read.csv(file.path(L2_dir, "ece_resist_resil_spei9_norm_L2.csv"))#calculated based on SPEI9 without including moderate events as part of the normal events 

names(plot)

# Only keep distinct rows in ece
#ece <- distinct(ece)
#ece_9<-distinct(ece_9)
ece_9_norm<-distinct(ece_9_norm)

# Change ex_year to year
#ece <- ece %>%
  #rename(year = ex_year)
#ece_9<-ece_9%>%
  #rename(year = ex_year)
ece_9_norm<-ece_9_norm%>%
  rename(year = ex_year,
         resistance_n=resistance,
         resilience_n=resilience)

# Merge plot with resistance and resilience
#plot_ece <- left_join(plot, ece)
#plot_ece_9<-left_join(plot, ece_9)
plot_ece_9<-left_join(plot, ece_9_norm)

# Make column with categories for high and low dominance
#plot_ece <- plot_ece %>%
  #mutate(dom_category = case_when(
    #Berger_Parker >= 0.5 ~ "high",
    #Berger_Parker < 0.5 ~ "low"
  #))

# Standardize column names
#plot_ece <- clean_names(plot_ece)
plot_ece_9 <- clean_names(plot_ece_9)


# Convert Inf values to NA for resilience for some KBS 2015 plots
# plot_ece$resilience[plot_ece$resilience == Inf] <- NA
# plot_ece_9$resilience[plot_ece_9$resilience == Inf] <- NA

# Add experiment column_SPEI6####
# plot_ece$experiment <- sub("nutnet.*", "nutnet", plot_ece$higher_order_organization)
# plot_ece$experiment <- sub("glbrc_scaleup.*", "glbrc_scaleup", plot_ece$experiment)
# plot_ece$experiment <- sub("glbrc_G10.*", "glbrc", plot_ece$experiment)
# plot_ece$experiment <- sub("glbrc_G9.*", "glbrc", plot_ece$experiment)
# plot_ece$experiment <- sub("mcse.*", "mcse", plot_ece$experiment)
# plot_ece$experiment <- sub("microplots.*", "microplots", plot_ece$experiment)
# plot_ece$experiment <- sub("Experiment 1.*", "Experiment 1", plot_ece$experiment)
# plot_ece$experiment <- sub("001d_A_fl", "001d_fl", plot_ece$experiment)
# plot_ece$experiment <- sub("001d_B_fl", "001d_fl", plot_ece$experiment)
# plot_ece$experiment <- sub("001d_C_fl", "001d_fl", plot_ece$experiment)
# plot_ece$experiment <- sub("001d_D_fl", "001d_fl", plot_ece$experiment)
# plot_ece$experiment <- sub("001d_A_tu", "001d_tu", plot_ece$experiment)
# plot_ece$experiment <- sub("001d_B_tu", "001d_tu", plot_ece$experiment)
# plot_ece$experiment <- sub("001d_C_tu", "001d_tu", plot_ece$experiment)
# plot_ece$experiment <- sub("001d_D_tu", "001d_tu", plot_ece$experiment)
# plot_ece$experiment <- sub("004a_A_fl", "004a_fl", plot_ece$experiment)
# plot_ece$experiment <- sub("004a_B_fl", "004a_fl", plot_ece$experiment)
# plot_ece$experiment <- sub("004a_C_fl", "004a_fl", plot_ece$experiment)
# plot_ece$experiment <- sub("004a_D_fl", "004a_fl", plot_ece$experiment)
# plot_ece$experiment <- sub("004a_A_tu", "004a_tu", plot_ece$experiment)
# plot_ece$experiment <- sub("004a_B_tu", "004a_tu", plot_ece$experiment)
# plot_ece$experiment <- sub("004a_C_tu", "004a_tu", plot_ece$experiment)
# plot_ece$experiment <- sub("004a_D_tu", "004a_tu", plot_ece$experiment)
# plot_ece$experiment <- sub("004b_A_fl", "004b_fl", plot_ece$experiment)
# plot_ece$experiment <- sub("004b_B_fl", "004b_fl", plot_ece$experiment)
# plot_ece$experiment <- sub("004b_C_fl", "004b_fl", plot_ece$experiment)
# plot_ece$experiment <- sub("004b_D_fl", "004b_fl", plot_ece$experiment)
# plot_ece$experiment <- sub("004b_A_tu", "004b_tu", plot_ece$experiment)
# plot_ece$experiment <- sub("004b_B_tu", "004b_tu", plot_ece$experiment)
# plot_ece$experiment <- sub("004b_C_tu", "004b_tu", plot_ece$experiment)
# plot_ece$experiment <- sub("004b_D_tu", "004b_tu", plot_ece$experiment)
# plot_ece$experiment <- sub("002d.*", "002d", plot_ece$experiment)
# plot_ece$experiment <- sub("Experiment 54.*", "Experiment 54", plot_ece$experiment)
# plot_ece$experiment <- sub("KNZ_WAT01.*", "KNZ_WAT01", plot_ece$experiment)
# plot_ece$experiment <- sub("002c.*", "002c", plot_ece$experiment)
# plot_ece$experiment <- sub("e061.*", "e061", plot_ece$experiment)
# plot_ece$experiment <- sub("e247.*", "e247", plot_ece$experiment)
# plot_ece$experiment <- sub("e245.*", "e245", plot_ece$experiment)
# plot_ece$experiment[plot_ece$experiment == "A"] <- "NGE"
# plot_ece$experiment[plot_ece$experiment == "B"] <- "NGE"
# plot_ece$experiment[plot_ece$experiment == "C"] <- "NGE"
# plot_ece$experiment[plot_ece$experiment == "D"] <- "NGE"
# plot_ece$experiment[plot_ece$experiment == "E"] <- "NGE"
# plot_ece$experiment[plot_ece$experiment == "F"] <- "NGE"

plot(plot_ece_9$prior_year_rich, plot_ece_9$richness)
cor.test(plot_ece_9$prior_year_rich, plot_ece_9$richness)

#define experiment column-SPEI9####
plot_ece_9$experiment <- sub("nutnet.*", "nutnet", plot_ece_9$higher_order_organization)
plot_ece_9$experiment <- sub("glbrc_scaleup.*", "glbrc_scaleup", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("glbrc_G10.*", "glbrc", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("glbrc_G9.*", "glbrc", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("mcse.*", "mcse", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("microplots.*", "microplots", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("Experiment 1.*", "Experiment 1", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("001d_A_fl", "001d_fl", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("001d_B_fl", "001d_fl", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("001d_C_fl", "001d_fl", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("001d_D_fl", "001d_fl", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("001d_A_tu", "001d_tu", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("001d_B_tu", "001d_tu", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("001d_C_tu", "001d_tu", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("001d_D_tu", "001d_tu", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("004a_A_fl", "004a_fl", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("004a_B_fl", "004a_fl", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("004a_C_fl", "004a_fl", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("004a_D_fl", "004a_fl", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("004a_A_tu", "004a_tu", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("004a_B_tu", "004a_tu", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("004a_C_tu", "004a_tu", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("004a_D_tu", "004a_tu", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("004b_A_fl", "004b_fl", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("004b_B_fl", "004b_fl", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("004b_C_fl", "004b_fl", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("004b_D_fl", "004b_fl", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("004b_A_tu", "004b_tu", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("004b_B_tu", "004b_tu", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("004b_C_tu", "004b_tu", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("004b_D_tu", "004b_tu", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("002d.*", "002d", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("Experiment 54.*", "Experiment 54", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("KNZ_WAT01.*", "KNZ_WAT01", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("002c.*", "002c", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("e061.*", "e061", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("e247.*", "e247", plot_ece_9$experiment)
plot_ece_9$experiment <- sub("e245.*", "e245", plot_ece_9$experiment)
plot_ece_9$experiment[plot_ece_9$experiment == "A"] <- "NGE"
plot_ece_9$experiment[plot_ece_9$experiment == "B"] <- "NGE"
plot_ece_9$experiment[plot_ece_9$experiment == "C"] <- "NGE"
plot_ece_9$experiment[plot_ece_9$experiment == "D"] <- "NGE"
plot_ece_9$experiment[plot_ece_9$experiment == "E"] <- "NGE"
plot_ece_9$experiment[plot_ece_9$experiment == "F"] <- "NGE"

# taake out 
# plot_ece_prior = plot_ece %>% select(-c(dominant_relative_abund_zero,richness,evar))
# plot_ece_prior = plot_ece_prior %>%  rename(dominant_relative_abund_zero = prior_year_dom)
# plot_ece_prior = plot_ece_prior %>%  rename(evar = prior_year_evar)
# plot_ece_prior = plot_ece_prior %>%  rename(richness =prior_year_rich)

plot_ece_9_prior = plot_ece_9 %>% select(-c(dominant_relative_abund_zero,richness,evar))
plot_ece_9_prior = plot_ece_9_prior %>%  rename(dominant_relative_abund_zero = prior_year_dom)
plot_ece_9_prior = plot_ece_9_prior %>%  rename(evar = prior_year_evar)
plot_ece_9_prior = plot_ece_9_prior %>%  rename(richness =prior_year_rich)


# Merge with metadata

# plot_ece_meta <- left_join(plot_ece_prior, meta)
# unique(plot_ece_meta$experiment)
plot_ece_9_meta <- left_join(plot_ece_9_prior, meta)
unique(plot_ece_9_meta$experiment)


# Remove NAs for non-extreme years
# plot_ece_rm_na <- plot_ece_meta %>%
#   drop_na(resistance)
plot_ece_9_rm_na <- plot_ece_9_meta %>%
  drop_na(resistance_n)
#checking each experiment/study
# plot_filter<-plot_ece_rm_na%>%
#   filter(experiment=="004b_fl")

# Make year a factor
#str(plot_ece_rm_na)
str(plot_ece_9_rm_na)
# plot_ece_rm_na$year <- as.factor(plot_ece_rm_na$year)
# plot_ece_rm_na$measurement_scale_cover <- as.factor(plot_ece_rm_na$measurement_scale_cover)
plot_ece_9_rm_na$year <- as.factor(plot_ece_9_rm_na$year)
plot_ece_9_rm_na$measurement_scale_cover <- as.factor(plot_ece_9_rm_na$measurement_scale_cover)


#subset to have control and nitrogen treated plots
plot_ece_9_cn <- plot_ece_9_rm_na %>%
  filter(nutrients_added == "NPK+" | nutrients_added == "N" | nutrients_added == "NP" | nutrients_added == "NPK" | nutrients_added == "NK" | nutrients_added == "NPK+Fence" | nutrients_added == "none" | nutrients_added == "no_fertilizer")%>%
  filter(disturbance!="disturbed")%>%
  filter(treatment != "irrigated")%>%
  filter(treatment != "altered")%>%
  filter(grazing!="grazed")%>%
  filter(treatment!="insecticide")%>%#remove plots with additional treatment
  #add column with recoded nitrogen treatment
  mutate(nitrogen = dplyr::recode(nutrients_added, "NPK+" = "N", "NP" = "N", "NPK" = "N", "NK" = "N", "NPK+Fence" = "N", "none" = "no_fertilizer"),
         nutrients_grouped = dplyr::recode(nutrients_added, "NPK+" = "Nplus", "NP" = "Nplus", "NPK" = "Nplus", "NK" = "Nplus", "NPK+Fence" = "Nplus", "none" = "no_fertilizer")
  ) %>%
  mutate(extreme_after=case_when(after_year_type=="Extreme wet"~"yes",
                                 after_year_type=="Extreme dry"~"yes",
                                 .default="no"),
         nitrogen=case_when(nitrogen=="no_fertilizer"~"ambient",
                            nitrogen=="N"~"nutrients"))#when the year after the extreme event is another extreme event

#check treatment groups
unique(plot_ece_9_cn$treatment)
unique(plot_ece_9_cn$nitrogen)
table(plot_ece_9_cn$nutrients_grouped)
table(plot_ece_9_rm_na$nutrients_added)

#remove resilience with extreme events after an extreme event
plot_ece_9_cn_rm <- plot_ece_9_cn %>%
  mutate(resilience_n = case_when(extreme_after == "yes" ~ NA,
                        extreme_after == "no"  ~ resilience_n))

##### SEM #####
library(piecewiseSEM, lavaan)
source( paste(getwd(), "/R/L2/SEM_custom_functions_L2.R", sep = ""))

names(plot_ece_9_cn_rm)

table ( plot_ece_9_cn_rm$nutrients_grouped )

# get list of experiments that manipulated N
manipulated_n_exp = plot_ece_meta %>% filter (!is.na(nitrogen_amount))  %>% distinct(experiment) %>% select(experiment)

# get subset file for SEM:
plot_ece_meta_selSEM_init = plot_ece_9_cn_rm %>% 
  select (site, experiment, uniqueid,spei9_category,spei9, prior_year_spei9, nutrients_grouped , evar,
          richness,dominant_relative_abund_zero,resistance_n,resilience_n, year) #%>%  # only rows we care about
  #filter (  nutrients_grouped != "N" )  # experiment %in%  manipulated_n_exp[,1] 

plot_ece_meta_selSEM = plot_ece_meta_selSEM_init[complete.cases(plot_ece_meta_selSEM_init ) ,]
plot_ece_meta_selSEM$spei9_category = as.factor(plot_ece_meta_selSEM$spei9_category)
plot_ece_meta_selSEM$spei9_abs = abs(plot_ece_meta_selSEM$spei9)

table(plot_ece_meta_selSEM$nutrients_grouped)
plot_ece_meta_selSEM$nut_dummy = case_when(plot_ece_meta_selSEM$nutrients_grouped =="no_fertilizer" ~ 0,
                                           plot_ece_meta_selSEM$nutrients_grouped =="Nplus" ~ 1,
                                           plot_ece_meta_selSEM$nutrients_grouped =="N" ~ 1)
dim(plot_ece_meta_selSEM)
unique(plot_ece_meta_selSEM$spei9_category)
table(plot_ece_meta_selSEM$nitrogen_amount)
unique(plot_ece_9_cn$spei9_category)
unique(plot_ece_9_cn$nitrogen_amount)

table (is.na( plot_ece_meta_selSEM ) )

names(plot_ece_meta_selSEM)

df <- plot_ece_meta_selSEM %>%
  mutate(#richness = scale(richness),
    # dominant_relative_abund_zero = scale(dominant_relative_abund_zero),
    # evar = scale(evar),
    # resistance = scale(resistance, center =F),
    # resilience = scale(resilience, center=F),
    log_resistance = log10(resistance_n),
    log_resilience = log10(resilience_n),
    # log_resistance = scale(resistance)),
    # log_resilience = scale(resilience),
    # spei9  = scale(case_when(spei9 < 0 ~ spei9*(-1),.default = spei9)), 
    site = factor(site),
    experiment = factor(experiment),
    uniqueid = factor(uniqueid),
    spei9_category = factor(spei9_category),
    year = factor(year)
  )

# SCALE!!
df[,c("richness", "dominant_relative_abund_zero", "evar", "spei9_abs", "prior_year_spei9",
      "resistance_n", "resilience_n","log_resistance", "log_resilience")] <- 
  scale(df[,c("richness", "dominant_relative_abund_zero", "evar", "spei9_abs", "prior_year_spei9",
              "resistance_n", "resilience_n","log_resistance", "log_resilience")])

#df[,c("richness", "dominant_relative_abund_zero", "evar", "spei9_abs")] <- 
#  scale(df[,c("richness", "dominant_relative_abund_zero", "evar", "spei9_abs")])

df$eventSite = paste(df$spei9_category, df$site)

names(df)
head(df)
plot(plot_ece_meta_selSEM$richness, df$richness)


df_dry <- df %>% filter(spei9_category %in% "Extreme dry")
df_wet <- df %>% filter(spei9_category %in% "Extreme wet")




##################################### #
# Full SEM (without legacy effect) ####
##################################### #


## piecewise SEM ####
#-------------------#

### multigroup ####
# this is our main model, but it doesn't yet include the legacy effect of past spei.
# Doesn't converge with year as random intercept, weird...
model1 <- psem(
  lmer(log_resistance ~ spei9_abs + richness + dominant_relative_abund_zero + evar + nut_dummy + (1|site/experiment/uniqueid),
       data = df ),
  lmer(log_resilience ~ spei9_abs + richness + dominant_relative_abund_zero + evar + nut_dummy + (1|site/experiment/uniqueid),
       data = df ) ,
  lmer(richness ~                                                                    nut_dummy + (1|site/experiment/uniqueid),
       data = df),
  lmer(dominant_relative_abund_zero ~                                                nut_dummy + (1|site/experiment/uniqueid),
       data = df), 
  lmer(evar ~                                                                        nut_dummy + (1|site/experiment/uniqueid),
       data = df), 
  
  richness %~~% evar, 
  richness %~~% dominant_relative_abund_zero, 
  evar %~~% dominant_relative_abund_zero, 
  log_resistance %~~% log_resilience,
  
  data = df
)

model1.summary <- summary(model1)
multigroup2(model1,  group = "spei9_category")
multigroup3(model1,  group = "spei9_category") # use this if multigroup2 doesn't work


model1_all_sum = summary(model1, intercepts = T)
model1_all_sum_paths = model1_all_sum$dTable %>% as.data.frame()
model1_all_sum$dTable$Independ.Claim
model1_all_sum_table = model1_all_sum$coefficients %>% as.data.frame()
model1_all_sum_table$model = "All_no_legacy"
model1_all_sum_paths$model = "All_no_legacy"

# Fix convergence error by removing site from random effects structure
model1.1 <- psem(
  lmer(log_resistance ~ spei9_abs + richness + dominant_relative_abund_zero + evar + nut_dummy + (1|experiment/uniqueid), data = df),
  lmer(log_resilience ~ spei9_abs + richness + dominant_relative_abund_zero + evar + nut_dummy + (1|experiment/uniqueid), data = df) ,
  lmer(richness ~ nut_dummy + (1|experiment/uniqueid), data = df),
  lmer(dominant_relative_abund_zero ~ nut_dummy + (1|experiment/uniqueid), data = df), 
  lmer(evar ~ nut_dummy + (1|experiment/uniqueid), data = df), 
  
  richness %~~% evar, 
  richness %~~% dominant_relative_abund_zero, 
  evar %~~% dominant_relative_abund_zero, 
  log_resistance %~~% log_resilience,
  
  data = df
)

lmer1 <- lmer(log_resistance ~ spei9_abs + richness + dominant_relative_abund_zero + evar + nut_dummy + (1|site/experiment/uniqueid), data = df )
lmer2 <- lmer(log_resilience ~ spei9_abs + richness + dominant_relative_abund_zero + evar + nut_dummy + (1|experiment), data = df )
lmer3 <- lmer(richness ~ nut_dummy + (1|site/experiment/uniqueid), data = df)
lmer4 <- lmer(dominant_relative_abund_zero ~ nut_dummy + (1|site/experiment/uniqueid), data = df)
lmer5 <- lmer(evar ~ nut_dummy + (1|site/experiment/uniqueid), data = df)

model1.1.summary <- summary(model1.1) 
plot(model1.1)
multigroup2(model1.1,  group = "spei9_category")
multigroup3(model1.1,  group = "spei9_category") # use this if multigroup2 doesn't work

model1.1_all_sum = summary(model1.1, intercepts = T)
model1.1_all_sum_paths = model1.1_all_sum$dTable %>% as.data.frame()
model1.1_all_sum$dTable$Independ.Claim
model1.1_all_sum_table = model1.1_all_sum$coefficients %>% as.data.frame()
model1.1_all_sum_table$model = "All_no_legacy"
model1.1_all_sum_paths$model = "All_no_legacy"




# Piecewise table for export
piecewise <- multigroup2(model1.1,  group = "spei9_category")

piecewise$group.coefs %>%
  bind_rows(.,.id = "EventType") %>%
  select(EventType, Response, Predictor, Estimate, Std.Error, Std.Estimate, DF, P.Value) %>%
  
  mutate(across(c(Estimate, Std.Error, Std.Estimate, DF, P.Value), as.numeric)) %>% 
  mutate(P_print = case_when(P.Value < 0.001 ~ "<0.001", .default = paste(round(P.Value, digits = 3))),
         Estimate_print = round(Estimate, digits = 4),
         Std.Error_print = round(Std.Error, digits = 4),
         Std.Estimate_print = round(Std.Estimate, digits = 4),
         DF_print = round(DF, digits = 4),
         stars = case_when(P.Value < 0.001 ~ "***",
                           P.Value > 0.001 & P.Value <0.01 ~ "**",
                           P.Value > 0.01  & P.Value < 0.05 ~ "*",
                           P.Value > 0.05 & P.Value <0.1 ~ ".",
                           .default = ""),
         Response = case_when(Response %in% "log_resistance"               ~ "Resistance",
                              Response %in% "log_resilience"               ~ "Resilience",
                              Response %in% "richness"                     ~ "Richness",
                              Response %in% "dominant_relative_abund_zero" ~ "Dominance",
                              Response %in% "evar"                         ~ "Evenness"),
         Predictor = case_when(Predictor %in% "spei9_abs"                    ~ "Event strength",
                               Predictor %in% "richness"                     ~ "Richness",
                               Predictor %in% "dominant_relative_abund_zero" ~ "Dominance",
                               Predictor %in% "evar"                         ~ "Evenness",
                               Predictor %in% "nut_dummy"                    ~ "Nutrient addition"),
         Response = factor(Response, levels = c("Resistance", "Resilience", "Richness", "Dominance", "Evenness")),
         Predictor = factor(Predictor, levels = c("Event strength", "Richness", "Dominance", "Evenness", "Nutrient addition")))  %>%
  
  mutate(P_print = cell_spec(P_print, bold = ifelse(P.Value < 0.05, TRUE, FALSE)) ) %>%
  mutate(P_print = gsub(P_print,
                        pattern = '<span style=" font-weight: bold; " >0</span>',
                        replacement = '<span style=" font-weight: bold;    " ><0.001</span>')) %>%
  select(EventType, Response, Predictor, Estimate_print, Std.Error_print, Std.Estimate_print, DF_print, P_print, stars) %>%
  pivot_wider(id_cols = c(Response, Predictor),
              values_from = c(Estimate_print, Std.Error_print, Std.Estimate_print, DF_print, P_print, stars),
              names_from = EventType) %>%
  select("Response", "Predictor",
         "Estimate_print_Extreme dry",  "Std.Error_print_Extreme dry", 
         "Std.Estimate_print_Extreme dry", "DF_print_Extreme dry",
         "P_print_Extreme dry", "stars_Extreme dry",
         
         "Estimate_print_Extreme wet", "Std.Error_print_Extreme wet", 
         "Std.Estimate_print_Extreme wet", "DF_print_Extreme wet",
         "P_print_Extreme wet",  "stars_Extreme wet") %>%
  arrange(Response, Predictor) %>%
  kbl(escape = F,
      caption = "Piecewise SEM",
      col.names = c("Response", "Predictor", rep(c("Est.", "Std. Err.", "Std. Est.", "DF", "P-value", " "), 2))) %>%
  add_header_above(c(" " = 2, "Extreme dry" = 6, "Extreme wet" = 6)) %>%
  kable_paper()






### dry ####
# this is the same model as model 1, but fitted to the dry subset of the data
# this is done to get the R2 and the covariances

model1_dry <- psem(
  lmer(log_resistance ~ spei9_abs + richness + dominant_relative_abund_zero + evar + nut_dummy + (1|site/experiment/uniqueid),
       data = df_dry ),
  lmer(log_resilience ~ spei9_abs + richness + dominant_relative_abund_zero + evar + nut_dummy + (1|site/experiment/uniqueid),
       data = df_dry) ,
  lmer(richness ~                                                                    nut_dummy + (1|site/experiment/uniqueid),
       data = df_dry),
  lmer(dominant_relative_abund_zero ~                                                nut_dummy + (1|site/experiment/uniqueid),
       data = df_dry), 
  lmer(evar ~                                                                        nut_dummy + (1|site/experiment/uniqueid),
       data = df_dry), 
  
  richness %~~% evar, 
  richness %~~% dominant_relative_abund_zero, 
  evar %~~% dominant_relative_abund_zero, 
  log_resistance %~~% log_resilience,
  
  data = df_dry
)

summary(model1_dry, intercepts = T)

model1_dry_sum = summary(model1_dry, intercepts = T)
model1_dry_sum_paths = model1_dry_sum$dTable %>% as.data.frame()
model1_dry_sum_table = model1_dry_sum$coefficients %>% as.data.frame()
model1_dry_sum_table$model = "Dry_no_legacy"
model1_dry_sum_paths$model = "Dry_no_legacy"

# Repeat but for model1.1
model1.1_dry <- psem(
  lmer(log_resistance ~ spei9_abs + richness + dominant_relative_abund_zero + evar + nut_dummy + (1|experiment/uniqueid),
       data = df_dry ),
  lmer(log_resilience ~ spei9_abs + richness + dominant_relative_abund_zero + evar + nut_dummy + (1|experiment/uniqueid),
       data = df_dry) ,
  lmer(richness ~       spei9_abs +                                                  nut_dummy + (1|experiment/uniqueid),
       data = df_dry),
  lmer(dominant_relative_abund_zero ~                                                nut_dummy + (1|experiment/uniqueid),
       data = df_dry), 
  lmer(evar ~           spei9_abs +                                                  nut_dummy + (1|experiment/uniqueid),
       data = df_dry), 
  
  richness %~~% evar, 
  richness %~~% dominant_relative_abund_zero, 
  evar %~~% dominant_relative_abund_zero, 
  log_resistance %~~% log_resilience,
  
  data = df_dry
)

summary(model1.1_dry, intercepts = T)

model1.1_dry_sum = summary(model1.1_dry, intercepts = T)
model1.1_dry_sum_paths = model1.1_dry_sum$dTable %>% as.data.frame()
model1.1_dry_sum_table = model1.1_dry_sum$coefficients %>% as.data.frame()
model1.1_dry_sum_table$model = "Dry_no_legacy"
model1.1_dry_sum_paths$model = "Dry_no_legacy"

### wet ####
# this is the same model as model 1, but fitted to the wet subset of the data
# this is done to get the R2 and the covariances

model1_wet <- psem(
  lmer(log_resistance ~ spei9_abs + richness + dominant_relative_abund_zero + evar + nut_dummy + (1|site/experiment/uniqueid),
       data = df_wet ),
  lmer(log_resilience ~ spei9_abs + richness + dominant_relative_abund_zero + evar + nut_dummy + (1|site/experiment/uniqueid),
       data = df_wet) ,
  lmer(richness ~                                                                    nut_dummy + (1|site/experiment/uniqueid),
       data = df_wet),
  lmer(dominant_relative_abund_zero ~                                                nut_dummy +  (1|site/experiment/uniqueid),
       data = df_wet), 
  lmer(evar ~                                                                        nut_dummy + (1|site/experiment/uniqueid),
       data = df_wet), 
  
  richness %~~% evar, 
  richness %~~% dominant_relative_abund_zero, 
  evar %~~% dominant_relative_abund_zero, 
  log_resistance %~~% log_resilience,
  
  data = df_wet
)

summary(model1_wet, intercepts = T)
model1_wet_sum = summary(model1_wet, intercepts = T)
model1_wet_sum_paths = model1_wet_sum$dTable %>% as.data.frame()
model1_wet_sum_table = model1_wet_sum$coefficients %>% as.data.frame()
model1_wet_sum_table$model = "Wet_no_legacy"
model1_wet_sum_paths$model = "Wet_no_legacy"

# Repeat for model 1.1
model1.1_wet <- psem(
  lmer(log_resistance ~ spei9_abs + richness + dominant_relative_abund_zero + evar + nut_dummy + (1|experiment/uniqueid),
       data = df_wet ),
  lmer(log_resilience ~ spei9_abs + richness + dominant_relative_abund_zero + evar + nut_dummy + (1|experiment/uniqueid),
       data = df_wet) ,
  lmer(richness ~       spei9_abs +                                                  nut_dummy + (1|experiment/uniqueid),
       data = df_wet),
  lmer(dominant_relative_abund_zero ~                                                nut_dummy +  (1|experiment/uniqueid),
       data = df_wet), 
  lmer(evar ~          spei9_abs +                                                   nut_dummy + (1|experiment/uniqueid),
       data = df_wet), 
  
  richness %~~% evar, 
  richness %~~% dominant_relative_abund_zero, 
  evar %~~% dominant_relative_abund_zero, 
  log_resistance %~~% log_resilience,
  
  data = df_wet
)

summary(model1.1_wet, intercepts = T)
model1.1_wet_sum = summary(model1.1_wet, intercepts = T)
model1.1_wet_sum_paths = model1.1_wet_sum$dTable %>% as.data.frame()
model1.1_wet_sum_table = model1.1_wet_sum$coefficients %>% as.data.frame()
model1.1_wet_sum_table$model = "Wet_no_legacy"
model1.1_wet_sum_paths$model = "Wet_no_legacy"

## lavaan ####
#------------#
model1_lavaan <- '
# effects
log_resistance               ~ nut_dummy + spei9_abs + richness + dominant_relative_abund_zero + evar
log_resilience               ~ nut_dummy + spei9_abs + richness + dominant_relative_abund_zero + evar
richness                     ~ nut_dummy 
dominant_relative_abund_zero ~ nut_dummy 
evar                         ~ nut_dummy 

# correlations (by default constrained)
richness                     ~~ c("sr_ev", "sr_ev") * evar
evar                         ~~ c("ev_do", "ev_do") * dominant_relative_abund_zero
dominant_relative_abund_zero ~~ c("do_sr", "do_sr") * richness
log_resistance               ~~ c("rs_rl", "rs_rl") * log_resilience'


model1_lavaan_fit <- sem(model1_lavaan, data = df, group = "spei9_category")

fitmeasures(model1_lavaan_fit,  c("df", "AIC", "pvalue", "RMSEA", "CFI", "SRMR", "IFI")) 

test_constraints(fit = model1_lavaan_fit, groups_constrain = c(1,2)) %>%
  mutate(stars = case_when(pValue < 0.001 ~ "***",
                           pValue > 0.001 & pValue < 0.01 ~ "**",
                           pValue > 0.01 & pValue < 0.05 ~ "*",
                           pValue > 0.05 & pValue < 0.1 ~".",
                           .default =""))



model1_lavaan_c <- '
# effects
log_resistance               ~ c("rs_nu", "rs_nu") * nut_dummy + c("rs_sp", "rs_sp") * spei9_abs +                       richness + c("rs_do", "rs_do") * dominant_relative_abund_zero + c("rs_ev", "rs_ev") * evar
log_resilience               ~ c("rl_nu", "rl_nu") * nut_dummy +                       spei9_abs +                       richness + c("rl_do", "rl_do") * dominant_relative_abund_zero + c("rl_ev", "rl_ev") * evar # resilience ~ spei is marginally significant, but if I constrain it, there is a difference between the final constrained and fully unconstrained model
richness                     ~                       nut_dummy 
dominant_relative_abund_zero ~ c("do_nu", "do_nu") * nut_dummy 
evar                         ~                       nut_dummy 

# correlations (by default constrained)
richness                     ~~ c("sr_ev", "sr_ev") * evar
evar                         ~~ c("ev_do", "ev_do") * dominant_relative_abund_zero
dominant_relative_abund_zero ~~ c("do_sr", "do_sr") * richness
log_resistance               ~~ c("rs_rl", "rs_rl") * log_resilience' 

model1_lavaan_fit_c <- sem(model1_lavaan_c, data = df, group = "spei9_category")

anova(model1_lavaan_fit, model1_lavaan_fit_c) 

fitmeasures(model1_lavaan_fit_c,  c("df", "AIC", "pvalue", "RMSEA", "CFI", "SRMR", "IFI")) 
# good fit when
# p value > 0.05
# rmsea   < 0.05
# cfi     > 0.96 (0.90)
# srmr    < 0.08

summary(model1_lavaan_fit_c, rsquare=T)




# add survey design to lavaan model to account for sampling structure. Note: it is not as fancy as the nested structure

library(lavaan.survey)

site.survey <- svydesign(ids=~uniqueid, prob=~1, data=df)
# site.survey <- svydesign(ids=~experiment + uniqueid, nest = TRUE, data=df) # convergence issues like with piecewise if the nested structure is used

fit.survey   <- lavaan.survey(model1_lavaan_fit, survey.design=site.survey)
fit.survey_c <- lavaan.survey(model1_lavaan_fit_c, survey.design=site.survey)


anova(fit.survey, fit.survey_c)

fitmeasures(fit.survey_c,  c("df", "AIC", "pvalue", "RMSEA", "CFI", "SRMR", "IFI")) 
summary(fit.survey_c)

# constrained lavaan table for export
parameterEstimates(fit.survey_c) %>%
  filter(op %in% "~") %>%
  rename("Response" = "lhs",
         "Predictor" = "rhs",
         "P.Value" = "pvalue",
         "Estimate" = "est", 
         "Std.Error" = "se",
         "EventType" = "group",
         "Constrained" = "label") %>%
  mutate(EventType = case_when(EventType == 1 ~ "Extreme dry",
                               EventType == 2 ~ "Extreme wet"))%>% 
  mutate(P_print = case_when(P.Value < 0.001 ~ "<0.001", .default = paste(round(P.Value, digits = 3))),
         Estimate_print = round(Estimate, digits = 4),
         Std.Error_print = round(Std.Error, digits = 4),
         z_print = round(z, digits = 4),
         stars = case_when(P.Value < 0.001 ~ "***",
                           P.Value > 0.001 & P.Value <0.01 ~ "**",
                           P.Value > 0.01  & P.Value < 0.05 ~ "*",
                           P.Value > 0.05 & P.Value <0.1 ~ ".",
                           .default = ""),
         Response = case_when(Response %in% "log_resistance"               ~ "Resistance",
                              Response %in% "log_resilience"               ~ "Resilience",
                              Response %in% "richness"                     ~ "Richness",
                              Response %in% "dominant_relative_abund_zero" ~ "Dominance",
                              Response %in% "evar"                         ~ "Evenness"),
         Predictor = case_when(Predictor %in% "spei9_abs"                    ~ "Event strength",
                               Predictor %in% "richness"                     ~ "Richness",
                               Predictor %in% "dominant_relative_abund_zero" ~ "Dominance",
                               Predictor %in% "evar"                         ~ "Evenness",
                               Predictor %in% "nut_dummy"                    ~ "Nutrient addition"),
         Response = factor(Response, levels = c("Resistance", "Resilience", "Richness", "Dominance", "Evenness")),
         Predictor = factor(Predictor, levels = c("Event strength", "Richness", "Dominance", "Evenness", "Nutrient addition")),
         Constrained = case_when(Constrained %in% "" ~ "", .default = "c"))  %>%
  
  mutate(P_print = cell_spec(P_print, bold = ifelse(P.Value < 0.05, TRUE, FALSE)) ) %>%
  mutate(P_print = gsub(P_print,
                        pattern = '<span style=" font-weight: bold; " >0</span>',
                        replacement = '<span style=" font-weight: bold;    " ><0.001</span>')) %>%
  select(EventType, Response, Constrained, Predictor, Estimate_print, Std.Error_print, z_print, P_print, stars) %>%
  pivot_wider(id_cols = c(Response, Predictor, Constrained),
              values_from = c(Estimate_print, Std.Error_print, z_print, P_print, stars),
              names_from = EventType) %>%
  select("Response", "Predictor", "Constrained",
         "Estimate_print_Extreme dry",  "Std.Error_print_Extreme dry", 
         "z_print_Extreme dry",
         "P_print_Extreme dry", "stars_Extreme dry",
         
         "Estimate_print_Extreme wet", "Std.Error_print_Extreme wet", 
         "z_print_Extreme wet",
         "P_print_Extreme wet",  "stars_Extreme wet") %>%
  arrange(Response, Predictor) %>%
  kbl(escape = F,
      caption = "lavaan SEM constrained",
      col.names = c("Response", "Predictor", "", rep(c("Est.", "Std. Err.", "z-value", "P-value", " "), 2))) %>%
  add_header_above(c(" " = 3, "Extreme dry" = 5, "Extreme wet" = 5)) %>%
  kable_paper()


# unconstrained lavaan table for export
parameterEstimates(fit.survey) %>%
  filter(op %in% "~") %>%
  rename("Response" = "lhs",
         "Predictor" = "rhs",
         "P.Value" = "pvalue",
         "Estimate" = "est", 
         "Std.Error" = "se",
         "EventType" = "group",
         "Constrained" = "label") %>%
  mutate(EventType = case_when(EventType == 1 ~ "Extreme dry",
                               EventType == 2 ~ "Extreme wet"))%>% 
  mutate(P_print = case_when(P.Value < 0.001 ~ "<0.001", .default = paste(round(P.Value, digits = 3))),
         Estimate_print = round(Estimate, digits = 4),
         Std.Error_print = round(Std.Error, digits = 4),
         z_print = round(z, digits = 4),
         stars = case_when(P.Value < 0.001 ~ "***",
                           P.Value > 0.001 & P.Value <0.01 ~ "**",
                           P.Value > 0.01  & P.Value < 0.05 ~ "*",
                           P.Value > 0.05 & P.Value <0.1 ~ ".",
                           .default = ""),
         Response = case_when(Response %in% "log_resistance"               ~ "Resistance",
                              Response %in% "log_resilience"               ~ "Resilience",
                              Response %in% "richness"                     ~ "Richness",
                              Response %in% "dominant_relative_abund_zero" ~ "Dominance",
                              Response %in% "evar"                         ~ "Evenness"),
         Predictor = case_when(Predictor %in% "spei9_abs"                    ~ "Event strength",
                               Predictor %in% "richness"                     ~ "Richness",
                               Predictor %in% "dominant_relative_abund_zero" ~ "Dominance",
                               Predictor %in% "evar"                         ~ "Evenness",
                               Predictor %in% "nut_dummy"                    ~ "Nutrient addition"),
         Response = factor(Response, levels = c("Resistance", "Resilience", "Richness", "Dominance", "Evenness")),
         Predictor = factor(Predictor, levels = c("Event strength", "Richness", "Dominance", "Evenness", "Nutrient addition")),
         Constrained = case_when(Constrained %in% "" ~ "", .default = "c"))  %>%
  
  mutate(P_print = cell_spec(P_print, bold = ifelse(P.Value < 0.05, TRUE, FALSE)) ) %>%
  mutate(P_print = gsub(P_print,
                        pattern = '<span style=" font-weight: bold; " >0</span>',
                        replacement = '<span style=" font-weight: bold;    " ><0.001</span>')) %>%
  select(EventType, Response, Constrained, Predictor, Estimate_print, Std.Error_print, z_print, P_print, stars) %>%
  pivot_wider(id_cols = c(Response, Predictor, Constrained),
              values_from = c(Estimate_print, Std.Error_print, z_print, P_print, stars),
              names_from = EventType) %>%
  select("Response", "Predictor", "Constrained",
         "Estimate_print_Extreme dry",  "Std.Error_print_Extreme dry", 
         "z_print_Extreme dry",
         "P_print_Extreme dry", "stars_Extreme dry",
         
         "Estimate_print_Extreme wet", "Std.Error_print_Extreme wet", 
         "z_print_Extreme wet",
         "P_print_Extreme wet",  "stars_Extreme wet") %>%
  arrange(Response, Predictor) %>%
  kbl(escape = F,
      caption = "lavaan SEM unconstrained",
      col.names = c("Response", "Predictor", "", rep(c("Est.", "Std. Err.", "z-value", "P-value", " "), 2))) %>%
  add_header_above(c(" " = 3, "Extreme dry" = 5, "Extreme wet" = 5)) %>%
  kable_paper()


################################## #
# Full SEM (with legacy effect) ####
################################## #

## piecewise SEM ####
#-------------------#

### multigroup ####
# this is our second model, it includes the legacy effect of past spei.

model2 <- psem(
  lmer(log_resistance ~ spei9_abs + richness + dominant_relative_abund_zero + evar + nut_dummy + prior_year_spei9 + (1|site/experiment/uniqueid),
       data = df ),
  lmer(log_resilience ~ spei9_abs + richness + dominant_relative_abund_zero + evar + nut_dummy + prior_year_spei9 + (1|site/experiment/uniqueid),
       data = df ) ,
  lmer(richness ~                                                                    nut_dummy +                    (1|site/experiment/uniqueid),
       data = df),
  lmer(dominant_relative_abund_zero ~                                                nut_dummy +                    (1|site/experiment/uniqueid),
       data = df), 
  lmer(evar ~                                                                        nut_dummy +                    (1|site/experiment/uniqueid),
       data = df), 
  
  richness %~~% evar, 
  richness %~~% dominant_relative_abund_zero, 
  evar %~~% dominant_relative_abund_zero, 
  log_resistance %~~% log_resilience,
  
  data = df
)

summary(model2)
multigroup2(model2,  group = "spei9_category")

model2_all_sum = summary(model2, intercepts = T)
model2_all_sum_paths = model2_all_sum$dTable %>% as.data.frame()
model2_all_sum_table = model2_all_sum$coefficients %>% as.data.frame()
model2_all_sum_table$model = "All_legacy"
model2_all_sum_paths$model = "All_legacy"




### dry ####
# this is the same model as model 2, but fitted to the dry subset of the data
# this is done to get the R2 and the covariances

model2_dry <- psem(
  lmer(log_resistance ~ spei9_abs + richness + dominant_relative_abund_zero + evar + nut_dummy + prior_year_spei9 + (1|site/experiment/uniqueid),
       data = df_dry ),
  lmer(log_resilience ~ spei9_abs + richness + dominant_relative_abund_zero + evar + nut_dummy + prior_year_spei9 + (1|site/experiment/uniqueid),
       data = df_dry) ,
  lmer(richness ~                                                                    nut_dummy +                    (1|site/experiment/uniqueid),
       data = df_dry),
  lmer(dominant_relative_abund_zero ~                                                nut_dummy +                    (1|site/experiment/uniqueid),
       data = df_dry), 
  lmer(evar ~                                                                        nut_dummy +                    (1|site/experiment/uniqueid),
       data = df_dry), 
  
  richness %~~% evar, 
  richness %~~% dominant_relative_abund_zero, 
  evar %~~% dominant_relative_abund_zero, 
  log_resistance %~~% log_resilience,
  
  data = df_dry
)

summary(model2_dry, intercepts = T)

model2_dry_sum = summary(model2_dry, intercepts = T)
model2_dry_sum_paths = model2_dry_sum$dTable %>% as.data.frame()
model2_dry_sum_table = model2_dry_sum$coefficients %>% as.data.frame()
model2_dry_sum_table$model = "Dry_legacy"
model2_dry_sum_paths$model = "Dry_legacy"


### wet ####
# this is the same model as model 2, but fitted to the wet subset of the data
# this is done to get the R2 and the covariances

model2_wet <- psem(
  lmer(log_resistance ~ spei9_abs + richness + dominant_relative_abund_zero + evar + nut_dummy + prior_year_spei9 + (1|site/experiment/uniqueid),
       data = df_wet ),
  lmer(log_resilience ~ spei9_abs + richness + dominant_relative_abund_zero + evar + nut_dummy + prior_year_spei9 + (1|site/experiment/uniqueid),
       data = df_wet) ,
  lmer(richness ~                                                                    nut_dummy +                    (1|site/experiment/uniqueid),
       data = df_wet),
  lmer(dominant_relative_abund_zero ~                                                nut_dummy +                    (1|site/experiment/uniqueid),
       data = df_wet), 
  lmer(evar ~                                                                        nut_dummy +                    (1|site/experiment/uniqueid),
       data = df_wet), 
  
  richness %~~% evar, 
  richness %~~% dominant_relative_abund_zero, 
  evar %~~% dominant_relative_abund_zero, 
  log_resistance %~~% log_resilience,
  
  data = df_wet
)

summary(model2_wet, intercepts = T)

model2_wet_sum = summary(model2_wet, intercepts = T)
model2_wet_sum$coefficients
model2_wet_sum_paths = model2_wet_sum$dTable %>% as.data.frame()
model2_wet_sum_table = model2_wet_sum$coefficients %>% as.data.frame()
model2_wet_sum_table$model = "Wet_legacy"
model2_wet_sum_paths$model = "Wet_legacy"


##### combine output tables 
alloutput = rbind (model1_all_sum_table, model1_dry_sum_table, model1_wet_sum_table, 
                   model2_all_sum_table, model2_dry_sum_table, model2_wet_sum_table)

alloutput_order = alloutput  %>% relocate(model) 
str(alloutput_order)
names(alloutput_order)[10] <- "Stars"
alloutput_order$Std.Error = as.numeric(alloutput_order$Std.Error)
head(alloutput_order )

#View(alloutput_order)
str(alloutput_order)
select_if(alloutput_order, is.numeric )
alloutput_order_round = alloutput_order %>% dplyr::select(where(is.numeric))  %>% round(3)
meta_semtable = data.frame (alloutput_order[,1:3] )
stars_semtable = data.frame (alloutput_order[,10] )

final = cbind (meta_semtable, alloutput_order_round,stars_semtable) %>% as.data.frame()
names(final)
names(final)[10] <- "Stars"
getwd()
#write.csv(final, file.path(L2_dir, "./SEM_coefficients_output_L2.csv"), row.names=F)


allpath = rbind (model1_all_sum_paths, model1_dry_sum_paths, model1_wet_sum_paths, 
                   model2_all_sum_paths, model2_dry_sum_paths, model2_wet_sum_paths)
head(allpath)
allpath_order = allpath  %>% relocate(model) 
names(allpath_order)[7] <- "Stars"
allpath_final = allpath_order %>% mutate(across(where(is.numeric), round, 3))
#write.csv(allpath_final, file.path(L2_dir, "./SEM_coefficients_paths_L2.csv"), row.names=F)

## lavaan ####
#------------#
model2_lavaan <- '
# effects
log_resistance               ~ nut_dummy + spei9_abs + richness + dominant_relative_abund_zero + evar + prior_year_spei9
log_resilience               ~ nut_dummy + spei9_abs + richness + dominant_relative_abund_zero + evar + prior_year_spei9
richness                     ~ c("sr_nu")*nut_dummy 
dominant_relative_abund_zero ~ c("do_nu")*nut_dummy 
evar                         ~ c("ev_nu")*nut_dummy 

# correlations (by default constrained)
richness                     ~~ c("sr_ev", "sr_ev") * evar
evar                         ~~ c("ev_do", "ev_do") * dominant_relative_abund_zero
dominant_relative_abund_zero ~~ c("do_sr", "do_sr") * richness
log_resistance               ~~ log_resilience'


model2_lavaan_fit <- sem(model2_lavaan, data = df, group = "spei9_category")

test_constraints(fit = model2_lavaan_fit, groups_constrain = c(1,2)) %>%
  mutate(stars = case_when(pValue < 0.001 ~ "***",
                           pValue > 0.001 & pValue < 0.01 ~ "**",
                           pValue > 0.01 & pValue < 0.05 ~ "*",
                           pValue > 0.05 & pValue < 0.1 ~".",
                           .default =""))

model2_lavaan_c <- '
# effects
log_resistance               ~ c("rs_nu", "rs_nu") * nut_dummy + c("rs_sp", "rs_sp") * spei9_abs + c("rs_sr", "rs_sr") * richness + c("rs_do", "rs_do") * dominant_relative_abund_zero + c("rs_ev", "rs_ev") * evar + prior_year_spei9
log_resilience               ~                       nut_dummy +                       spei9_abs + c("rl_sr", "rl_sr") * richness + c("rl_do", "rl_do") * dominant_relative_abund_zero + c("rl_ev", "rl_ev") * evar + prior_year_spei9# resilience ~ spei is marginally significant, but if I constrain it, there is a difference between the final constrained and fully unconstrained model
richness                     ~ c("sr_nu")*nut_dummy 
dominant_relative_abund_zero ~ c("do_nu")*nut_dummy 
evar                         ~ c("ev_nu")*nut_dummy 

# correlations (by default constrained)
richness                     ~~ c("sr_ev", "sr_ev") * evar
evar                         ~~ c("ev_do", "ev_do") * dominant_relative_abund_zero
dominant_relative_abund_zero ~~ c("do_sr", "do_sr") * richness
log_resistance               ~~ log_resilience' 

model2_lavaan_fit_c <- sem(model2_lavaan_c, data = df, group = "spei9_category")

anova(model2_lavaan_fit, model2_lavaan_fit_c) 

fitmeasures(model2_lavaan_fit,  c("df", "AIC", "pvalue", "RMSEA", "CFI", "SRMR")) 
# good fit when
# p value > 0.05
# rmsea   < 0.05
# cfi     > 0.96 (0.90)
# srmr    < 0.08

summary(model2_lavaan_fit_c, rsquare=T)



########################################### #
# Simple SEM (no SPEI, no legacy effect) ####
########################################### #

## piecewise SEM ####
#-------------------#

### multigroup ####
# this is our main model, but it doesn't yet include the legacy effect of past spei.

model0 <- psem(
  lmer(log_resistance ~  richness + dominant_relative_abund_zero + evar + nut_dummy + (1|site/experiment/uniqueid),
       data = df ),
  lmer(log_resilience ~  richness + dominant_relative_abund_zero + evar + nut_dummy + (1|site/experiment/uniqueid),
       data = df ) ,
  lmer(richness ~                                                         nut_dummy + (1|site/experiment/uniqueid),
       data = df),
  lmer(dominant_relative_abund_zero ~                                     nut_dummy + (1|site/experiment/uniqueid),
       data = df), 
  lmer(evar ~                                                             nut_dummy + (1|site/experiment/uniqueid),
       data = df), 
  
  richness %~~% evar, 
  richness %~~% dominant_relative_abund_zero, 
  evar %~~% dominant_relative_abund_zero, 
  log_resistance %~~% log_resilience,
  
  data = df
)

summary(model0)
multigroup2(model0,  group = "spei9_category")



### dry ####
# this is the same model as model 0, but fitted to the dry subset of the data
# this is done to get the R2 and the covariances

model0_dry <- psem(
  lmer(log_resistance ~ richness + dominant_relative_abund_zero + evar + nut_dummy + (1|site/experiment/uniqueid),
       data = df_dry ),
  lmer(log_resilience ~ richness + dominant_relative_abund_zero + evar + nut_dummy + (1|site/experiment/uniqueid),
       data = df_dry) ,
  lmer(richness ~                                                        nut_dummy + (1|site/experiment/uniqueid),
       data = df_dry),
  lmer(dominant_relative_abund_zero ~                                    nut_dummy + (1|site/experiment/uniqueid),
       data = df_dry), 
  lmer(evar ~                                                            nut_dummy + (1|site/experiment/uniqueid),
       data = df_dry), 
  
  richness %~~% evar, 
  richness %~~% dominant_relative_abund_zero, 
  evar %~~% dominant_relative_abund_zero, 
  log_resistance %~~% log_resilience,
  
  data = df_dry
)

summary(model0_dry, intercepts = T)


### wet ####
# this is the same model as model 0, but fitted to the wet subset of the data
# this is done to get the R2 and the covariances

model0_wet <- psem(
  lmer(log_resistance ~ richness + dominant_relative_abund_zero + evar + nut_dummy + (1|site/experiment/uniqueid),
       data = df_wet ),
  lmer(log_resilience ~ richness + dominant_relative_abund_zero + evar + nut_dummy + (1|site/experiment/uniqueid),
       data = df_wet) ,
  lmer(richness ~                                                        nut_dummy + (1|site/experiment/uniqueid),
       data = df_wet),
  lmer(dominant_relative_abund_zero ~                                    nut_dummy +  (1|site/experiment/uniqueid),
       data = df_wet), 
  lmer(evar ~                                                            nut_dummy + (1|site/experiment/uniqueid),
       data = df_wet), 
  
  richness %~~% evar, 
  richness %~~% dominant_relative_abund_zero, 
  evar %~~% dominant_relative_abund_zero, 
  log_resistance %~~% log_resilience,
  
  data = df_wet
)

summary(model0_wet, intercepts = T)


## lavaan ####
#------------#
model0_lavaan <- '
# effects
log_resistance               ~ nut_dummy + richness + dominant_relative_abund_zero + evar
log_resilience               ~ nut_dummy + richness + dominant_relative_abund_zero + evar
richness                     ~ c("sr_nu")*nut_dummy 
dominant_relative_abund_zero ~ c("do_nu")*nut_dummy 
evar                         ~ c("ev_nu")*nut_dummy 

# correlations (by default constrained)
richness                     ~~ c("sr_ev", "sr_ev") * evar
evar                         ~~ c("ev_do", "ev_do") * dominant_relative_abund_zero
dominant_relative_abund_zero ~~ c("do_sr", "do_sr") * richness
log_resistance               ~~ log_resilience'


model0_lavaan_fit <- sem(model0_lavaan, data = df, group = "spei9_category")

test_constraints(fit = model0_lavaan_fit, groups_constrain = c(1,2)) %>%
  mutate(stars = case_when(pValue < 0.001 ~ "***",
                           pValue > 0.001 & pValue < 0.01 ~ "**",
                           pValue > 0.01 & pValue < 0.05 ~ "*",
                           pValue > 0.05 & pValue < 0.1 ~".",
                           .default =""))

model0_lavaan_c <- '
# effects
log_resistance               ~ c("rs_nu", "rs_nu") * nut_dummy + c("rs_sr", "rs_sr") * richness +                       dominant_relative_abund_zero + c("rs_ev", "rs_ev") * evar # resistance ~ dominance marginally significant -> but total result much better if unconstr.
log_resilience               ~                       nut_dummy +                       richness + c("rl_do", "rl_do") * dominant_relative_abund_zero + c("rl_ev", "rl_ev") * evar # resilience ~richness marginally significant -> but total result mucht better if unconstr.
richness                     ~ c("sr_nu")*nut_dummy 
dominant_relative_abund_zero ~ c("do_nu")*nut_dummy 
evar                         ~ c("ev_nu")*nut_dummy 

# correlations (by default constrained)
richness                     ~~ c("sr_ev", "sr_ev") * evar
evar                         ~~ c("ev_do", "ev_do") * dominant_relative_abund_zero
dominant_relative_abund_zero ~~ c("do_sr", "do_sr") * richness
log_resistance               ~~ log_resilience' 

model0_lavaan_fit_c <- sem(model0_lavaan_c, data = df, group = "spei9_category")

anova(model0_lavaan_fit, model0_lavaan_fit_c) # whopsie

fitmeasures(model0_lavaan_fit,  c("df", "AIC", "pvalue", "RMSEA", "CFI", "SRMR")) 
# good fit when
# p value > 0.05
# rmsea   < 0.05
# cfi     > 0.96 (0.90)
# srmr    < 0.08

summary(model0_lavaan_fit_c, rsquare=T)

################################################################################
################################################################################

# old (some overlap with before) #####
seraina_resist <- psem(
  lm(log_resistance ~ spei9_abs +richness + dominant_relative_abund_zero + evar + nut_dummy,
     data = df ),
  lm(richness ~ spei9_abs +nut_dummy ,
     data = df),
  lm(dominant_relative_abund_zero ~spei9_abs +nut_dummy ,
     data = df), 
  lm(evar ~spei9_abs +nut_dummy ,
       data = df), 
  lm(log_resilience~ spei9_abs + richness + dominant_relative_abund_zero + evar + nut_dummy,
     data = df) 
)

summary(seraina_resist)
cerror(richness %~~% evar, seraina_resist, df )
seraina_resist2 <- update(seraina_resist, richness %~~% evar)
seraina_resist3 <- update(seraina_resist2, richness %~~% dominant_relative_abund_zero)
seraina_resist4 <- update(seraina_resist3, evar %~~% dominant_relative_abund_zero)
seraina_resist5 <- update(seraina_resist4, log_resilience %~~%  log_resistance)

summary(seraina_resist5) # this works !!!! if 
plot(seraina_resist4)
multi = multigroup2(seraina_resist4, group = "spei9_category")
multi

df$eventSite
table(df$eventSite)
multigroup2(seraina_resist5, group = "eventSite")



seraina_nosem <- psem(
  lmer(log_resistance ~ richness + dominant_relative_abund_zero + evar +nut_dummy+(1|site/experiment/uniqueid),
       data = df ),
  lmer(richness ~ nut_dummy+(1|site/experiment/uniqueid),
       data = df),
  lmer(dominant_relative_abund_zero ~ nut_dummy +(1|site/experiment/uniqueid),
       data = df), 
  lmer(evar ~ nut_dummy +(1|site/experiment/uniqueid),
       data = df), 
  
  lmer(log_resilience ~richness + dominant_relative_abund_zero +evar + nut_dummy +  (1|site/experiment/uniqueid),
       data = df ) ,
  
  #richness %~~% evar, 
  richness %~~% dominant_relative_abund_zero, 
  #evar %~~% dominant_relative_abund_zero, 
  log_resistance %~~% log_resilience,
  
  data = df
)

summary(seraina_nosem )
multigroup2(seraina_nosem,  group = "spei9_category")
hist(df$spei9_abs)

# this is the model !!!!!!
seraina_all <- psem(
  lmer(log_resistance ~ spei9_abs +richness + dominant_relative_abund_zero + evar +nut_dummy+ +  (1|site/experiment/uniqueid),
       data = df ),
  lmer(richness ~ nut_dummy +  (1|site/experiment/uniqueid),
       data = df),
  lmer(dominant_relative_abund_zero ~ nut_dummy +   (1|site/experiment/uniqueid),
       data = df), 
  lmer(evar ~ nut_dummy +  (1|site/experiment/uniqueid),
       data = df), 
  
  lmer(log_resilience ~ spei9_abs + richness + dominant_relative_abund_zero +evar + nut_dummy +  (1|site/experiment/uniqueid),
       data = df ) ,
  
  richness %~~% evar, 
  richness %~~% dominant_relative_abund_zero, 
  evar %~~% dominant_relative_abund_zero, 
  log_resistance %~~% log_resilience,
  
  data = df
)
length(unique(df$experiment))


seraina_legacy <- psem(
  lmer(log_resistance ~ spei9_abs +richness + dominant_relative_abund_zero + evar +nut_dummy+ +  (1|site/experiment/uniqueid),
       data = df ),
  lmer(richness ~ nut_dummy +  (1|site/experiment/uniqueid),
       data = df),
  lmer(dominant_relative_abund_zero ~ nut_dummy +   (1|site/experiment/uniqueid),
       data = df), 
  lmer(evar ~ nut_dummy +  (1|site/experiment/uniqueid),
       data = df), 
  
  lmer(log_resilience ~ spei9_abs + richness + dominant_relative_abund_zero +evar + nut_dummy +  (1|site/experiment/uniqueid),
       data = df ) ,
  
  richness %~~% evar, 
  richness %~~% dominant_relative_abund_zero, 
  evar %~~% dominant_relative_abund_zero, 
  log_resistance %~~% log_resilience,
  
  data = df
)


 summary(seraina_all)

plot(seraina_all)
multigroup2(seraina_all, group = "spei9_category")
plot(multigroup$group.coefs)

df$eventSite
multigroup2(seraina_all, group = "eventSite")

#  separate models ....


seraina_RESISTANCE <- psem(
  lmer(log_resistance ~ richness + dominant_relative_abund_zero + evar +nut_dummy+(1|site/experiment/uniqueid),
       data = df ),
  lmer(richness ~ nut_dummy+(1|site/experiment/uniqueid),
       data = df),
  lmer(dominant_relative_abund_zero ~ nut_dummy +(1|site/experiment/uniqueid),
       data = df), 
  lmer(evar ~ nut_dummy +(1|site/experiment/uniqueid),
       data = df), 
  

  richness %~~% evar, 
  richness %~~% dominant_relative_abund_zero, 
  evar %~~% dominant_relative_abund_zero,
  
  data = df
)

summary(seraina_RESISTANCE )

############## dry 

df_dry = df %>% filter(spei9_category=="Extreme dry")

seraina_dry <- psem(
  lmer(log_resistance ~ spei9_abs +richness + dominant_relative_abund_zero + evar +nut_dummy+(1|site/experiment/uniqueid),
       data = df_dry ),
  lmer(richness ~ spei9_abs + nut_dummy+(1|site/experiment/uniqueid),
       data = df_dry),
  lmer(dominant_relative_abund_zero ~ spei9_abs +nut_dummy +(1|site/experiment/uniqueid),
       data = df_dry), 
  lmer(evar ~ spei9_abs + nut_dummy +(1|site/experiment/uniqueid),
       data = df_dry), 
  
  lmer(log_resilience ~ spei9_abs + richness + dominant_relative_abund_zero +evar + nut_dummy +  (1|site/experiment/uniqueid),
       data = df_dry ) ,
  
  richness %~~% evar,
  richness %~~% dominant_relative_abund_zero, 
 evar %~~% dominant_relative_abund_zero,
  log_resistance %~~% log_resilience,
  
  data = df_dry
)


summary(seraina_dry )


plot(seraina_dry )

# dropped 
seraina_dry_dropped <- psem(
  lmer(log_resistance ~ spei9 +richness +  evar +nut_dummy+(1|site/experiment/uniqueid),
       data = df_dry ),
  lmer(richness ~ spei9 + nut_dummy+(1|site/experiment/uniqueid),
       data = df_dry),
  lmer(evar ~ spei9 + nut_dummy +(1|site/experiment/uniqueid),
       data = df_dry), 
  
  lmer(log_resilience ~ spei9 + richness +evar + nut_dummy +  (1|site/experiment/uniqueid),
       data = df_dry ) ,
  
  richness %~~% evar, 
  log_resistance %~~% log_resilience,
  
  data = df_dry
)

summary(seraina_dry_dropped )



# no dominance 
seraina_dry2 <- psem(
  lmer(resistance ~ spei9 +richness +  evar +nut_dummy+(1|site/experiment/uniqueid),
       data = df_dry ),
  lmer(richness ~ spei9 + nut_dummy+(1|site/experiment/uniqueid),
       data = df_dry),
  lmer(evar ~ spei9 + nut_dummy +(1|site/experiment/uniqueid),
       data = df_dry), 
  
  lmer(resilience ~ spei9 + richness + evar + nut_dummy +  (1|site/experiment/uniqueid),
       data = df_dry ) ,
  
  richness %~~% evar, 
  resistance %~~% resilience,
  
  data = df_dry
)
summary(seraina_dry2)
plot(seraina_dry2)
plot(df_dry$richness, df_dry$evar)

# no spei
seraina_dry3 <- psem(
  lmer(log_resistance ~ richness + dominant_relative_abund_zero + evar +nut_dummy+(1|site/experiment/uniqueid),
       data = df_dry ),
  lmer(richness ~ nut_dummy+(1|site/experiment/uniqueid),
       data = df_dry),
  lmer(dominant_relative_abund_zero ~ nut_dummy +(1|site/experiment/uniqueid),
       data = df_dry), 
  lmer(evar ~ nut_dummy +(1|site/experiment/uniqueid),
       data = df_dry), 
  
  lmer(log_resilience ~richness + dominant_relative_abund_zero +evar + nut_dummy +  (1|site/experiment/uniqueid),
       data = df_dry ) ,
  
  #richness %~~% evar, 
  richness %~~% dominant_relative_abund_zero, 
  #evar %~~% dominant_relative_abund_zero, 
  log_resistance %~~% log_resilience,
  
  data = df_dry
)

summary(seraina_dry3 )
plot(seraina_dry3 )

############# wet
df_wet =df %>% filter(spei9_category=="Extreme wet")

seraina_wet <- psem(
  lmer(log_resistance ~ spei9_abs +richness + dominant_relative_abund_zero +evar + nut_dummy+(1|site/experiment/uniqueid),
       data = df_wet ),
  lmer(richness ~ spei9_abs +nut_dummy+(1|site/experiment/uniqueid),
       data = df_wet),
  lmer(dominant_relative_abund_zero ~spei9_abs +nut_dummy +(1|site/experiment/uniqueid),
       data = df_wet), 
  lmer(evar ~spei9_abs+nut_dummy +(1|site/experiment/uniqueid),
       data = df_wet), 
  
  lmer(log_resilience ~ spei9_abs+richness + dominant_relative_abund_zero + nut_dummy+ evar +(1|site/experiment/uniqueid),
       data = df_wet ) ,
  
  richness %~~% evar,
  richness %~~% dominant_relative_abund_zero, 
  evar %~~% dominant_relative_abund_zero,
  log_resistance %~~% log_resilience,
  
  data = df_wet
)


summary(seraina_wet )

plot(seraina_wet )


seraina_wet3 <- psem(
  lmer(resistance ~ richness + dominant_relative_abund_zero +evar + nut_dummy+(1|site/experiment/uniqueid),
       data = df_wet ),
  lmer(richness ~ nut_dummy+(1|site/experiment/uniqueid),
       data = df_wet),
  lmer(dominant_relative_abund_zero ~nut_dummy +(1|site/experiment/uniqueid),
       data = df_wet), 
  lmer(evar ~nut_dummy +(1|site/experiment/uniqueid),
       data = df_wet), 
  
  lmer(resilience ~ richness + dominant_relative_abund_zero + nut_dummy+ evar +(1|site/experiment/uniqueid),
       data = df_wet ) ,
  
  richness %~~% evar,
  richness %~~% dominant_relative_abund_zero, 
  evar %~~% dominant_relative_abund_zero,
  resistance %~~% resilience,
  
  data = df_wet
)

summary(seraina_wet3 )
plot(seraina_wet3 )
#multigroup2(list(A=seraina_wet, B = seraina_dry))


############## cdr
df$site
df_cdr = df %>% filter(site=="CDR")

seraina_cdr <- psem(
  lmer(log_resistance ~ spei9_abs +richness + dominant_relative_abund_zero + evar +nut_dummy+(1|experiment/uniqueid),
       data = df_cdr ),
  lmer(richness ~ spei9_abs + nut_dummy+(1|experiment/uniqueid),
       data = df_cdr),
  lmer(dominant_relative_abund_zero ~ spei9_abs +nut_dummy +(1|experiment/uniqueid),
       data = df_cdr), 
  lmer(evar ~ spei9_abs + nut_dummy +(1|experiment/uniqueid),
       data = df_cdr), 
  
  lmer(log_resilience ~ spei9_abs + richness + dominant_relative_abund_zero +evar + nut_dummy +  (1|experiment/uniqueid),
       data = df_cdr ) ,
  
  #richness %~~% evar, 
  richness %~~% dominant_relative_abund_zero, 
  # evar %~~% dominant_relative_abund_zero, 
  log_resistance %~~% log_resilience,
  
  data = df_cdr
)


summary(seraina_cdr )


multigroup2(seraina_cdr , group = "spei9_category")


############## knz
df$site
df_knz = df %>% filter(site=="knz")

seraina_knz <- psem(
  lmer(log_resistance ~ spei9_abs +richness + dominant_relative_abund_zero + evar +nut_dummy+(1|experiment/uniqueid),
       data = df_knz ),
  lmer(richness ~ spei9_abs + nut_dummy+(1|experiment/uniqueid),
       data = df_knz),
  lmer(dominant_relative_abund_zero ~ spei9_abs +nut_dummy +(1|experiment/uniqueid),
       data = df_knz), 
  lmer(evar ~ spei9_abs + nut_dummy +(1|experiment/uniqueid),
       data = df_knz), 
  
  lmer(log_resilience ~ spei9_abs + richness + dominant_relative_abund_zero +evar + nut_dummy +  (1|experiment/uniqueid),
       data = df_knz ) ,
  
  #richness %~~% evar, 
  richness %~~% dominant_relative_abund_zero, 
  # evar %~~% dominant_relative_abund_zero, 
  log_resistance %~~% log_resilience,
  
  data = df_knz
)


summary(seraina_knz )


multigroup2(seraina_knz , group = "spei9_category")












# matt's code below: 
# check distribution of data
plot_ece_control %>% 
  gather(key = Metric, value = Value, richness, 
         evar, dominant_relative_abund, dominant_relative_abund_zero,
         resistance) %>% 
  dplyr::select(Metric, Value) %>%
  ggplot(aes(x = Value))+
  geom_histogram()+
  facet_grid(.~Metric, scales = "free")

dim(plot_ece_control)

plot_ece_control_noNA <- plot_ece_control %>%
  dplyr::select(richness,
                evar, dominant_relative_abund, dominant_relative_abund_zero,
                resistance, fire_frequency) %>%
  drop_na() %>%
  data.frame

dim(plot_ece_control_noNA)


model.1 <- psem(
  lm(log10(resistance) ~ richness + dominant_relative_abund_zero + fire_frequency,
     data = plot_ece_control_noNA),
  lm(richness ~ fire_frequency,
     data = plot_ece_control_noNA),
  lm(dominant_relative_abund_zero ~ fire_frequency,
     data = plot_ece_control_noNA),
  data = plot_ece_control_noNA
)

summary(model.1, .progressBar = F) # suggests adding a path between dominant_relative_abund_zero and richness

plot(model.1)

plot_ece_control_noNA %>% 
  
  ggplot(., aes(x = richness, y = dominant_relative_abund_zero))+
  
  geom_point()+
  
  geom_smooth(method = "lm")+
  theme_bw()

summary(lm(dominant_relative_abund_zero ~ richness,
           data = plot_ece_control_noNA))

model.2 <- psem(
  lm(log10(resistance) ~ richness + dominant_relative_abund_zero + fire_frequency,
     data = plot_ece_control_noNA),
  lm(richness ~ fire_frequency,
     data = plot_ece_control_noNA),
  lm(dominant_relative_abund_zero ~ richness+ fire_frequency,
     data = plot_ece_control_noNA),
  data = plot_ece_control_noNA
)

summary(model.2, .progressBar = F) # suggests adding a path between dominant_relative_abund_zero and richness

plot(model.2)

model.3 <- psem(
  lm(log10(resistance) ~ richness + dominant_relative_abund_zero,
     data = plot_ece_control_noNA),
  lm(richness ~ fire_frequency,
     data = plot_ece_control_noNA),
  lm(dominant_relative_abund_zero ~ richness,
     data = plot_ece_control_noNA),
  data = plot_ece_control_noNA
)

summary(model.3, .progressBar = F) # suggests adding a path between dominant_relative_abund_zero and richness

plot(model.3)

model.4 <- psem(
  lm(log10(resistance) ~ richness + dominant_relative_abund_zero + fire_frequency,
     data = plot_ece_control_noNA),
  lm(richness ~ fire_frequency,
     data = plot_ece_control_noNA),
  lm(dominant_relative_abund_zero ~ fire_frequency + richness,
     data = plot_ece_control_noNA),
  data = plot_ece_control_noNA
)

summary(model.4, .progressBar = F) # suggests adding a path between dominant_relative_abund_zero and richness

plot(model.4)

anova(model.2, model.3)
