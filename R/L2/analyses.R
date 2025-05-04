# TITLE:        LTER Grassland Rock: Create core analyses 
# AUTHORS:      Ashley Darst, Joshua Ajowele
# COLLABORATORS:  
# DATA INPUT:   Data imported as csv files from shared Google drive L2 folder
# DATA OUTPUT:  Core analyses
# PROJECT:      LTER Grassland Rock
# DATE:         October 2024 , last updated: May 2025

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
library(patchwork)

# Set working directory 
L0_dir <- Sys.getenv("L0DIR")
L1_dir <- Sys.getenv("L1DIR")
L2_dir <- Sys.getenv("L2DIR")
list.files(L2_dir)

# Read in CSV files
plot <- read.csv(file.path(L2_dir, "plot_metrics_SPEI_diversity.csv"))
ece <- read.csv(file.path(L2_dir, "ece_resist_resil.csv"))#based on SPEI6
meta <- read.csv(file.path(L2_dir, "metadata.csv"))
ece_9<-read.csv(file.path(L2_dir, "ece_resist_resil_spei9.csv"))#calculated based on SPEI9

# Only keep distinct rows in ece
ece <- distinct(ece)
ece_9<-distinct(ece_9)

# Change ex_year to year
ece <- ece %>%
  rename(year = ex_year)
ece_9<-ece_9%>%
  rename(year = ex_year)

# Merge plot with resistance and resilience
plot_ece <- left_join(plot, ece)
plot_ece_9<-left_join(plot, ece_9)

# Make column with categories for high and low dominance
plot_ece <- plot_ece %>%
  mutate(dom_category = case_when(
    Berger_Parker >= 0.5 ~ "high",
    Berger_Parker < 0.5 ~ "low"
  ))

# Standardize column names
plot_ece <- clean_names(plot_ece)
plot_ece_9 <- clean_names(plot_ece_9)


# Convert Inf values to NA for resilience for some KBS 2015 plots
plot_ece$resilience[plot_ece$resilience == Inf] <- NA
plot_ece_9$resilience[plot_ece_9$resilience == Inf] <- NA

# Add experiment column_SPEI6####
plot_ece$experiment <- sub("nutnet.*", "nutnet", plot_ece$higher_order_organization)
plot_ece$experiment <- sub("glbrc_scaleup.*", "glbrc_scaleup", plot_ece$experiment)
plot_ece$experiment <- sub("glbrc_G10.*", "glbrc", plot_ece$experiment)
plot_ece$experiment <- sub("glbrc_G9.*", "glbrc", plot_ece$experiment)
plot_ece$experiment <- sub("mcse.*", "mcse", plot_ece$experiment)
plot_ece$experiment <- sub("microplots.*", "microplots", plot_ece$experiment)
plot_ece$experiment <- sub("Experiment 1.*", "Experiment 1", plot_ece$experiment)
plot_ece$experiment <- sub("001d_A_fl", "001d_fl", plot_ece$experiment)
plot_ece$experiment <- sub("001d_B_fl", "001d_fl", plot_ece$experiment)
plot_ece$experiment <- sub("001d_C_fl", "001d_fl", plot_ece$experiment)
plot_ece$experiment <- sub("001d_D_fl", "001d_fl", plot_ece$experiment)
plot_ece$experiment <- sub("001d_A_tu", "001d_tu", plot_ece$experiment)
plot_ece$experiment <- sub("001d_B_tu", "001d_tu", plot_ece$experiment)
plot_ece$experiment <- sub("001d_C_tu", "001d_tu", plot_ece$experiment)
plot_ece$experiment <- sub("001d_D_tu", "001d_tu", plot_ece$experiment)
plot_ece$experiment <- sub("004a_A_fl", "004a_fl", plot_ece$experiment)
plot_ece$experiment <- sub("004a_B_fl", "004a_fl", plot_ece$experiment)
plot_ece$experiment <- sub("004a_C_fl", "004a_fl", plot_ece$experiment)
plot_ece$experiment <- sub("004a_D_fl", "004a_fl", plot_ece$experiment)
plot_ece$experiment <- sub("004a_A_tu", "004a_tu", plot_ece$experiment)
plot_ece$experiment <- sub("004a_B_tu", "004a_tu", plot_ece$experiment)
plot_ece$experiment <- sub("004a_C_tu", "004a_tu", plot_ece$experiment)
plot_ece$experiment <- sub("004a_D_tu", "004a_tu", plot_ece$experiment)
plot_ece$experiment <- sub("004b_A_fl", "004b_fl", plot_ece$experiment)
plot_ece$experiment <- sub("004b_B_fl", "004b_fl", plot_ece$experiment)
plot_ece$experiment <- sub("004b_C_fl", "004b_fl", plot_ece$experiment)
plot_ece$experiment <- sub("004b_D_fl", "004b_fl", plot_ece$experiment)
plot_ece$experiment <- sub("004b_A_tu", "004b_tu", plot_ece$experiment)
plot_ece$experiment <- sub("004b_B_tu", "004b_tu", plot_ece$experiment)
plot_ece$experiment <- sub("004b_C_tu", "004b_tu", plot_ece$experiment)
plot_ece$experiment <- sub("004b_D_tu", "004b_tu", plot_ece$experiment)
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

# Merge with metadata
plot_ece_meta <- left_join(plot_ece, meta)
unique(plot_ece_meta$experiment)
plot_ece_9_meta <- left_join(plot_ece_9, meta)
unique(plot_ece_9_meta$experiment)


# Remove NAs for non-extreme years
plot_ece_rm_na <- plot_ece_meta %>%
  drop_na(resistance)
plot_ece_9_rm_na <- plot_ece_9_meta %>%
  drop_na(resistance)
#checking each experiment/study
plot_filter<-plot_ece_rm_na%>%
  filter(experiment=="004b_fl")

# Make year a factor
str(plot_ece_rm_na)
str(plot_ece_9_rm_na)
plot_ece_rm_na$year <- as.factor(plot_ece_rm_na$year)
plot_ece_rm_na$measurement_scale_cover <- as.factor(plot_ece_rm_na$measurement_scale_cover)
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
         )

#check treatment groups
unique(plot_ece_9_cn$treatment)
unique(plot_ece_9_cn$nitrogen)
table(plot_ece_9_cn$nutrients_grouped)
table(plot_ece_9_rm_na$nutrients_added)

#analysis of nutrient and climate event category on resistance####
resis_nitro<-lmer(log(resistance)~nitrogen*spei9_category+(1|site/experiment/uniqueid)
                  +(1|year), data=plot_ece_9_cn)
anova(resis_nitro)
summary(resis_nitro)
simres <- simulateResiduals(resis_nitro)
plot(simres)

resis_N_fig<-plot_ece_9_cn%>%
  ggplot(aes(x=spei9_category, y=log(resistance), col=nitrogen))+
  stat_summary(fun.data=mean_cl_boot, position=position_dodge(0.2))+
  theme_bw()

#including prior spei9 as covariate####
resis_prior<-lmer(log(resistance)~richness*dominant_relative_abund_zero+prior_year_spei9+nitrogen+richness:nitrogen+evar+
                 dominant_relative_abund_zero:nitrogen+richness:spei9_category+dominant_relative_abund_zero:spei9_category+
                 evar:nitrogen+evar:spei9_category+spei9_category+nitrogen:spei9_category+(1|site/experiment/uniqueid)
               +(1|year), data=plot_ece_9_cn_prior, REML=F)
summary(resis_prior)
anova(resis_prior)

resis_prior1<-lmer(log(resistance)~richness*dominant_relative_abund_zero+prior_year_spei9+nitrogen+richness:nitrogen+evar+
                    dominant_relative_abund_zero:nitrogen+richness:spei9_category+dominant_relative_abund_zero:spei9_category+
                    evar:nitrogen+evar:spei9_category+spei9_category+(1|site/experiment/uniqueid)
                  +(1|year), data=plot_ece_9_cn_prior, REML=F)
summary(resis_prior1)
anova(resis_prior1)

resis_prior2<-lmer(log(resistance)~richness*dominant_relative_abund_zero+prior_year_spei9+nitrogen+evar+
                     dominant_relative_abund_zero:nitrogen+richness:spei9_category+dominant_relative_abund_zero:spei9_category+
                     evar:nitrogen+evar:spei9_category+spei9_category+(1|site/experiment/uniqueid)
                   +(1|year), data=plot_ece_9_cn_prior, REML=F)
summary(resis_prior2)
anova(resis_prior2)

resis_prior3<-lmer(log(resistance)~richness*dominant_relative_abund_zero+prior_year_spei9+nitrogen+evar+
                     dominant_relative_abund_zero:nitrogen+richness:spei9_category+dominant_relative_abund_zero:spei9_category+
                     evar:nitrogen+spei9_category+(1|site/experiment/uniqueid)
                   +(1|year), data=plot_ece_9_cn_prior, REML=F)
summary(resis_prior3)
anova(resis_prior3)

resis_prior4<-lmer(log(resistance)~richness+dominant_relative_abund_zero+prior_year_spei9+nitrogen+evar+
                     dominant_relative_abund_zero:nitrogen+richness:spei9_category+dominant_relative_abund_zero:spei9_category+
                     evar:nitrogen+spei9_category+(1|site/experiment/uniqueid)
                   +(1|year), data=plot_ece_9_cn_prior, REML=F)
summary(resis_prior4)
anova(resis_prior4)

resis_prior5<-lmer(log(resistance)~richness+dominant_relative_abund_zero+prior_year_spei9+nitrogen+evar+
                     dominant_relative_abund_zero:nitrogen+richness:spei9_category+dominant_relative_abund_zero:spei9_category+
                     spei9_category+(1|site/experiment/uniqueid)
                   +(1|year), data=plot_ece_9_cn_prior, REML=F)
anova(resis_prior5)
resis_prior6<-lmer(log(resistance)~richness+dominant_relative_abund_zero+prior_year_spei9+nitrogen+evar+
                     richness:spei9_category+dominant_relative_abund_zero:spei9_category+
                     spei9_category+(1|site/experiment/uniqueid)
                   +(1|year), data=plot_ece_9_cn_prior, REML=F)
anova(resis_prior6)
resis_prior7<-lmer(log(resistance)~richness+dominant_relative_abund_zero+prior_year_spei9+nitrogen+evar+
                     dominant_relative_abund_zero:spei9_category+
                     spei9_category+(1|site/experiment/uniqueid)
                   +(1|year), data=plot_ece_9_cn_prior, REML=F)
anova(resis_prior7)#no major changes from model without prior year spei9

#Analysis of resistance with control and nitrogen####
#rename prior year community metrics and create a new dataframe
plot_ece_9_cn_prior<-plot_ece_9_cn%>%
  select(-richness, -evar, -dominant_relative_abund, -dominant_relative_abund_zero)%>%
  rename(richness=prior_year_rich,
         dominant_relative_abund_zero= prior_year_dom_zero,
         dominant_relative_abund=prior_year_dom,
         evar=prior_year_evar)


resis_cn<-lmer(log(resistance)~richness*dominant_relative_abund_zero+nitrogen+richness:nitrogen+evar+
                 dominant_relative_abund_zero:nitrogen+richness:spei9_category+dominant_relative_abund_zero:spei9_category+
                    evar:nitrogen+evar:spei9_category+spei9_category+nitrogen:spei9_category+(1|site/experiment/uniqueid)
                  +(1|year), data=plot_ece_9_cn_prior, REML=F)
summary(resis_cn)
anova(resis_cn)
#model update
resis_cn1<-lmer(log(resistance)~richness*dominant_relative_abund_zero+nitrogen+richness:nitrogen+evar+
                  dominant_relative_abund_zero:nitrogen+richness:spei9_category+dominant_relative_abund_zero:spei9_category+
                  evar:nitrogen+evar:spei9_category+spei9_category+(1|site/experiment/uniqueid)
                +(1|year), data=plot_ece_9_cn_prior, REML=F)
anova(resis_cn1)
#model update
resis_cn2<-lmer(log(resistance)~richness*dominant_relative_abund_zero+nitrogen+richness:nitrogen+evar+
                  dominant_relative_abund_zero:nitrogen+richness:spei9_category+dominant_relative_abund_zero:spei9_category+
                  evar:spei9_category+spei9_category+(1|site/experiment/uniqueid)
                +(1|year), data=plot_ece_9_cn_prior, REML=F)
anova(resis_cn2)
#model update
resis_cn3<-lmer(log(resistance)~richness*dominant_relative_abund_zero+nitrogen+richness:nitrogen+evar+
                  dominant_relative_abund_zero:spei9_category+richness:spei9_category+
                  evar:spei9_category+spei9_category+(1|site/experiment/uniqueid)
                +(1|year), data=plot_ece_9_cn_prior, REML=F)
anova(resis_cn3)
#model update
resis_cn4<-lmer(log(resistance)~richness*dominant_relative_abund_zero+nitrogen+evar+
                  richness:spei9_category+dominant_relative_abund_zero:spei9_category+
                  spei9_category+(1|site/experiment/uniqueid)
                +(1|year), data=plot_ece_9_cn_prior, REML=F)
anova(resis_cn4)
#model update
resis_cn5<-lmer(log(resistance)~richness+dominant_relative_abund_zero+nitrogen+evar+
                  richness:spei9_category+dominant_relative_abund_zero:spei9_category+
                  spei9_category+(1|site/experiment/uniqueid)
                +(1|year), data=plot_ece_9_cn_prior, REML=F)
anova(resis_cn5) 
#model update
resis_cn6<-lmer(log(resistance)~richness+dominant_relative_abund_zero+nitrogen+evar+
                  dominant_relative_abund_zero:spei9_category+
                  spei9_category+(1|site/experiment/uniqueid)
                +(1|year), data=plot_ece_9_cn_prior, REML=F)
anova(resis_cn6)
#model update
resis_cn7<-lmer(log(resistance)~richness+dominant_relative_abund_zero+nitrogen+
                  dominant_relative_abund_zero:spei9_category+
                  spei9_category+(1|site/experiment/uniqueid)
                +(1|year), data=plot_ece_9_cn_prior, REML=F)
anova(resis_cn7)
#model update
resis_cn8<-lmer(log(resistance)~richness+dominant_relative_abund_zero+
                  dominant_relative_abund_zero:spei9_category+
                  spei9_category+(1|site/experiment/uniqueid)
                +(1|year), data=plot_ece_9_cn_prior, REML=F)
anova(resis_cn8)

#model selection with likelihood ratio
anova(resis_cn,resis_cn1,resis_cn2,resis_cn3,resis_cn4,resis_cn5,resis_cn6,resis_cn7)
#resis_cn7 is the best model
#refit best model with REML
resis_cn9<-lmer(log(resistance)~richness+dominant_relative_abund_zero+nitrogen+evar+
                  dominant_relative_abund_zero:spei9_category+
                  spei9_category+(1|site/experiment/uniqueid)
                +(1|year), data=plot_ece_9_cn_prior)

anova(resis_cn9)
summary(resis_cn9)
simres <- simulateResiduals(resis_cn9)
plot(simres)

#plot best model for resistance####
res_rich_plot<-ggpredict(model = resis_cn9, terms = c("richness", "spei9_category"), back_transform = F ) %>%
  plot(show_data = TRUE)+
  labs(x="richness")
ggpredict(model = resis_cn9, terms = "evar", back_transform = F) %>%
  plot(show_data = TRUE)+
  labs(x="evenness")
ggpredict(model = resis_cn9, terms = "spei9_category", back_transform = F) %>%
  plot()+
  labs(x="Spei9 category")
resis_domin_plot<-ggpredict(model = resis_cn9, terms =  c("dominant_relative_abund_zero", "spei9_category"), back_transform = F) %>%
  plot(show_data = TRUE)+
  labs(x="Relative abundance of dominant species")
ggpredict(model = resis_cn9, terms =  c("dominant_relative_abund_zero", "spei9_category"), back_transform = F) %>%
  plot(show_data = TRUE)+
  labs(x="Relative abundance of dominant species")


#create effect size plot from the model
#modified from Seraina
resis_cn8_std <- update(resis_cn9, 
                             data = plot_ece_9_cn_prior %>% 
                               mutate(richness = scale(richness),
                                      evar = scale(evar),
                                      dominant_relative_abund_zero = scale(dominant_relative_abund_zero),
                               ))

resis_estim_plot<-coef(summary(resis_cn8_std)) %>%
  data.frame() %>%
  rownames_to_column("Variable") %>%
  rename("SE" = "Std..Error",
         "p" = "Pr...t..") %>%
  # order things neatly
  mutate(Variable = factor(Variable,
                           levels = c("(Intercept)",
                                      "richness",
                                      "dominant_relative_abund_zero",
                                      "evar",
                                      "nitrogenno_fertilizer",
                                      "spei9_categoryExtreme wet",
                                      "dominant_relative_abund_zero:spei9_categoryExtreme wet")),
         # cosmetics in the names
         Variable_labels = factor(case_when(Variable %in% "(Intercept)" ~ "intercept",
                                            Variable %in% "evar" ~ "evenness",
                                            Variable %in% "dominant_relative_abund_zero" ~ "dominant abundance",
                                            Variable %in% "nitrogenno_fertilizer" ~ "no fertilizer",
                                            Variable %in% "spei9_categoryExtreme wet" ~ "SPEI9: extreme wet",
                                            Variable %in% "dominant_relative_abund_zero:spei9_categoryExtreme wet" ~ "dominance x SPEI9: extreme wet",
                                            .default = Variable)),
         Variable_labels = fct_reorder(Variable_labels, as.numeric(Variable)),
         # get significances and their plotting location (perhaps necessary to play with the +/- offset)
         stars = case_when(p < 0.001 ~ "***",
                           p > 0.001 & p <0.01 ~ "**",
                           p > 0.01  & p < 0.05 ~ "*",
                           p > 0.05 & p <0.1 ~ ".",
                           .default = ""),
         stars_location = case_when(Estimate < 0 ~ Estimate - SE - 0.1,
                                    Estimate > 0 ~ Estimate + SE + 0.1)
  ) %>%
  filter(Variable!="(Intercept)")%>%
  
  ggplot(aes(y = Variable_labels, x = Estimate)) +
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  labs(title = "log(resistance)",
       x = "estimate \u00B1 se") +
  scale_y_discrete(limits = rev) +
  geom_vline(xintercept = 0) +
  geom_point()+
  geom_errorbarh(aes(xmin = Estimate-SE, xmax = Estimate+SE), height = 0.2) +
  geom_text(aes(x = stars_location, label = stars))

#nitrogen and climate event effect on resilience####
resil_nitro<-lmer(log(resilience)~nitrogen*spei9_category+(1|site/experiment/uniqueid)
                  +(1|year), data=plot_ece_9_cn)
anova(resil_nitro)
summary(resil_nitro)
simres <- simulateResiduals(resil_nitro)
plot(simres)

resil_N_fig<-plot_ece_9_cn%>%
  ggplot(aes(x=spei9_category, y=log(resilience), col=nitrogen))+
  stat_summary(fun.data=mean_cl_boot, position=position_dodge(0.2))+
  theme_bw()


#resilience including prior year spei9
resil_prior<-lmer(log(resilience)~richness*dominant_relative_abund_zero+prior_year_spei9+nitrogen+richness:nitrogen+evar+
                 dominant_relative_abund_zero:nitrogen+richness:spei9_category+dominant_relative_abund_zero:spei9_category+
                 evar:nitrogen+evar:spei9_category+spei9_category+nitrogen:spei9_category+(1|site/experiment/uniqueid)
               +(1|year), data=plot_ece_9_cn_prior, REML=F)
summary(resil_prior)
anova(resil_prior)
resil_prior1<-lmer(log(resilience)~richness*dominant_relative_abund_zero+prior_year_spei9+nitrogen+richness:nitrogen+evar+
                    dominant_relative_abund_zero:nitrogen+richness:spei9_category+dominant_relative_abund_zero:spei9_category+
                    evar:spei9_category+spei9_category+nitrogen:spei9_category+(1|site/experiment/uniqueid)
                  +(1|year), data=plot_ece_9_cn_prior, REML=F)
summary(resil_prior1)
resil_prior2<-lmer(log(resilience)~richness+dominant_relative_abund_zero+prior_year_spei9+nitrogen+richness:nitrogen+evar+
                     dominant_relative_abund_zero:nitrogen+richness:spei9_category+dominant_relative_abund_zero:spei9_category+
                     evar:spei9_category+spei9_category+nitrogen:spei9_category+(1|site/experiment/uniqueid)
                   +(1|year), data=plot_ece_9_cn_prior, REML=F)
summary(resil_prior2)
resil_prior3<-lmer(log(resilience)~richness+dominant_relative_abund_zero+prior_year_spei9+nitrogen+evar+
                     dominant_relative_abund_zero:nitrogen+richness:spei9_category+dominant_relative_abund_zero:spei9_category+
                     evar:spei9_category+spei9_category+nitrogen:spei9_category+(1|site/experiment/uniqueid)
                   +(1|year), data=plot_ece_9_cn_prior, REML=F)
summary(resil_prior3)
resil_prior4<-lmer(log(resilience)~richness+dominant_relative_abund_zero+prior_year_spei9+nitrogen+evar+
                     dominant_relative_abund_zero:nitrogen+dominant_relative_abund_zero:spei9_category+
                     evar:spei9_category+spei9_category+nitrogen:spei9_category+(1|site/experiment/uniqueid)
                   +(1|year), data=plot_ece_9_cn_prior, REML=F)
summary(resil_prior4)
resil_prior5<-lmer(log(resilience)~richness+dominant_relative_abund_zero+prior_year_spei9+nitrogen+evar+
                     dominant_relative_abund_zero:nitrogen+dominant_relative_abund_zero:spei9_category+
                     spei9_category+nitrogen:spei9_category+(1|site/experiment/uniqueid)
                   +(1|year), data=plot_ece_9_cn_prior, REML=F)
summary(resil_prior5)
resil_prior6<-lmer(log(resilience)~richness*dominant_relative_abund_zero*spei9_category+prior_year_spei9+nitrogen+evar+
                     dominant_relative_abund_zero:nitrogen+dominant_relative_abund_zero:spei9_category+
                     nitrogen:spei9_category+(1|site/experiment/uniqueid)
                   +(1|year), data=plot_ece_9_cn_prior, REML=F)
summary(resil_prior6)
resil_prior7<-lmer(log(resilience)~richness*dominant_relative_abund_zero*spei9_category+prior_year_spei9+nitrogen+evar+
                     dominant_relative_abund_zero:spei9_category+
                     nitrogen:spei9_category+(1|site/experiment/uniqueid)
                   +(1|year), data=plot_ece_9_cn_prior, REML=F)
summary(resil_prior7)
anova(resil_prior7)
resil_prior8<-lmer(log(resilience)~richness*dominant_relative_abund_zero*spei9_category+prior_year_spei9+nitrogen+evar+
                     dominant_relative_abund_zero:spei9_category+
                     (1|site/experiment/uniqueid)
                   +(1|year), data=plot_ece_9_cn_prior, REML=F)
summary(resil_prior8)
anova(resil_prior8)#may need to investigate further-three way interaction similar to model without prior year spei9
anova(resil_prior8,resil_prior7,resil_prior6,resil_prior5,resil_prior4,resil_prior3,resil_prior2,resil_prior1,resil_prior)

#Analysis of resilience with control and nitrogen####
resil_cn<-lmer(log(resilience)~richness*dominant_relative_abund_zero+nitrogen+richness:nitrogen+evar+
                 dominant_relative_abund_zero:nitrogen+richness:spei9_category+dominant_relative_abund_zero:spei9_category+
                 evar:nitrogen+evar:spei9_category+spei9_category+nitrogen:spei9_category+(1|site/experiment/uniqueid)
               +(1|year), data=plot_ece_9_cn_prior, REML=F)
summary(resil_cn)
anova(resil_cn)
#model update
resil_cn1<-lmer(log(resilience)~richness+dominant_relative_abund_zero+nitrogen+richness:nitrogen+evar+
                  dominant_relative_abund_zero:nitrogen+richness:spei9_category+dominant_relative_abund_zero:spei9_category+
                  evar:nitrogen+evar:spei9_category+spei9_category+nitrogen:spei9_category+(1|site/experiment/uniqueid)
                +(1|year), data=plot_ece_9_cn_prior, REML=F)
anova(resil_cn1)
#model update
resil_cn2<-lmer(log(resilience)~richness+dominant_relative_abund_zero+nitrogen+richness:nitrogen+evar+
                  dominant_relative_abund_zero:nitrogen+richness:spei9_category+dominant_relative_abund_zero:spei9_category+
                  +evar:spei9_category+spei9_category+nitrogen:spei9_category+(1|site/experiment/uniqueid)
                +(1|year), data=plot_ece_9_cn_prior, REML=F)
anova(resil_cn2)
#model update
resil_cn3<-lmer(log(resilience)~richness*dominant_relative_abund_zero+nitrogen+richness:nitrogen+evar+
                  dominant_relative_abund_zero:nitrogen+richness:spei9_category+dominant_relative_abund_zero:spei9_category+
                  +spei9_category+nitrogen:spei9_category+(1|site/experiment/uniqueid)
                +(1|year), data=plot_ece_9_cn_prior, REML=F)
anova(resil_cn3)
#model update
resil_cn4<-lmer(log(resilience)~richness*dominant_relative_abund_zero+nitrogen+evar+
                  richness:spei9_category+evar:spei9_category+spei9_category+(1|site/experiment/uniqueid)
                +(1|year), data=plot_ece_9_cn_prior, REML=F)
anova(resil_cn4)
summary(resil_cn4)
#model update
resil_cn5<-lmer(log(resilience)~richness*dominant_relative_abund_zero+nitrogen+evar+
                  richness:spei9_category+spei9_category+(1|site/experiment/uniqueid)
                +(1|year), data=plot_ece_9_cn_prior, REML=F)
anova(resil_cn5)
summary(resil_cn5)

#model update
resil_cn6<-lmer(log(resilience)~richness*dominant_relative_abund_zero*spei9_category+nitrogen+evar+
                  +(1|site/experiment/uniqueid)
                +(1|year), data=plot_ece_9_cn_prior, REML=F)
anova(resil_cn6)
summary(resil_cn6)
#model selection with likelihood ratio
anova(resil_cn6,resil_cn5,resil_cn4,resil_cn3,resil_cn2,resil_cn1,resil_cn)
#since the model output multiple two way interaction, we decided to include a meaningful three way interaction
#model update
resil_cn7<-lmer(log(resilience)~richness*dominant_relative_abund_zero*spei9_category+nitrogen+evar+
       (1|site/experiment/uniqueid)
     +(1|year), data=plot_ece_9_cn_prior, REML=F)
anova(resil_cn7)
#model selection
anova(resil_cn7, resil_cn6)

#refit with REML
resil_cn8<-lmer(log(resilience)~richness*dominant_relative_abund_zero*spei9_category+nitrogen+evar+
                  (1|site/experiment/uniqueid)
                +(1|year), data=plot_ece_9_cn_prior)

anova(resil_cn8)
summary(resil_cn8)
#model diagnostics
simres <- simulateResiduals(resil_cn8)
plot(simres)
#include resistance as a covariate in the model to test if things change
resil_cn_test<-lmer(log(resilience)~richness*dominant_relative_abund_zero*spei9_category+nitrogen+evar+
                  resistance+(1|site/experiment/uniqueid)
                +(1|year), data=plot_ece_9_cn_prior)
anova(resil_cn_test)#similar with or without resistance

#plot best model for reilience####
ggpredict(model = resil_cn8, terms = c("richness", "dominant_relative_abund_zero","spei9_category"), back_transform = F ) %>%
  plot(show_data = TRUE)+
  labs(x="richness")
resil_model_plot<-ggpredict(model = resil_cn8, terms = c("richness","spei9_category","nitrogen", "dominant_relative_abund_zero"), back_transform = F ) %>%
  plot(show_data = TRUE)+
  labs(x="richness")
ggpredict(model = resil_cn8, terms = c("richness", "spei9_category","nitrogen"), back_transform = F ) %>%
  plot(show_data = TRUE)+
  labs(x="richness")
ggpredict(model = resil_cn8, terms = c("dominant_relative_abund_zero","spei9_category","nitrogen"), back_transform = F ) %>%
  plot(show_data = TRUE)+
  labs(x="Relative abundance of dominant species")
ggpredict(model = resil_cn8, terms = "spei9_category", back_transform = F) %>%
  plot()+
  labs(x="Spei9 category")
ggpredict(model = resil_cn8, terms = "dominant_relative_abund_zero", back_transform = F) %>%
  plot(show_data = TRUE)+
  labs(x="Relative abundance of dominant species")

#combine resistance and resilience figures####
res_rich_plot + resis_domin_plot& plot_annotation(tag_levels = 'A')
resil_model_plot& plot_annotation(tag_levels = 'A')

#refit model with prior year before event####
resil_cn9<-lmer(log(resilience)~prior_year_rich*prior_year_dom_zero*spei9_category+
                       nitrogen+prior_year_evar+(1|site/experiment/uniqueid)
                     +(1|year), data=plot_ece_9_cn)
summary(resil_cn9)
anova(resil_cn9)
simres <- simulateResiduals(resil_cn9)
plot(simres)
check_model(resis_spei9_m4)
#create effect size plot from resilience model
resil_cn8_std <- update(resil_cn8, 
                        data = plot_ece_9_cn_prior %>% 
                          mutate(richness = scale(richness),
                                 evar = scale(evar),
                                 dominant_relative_abund_zero = scale(dominant_relative_abund_zero),
                          ))

resil_estim_plot<-coef(summary(resil_cn8_std)) %>%
  data.frame() %>%
  rownames_to_column("Variable") %>%
  rename("SE" = "Std..Error",
         "p" = "Pr...t..") %>%
  # order things neatly
  mutate(Variable = factor(Variable,
                           levels = c("(Intercept)",
                                      "richness",
                                      "dominant_relative_abund_zero",
                                      "evar",
                                      "nitrogenno_fertilizer",
                                      "spei9_categoryExtreme wet",
                                      "richness:spei9_categoryExtreme wet",
                                      "richness:dominant_relative_abund_zero",
                                      "dominant_relative_abund_zero:spei9_categoryExtreme wet",
                                      "richness:dominant_relative_abund_zero:spei9_categoryExtreme wet")),
         # cosmetics in the names
         Variable_labels = factor(case_when(Variable %in% "(Intercept)" ~ "intercept",
                                            Variable %in% "evar" ~ "evenness",
                                            Variable %in% "dominant_relative_abund_zero" ~ "dominant abundance",
                                            Variable %in% "nitrogenno_fertilizer" ~ "no fertilizer",
                                            Variable %in% "spei9_categoryExtreme wet" ~ "SPEI9: extreme wet",
                                            Variable %in% "richness:dominant_relative_abund_zero" ~ "richness x dominant abundance",
                                            Variable %in% "dominant_relative_abund_zero:spei9_categoryExtreme wet" ~ "dominant abundnace x SPEI9: extreme wet",
                                            Variable %in% "richness:spei9_categoryExtreme wet" ~ "richness x SPEI9: extreme wet",
                                            Variable %in% "richness:dominant_relative_abund_zero:spei9_categoryExtreme wet" ~ "richness x dominant abundance x SPEI9: extreme wet",
                                            .default = Variable)),
         Variable_labels = fct_reorder(Variable_labels, as.numeric(Variable)),
         # get significances and their plotting location (perhaps necessary to play with the +/- offset)
         stars = case_when(p < 0.001 ~ "***",
                           p > 0.001 & p <0.01 ~ "**",
                           p > 0.01  & p < 0.05 ~ "*",
                           p > 0.05 & p <0.1 ~ ".",
                           .default = ""),
         stars_location = case_when(Estimate < 0 ~ Estimate - SE - 0.1,
                                    Estimate > 0 ~ Estimate + SE + 0.1)
  ) %>%
  filter(Variable!="(Intercept)")%>%
  
  ggplot(aes(y = Variable_labels, x = Estimate)) +
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  labs(title = "log(resilience)",
       x = "estimate \u00B1 se") +
  scale_y_discrete(limits = rev) +
  geom_vline(xintercept = 0) +
  geom_point()+
  geom_errorbarh(aes(xmin = Estimate-SE, xmax = Estimate+SE), height = 0.2) +
  geom_text(aes(x = stars_location, label = stars))

#combine resistance and resilience effect size figure####
resis_estim_plot+resil_estim_plot & plot_annotation(tag_levels = 'A')

#runing model for just control using community metrics from year before event#### 
#filetr for only control plots
plot_ece_9_c_prior<-plot_ece_9_cn_prior%>%
  filter(nitrogen=="no_fertilizer")

#analysis of nutrient and climate event category on resistance####
resis_c<-lmer(log(resistance)~spei9_category+(1|site/experiment/uniqueid)
                  +(1|year), data=plot_ece_9_c_prior)
anova(resis_c)
summary(resis_c)
simres <- simulateResiduals(resis_c)
plot(simres)

resis_c_fig<-plot_ece_9_c_prior%>%
  ggplot(aes(x=spei9_category, y=log(resistance)))+
  stat_summary(fun.data=mean_cl_boot, position=position_dodge(0.2))+
  theme_bw()





# Subset to only have control plots####
plot_ece_control <- plot_ece_rm_na %>%
  filter(treatment == "control")
plot_ece_9_control <- plot_ece_9_cn %>%
  filter(nitrogen == "no_fertilizer")

#examine extreme event category
plot_wet<-plot_ece_control%>%
  filter(spei6_category=="Extreme wet")

# Look if measurment scale cover matters for resilience
plot_ece_control %>%
  ggplot(aes(x = measurement_scale_cover, y = log(resilience))) +
  geom_point()
# Look if measurment scale cover matters for richness
plot_ece_control %>%
  ggplot(aes(x = measurement_scale_cover, y = richness)) +
  geom_point()

# Try to scale richness and BP by experiment for control plots only
plot_ece_control <- plot_ece_control %>%
  group_by(experiment, year) %>%
  mutate(richness_scaled = c(scale(richness)),
         berger_parker_scaled = c(scale(berger_parker)),
         evar_scaled=c(scale(evar)),
         dominant_relative_abund_scaled=c(scale(dominant_relative_abund)),
         dominant_relative_abund_zero_scaled=c(scale(dominant_relative_abund_zero)))

# Look if measurment scale cover matters for richness after scaling
plot_ece_control %>%
  ggplot(aes(x = measurement_scale_cover, y = richness_scaled)) +
  geom_point()#no longer matters
rich_area<-lm(richness_scaled~measurement_scale_cover, data=plot_ece_control)
summary(rich_area) 
simres <- simulateResiduals(rich_area)
plot(simres)#looks good

#analysis using SPEI9####
#Resistance
#interaction based on hypothesis
resis_spei9<-lmer(log(resistance)~richness*dominant_relative_abund_zero+evar+
                              richness:spei9_category+dominant_relative_abund_zero:spei9_category+
                              evar:spei9_category+spei9_category+(1|site/experiment/uniqueid)
                            +(1|year), data=plot_ece_9_control, REML=F)
summary(resis_spei9)
anova(resis_spei9)
#without interaction of main predictors
resis_spei9_m1<-lmer(log(resistance)~richness+evar+dominant_relative_abund_zero+
                             richness:spei9_category+dominant_relative_abund_zero:spei9_category+
                             evar:spei9_category+spei9_category+(1|site/experiment/uniqueid)
                           +(1|year), data=plot_ece_9_control, REML=F)
summary(resis_spei9_m1)
anova(resis_spei9_m1)

#remove non-significant interaction
resis_spei9_m2<-lmer(log(resistance)~richness+evar+dominant_relative_abund_zero+
                              richness:spei9_category+
                              evar:spei9_category+spei9_category+(1|site/experiment/uniqueid)
                            +(1|year), data=plot_ece_9_control, REML=F)
anova(resis_spei9_m2)
summary(resis_spei9_m2)

#remove other non-significant interactions
resis_spei9_m3<-lmer(log(resistance)~richness+evar+dominant_relative_abund_zero+
                              richness:spei9_category+spei9_category+(1|site/experiment/uniqueid)
                            +(1|year), data=plot_ece_9_control, REML=F)
anova(resis_spei9_m3)
anova(resis_spei9_m3, resis_spei9_m2,resis_spei9_m1,resis_spei9)#model selection
#refit additive model with REML
resis_spei9_m4<-lmer(log(resistance)~richness+evar+dominant_relative_abund_zero+
                              richness:spei9_category+spei9_category+(1|site/experiment/uniqueid)
                            +(1|year), data=plot_ece_9_control)
summary(resis_spei9_m4)
anova(resis_spei9_m4)
simres <- simulateResiduals(resis_spei9_m4)
plot(simres)
check_model(resis_spei9_m4)

#plot best model
ggpredict(model = resis_spei9_m4, terms = c("richness", "spei9_category"), back_transform = F ) %>%
  plot(show_data = TRUE)+
  labs(x="richness")
ggpredict(model = resis_spei9_m4, terms = "richness", back_transform = F) %>%
  plot(show_data = TRUE)+
  labs(x="richness")
ggpredict(model = resis_spei9_m4, terms = "dominant_relative_abund_zero", back_transform = F) %>%
  plot(show_data = TRUE)+
  labs(x="Relative abundance of dominant species")

# Plot standardized effect sizes # This seems weird...
plot_model(resis_spei9_m4,
           type = "std",
           rm.terms = "measurement_scale_cover [0.3,0.4,1,10]",
           vline.color = "black",
           sort.est = TRUE,
           ci.lvl = 0.95
) + ylim(-0.05, 0.05) + theme_bw()


# Serainas version
# I think sjPlot does average across factor levels of other variables. I am not sure how to correctly do that manually
# note: this is only relevant for richness, as there is the reichness times SPEI interaction
resis_spei9_m4_std <- update(resis_spei9_m4, 
                             data = plot_ece_9_control %>% 
                               mutate(richness = scale(richness),
                                      evar = scale(evar),
                                      dominant_relative_abund_zero = scale(dominant_relative_abund_zero),
                               ))

coef(summary(resis_spei9_m4_std)) %>%
  data.frame() %>%
  rownames_to_column("Variable") %>%
  rename("SE" = "Std..Error",
         "p" = "Pr...t..") %>%
  # order things neatly
  mutate(Variable = factor(Variable,
                           levels = c("(Intercept)",
                                      "richness",
                                      "evar",
                                      "dominant_relative_abund_zero",
                                      "spei9_categoryExtreme wet",
                                      "richness:spei9_categoryExtreme wet")),
         # cosmetics in the names
         Variable_labels = factor(case_when(Variable %in% "(Intercept)" ~ "intercept",
                                            Variable %in% "evar" ~ "evenness",
                                            Variable %in% "dominant_relative_abund_zero" ~ "dominant abundance",
                                            Variable %in% "spei9_categoryExtreme wet" ~ "SPEI9: extreme wet",
                                            Variable %in% "richness:spei9_categoryExtreme wet" ~ "richness x SPEI9: extreme wet",
                                            .default = Variable)),
         Variable_labels = fct_reorder(Variable_labels, as.numeric(Variable)),
         # get significances and their plotting location (perhaps necessary to play with the +/- offset)
         stars = case_when(p < 0.001 ~ "***",
                           p > 0.001 & p <0.01 ~ "**",
                           p > 0.01  & p < 0.05 ~ "*",
                           p > 0.05 & p <0.1 ~ ".",
                           .default = ""),
         stars_location = case_when(Estimate < 0 ~ Estimate - SE - 0.1,
                                    Estimate > 0 ~ Estimate + SE + 0.1)
  ) %>%filter(Variable_labels!="intercept")%>%
  
  ggplot(aes(y = Variable_labels, x = Estimate)) +
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  labs(title = "log(resistance)",
       x = "estimate \u00B1 se") +
  scale_y_discrete(limits = rev) +
  geom_vline(xintercept = 0) +
  geom_col()+
  geom_errorbarh(aes(xmin = Estimate-SE, xmax = Estimate+SE), height = 0.2) +
  geom_text(aes(x = stars_location, label = stars))


#Resilience
#interaction based on hypothesis
resil_spei9<-lmer(log(resilience)~richness*dominant_relative_abund_zero+evar+
                    richness:spei9_category+dominant_relative_abund_zero:spei9_category+
                    evar:spei9_category+spei9_category+(1|site/experiment/uniqueid)
                  +(1|year), data=plot_ece_9_control, REML=F)
summary(resil_spei9)
anova(resil_spei9)
#without interaction of main predictors
resil_spei9_m1<-lmer(log(resilience)~richness+evar+dominant_relative_abund_zero+
                       richness:spei9_category+dominant_relative_abund_zero:spei9_category+
                       evar:spei9_category+spei9_category+(1|site/experiment/uniqueid)
                     +(1|year), data=plot_ece_9_control, REML=F)
summary(resil_spei9_m1)
anova(resil_spei9_m1)

#remove non-significant interaction
resil_spei9_m2<-lmer(log(resilience)~richness+evar+dominant_relative_abund_zero+
                       richness:spei9_category+
                       evar:spei9_category+spei9_category+(1|site/experiment/uniqueid)
                     +(1|year), data=plot_ece_9_control, REML=F)
anova(resil_spei9_m2)
summary(resil_spei9_m2)

anova(resil_spei9_m2, resil_spei9_m1,resil_spei9)#model selection
#refit additive model with REML
resil_spei9_m3<-lmer(log(resilience)~richness+evar+dominant_relative_abund_zero+
                       richness:spei9_category+
                       evar:spei9_category+spei9_category+(1|site/experiment/uniqueid)
                     +(1|year), data=plot_ece_9_control)
summary(resil_spei9_m3)
anova(resil_spei9_m3)
simres <- simulateResiduals(resil_spei9_m3)
plot(simres)
check_model(resil_spei9_m3)

#plot best model
ggpredict(model = resil_spei9_m3, terms = c("richness", "spei9_category"), back_transform = F ) %>%
  plot(show_data = TRUE)+
  labs(x="richness")
ggpredict(model = resil_spei9_m3, terms = c("evar", "spei9_category"), back_transform = F) %>%
  plot(show_data = TRUE)+
  labs(x="evenness")
ggpredict(model = resil_spei9_m3, terms = "dominant_relative_abund_zero", back_transform = F) %>%
  plot(show_data = TRUE)+
  labs(x="Relative abundance of dominant species")

# Plot standardized effect sizes
plot_model(resil_spei9_m3,
           type = "std",
           rm.terms = "measurement_scale_cover [0.3,0.4,1,10]",
           vline.color = "black",
           sort.est = TRUE,
           ci.lvl = 0.95
) + ylim(-0.05, 0.05) + theme_bw()


# Serainas version
# I think sjPlot does average across factor levels of other variables. I am not sure how to correctly do that manually
# note: this is only relevant for richness, as there is the reichness times SPEI interaction
resil_spei9_m3_std <- update(resil_spei9_m3, 
                             data = plot_ece_9_control %>% 
                               mutate(richness = scale(richness),
                                      evar = scale(evar),
                                      dominant_relative_abund_zero = scale(dominant_relative_abund_zero),
                               ))

coef(summary(resil_spei9_m3_std)) %>%
  data.frame() %>%
  rownames_to_column("Variable") %>%
  rename("SE" = "Std..Error",
         "p" = "Pr...t..") %>%
  # order things neatly
  mutate(Variable = factor(Variable,
                           levels = c("(Intercept)",
                                      "richness",
                                      "evar",
                                      "dominant_relative_abund_zero",
                                      "spei9_categoryExtreme wet",
                                      "richness:spei9_categoryExtreme wet",
                                      "evar:spei9_categoryExtreme wet")),
         # cosmetics in the names
         Variable_labels = factor(case_when(Variable %in% "(Intercept)" ~ "intercept",
                                            Variable %in% "evar" ~ "evenness",
                                            Variable %in% "dominant_relative_abund_zero" ~ "dominant abundance",
                                            Variable %in% "spei9_categoryExtreme wet" ~ "SPEI9: extreme wet",
                                            Variable %in% "richness:spei9_categoryExtreme wet" ~ "richness x SPEI9: extreme wet",
                                            Variable %in% "evar:spei9_categoryExtreme wet" ~ "evenness x SPEI9: extreme wet",
                                            .default = Variable)),
         Variable_labels = fct_reorder(Variable_labels, as.numeric(Variable)),
         # get significances and their plotting location (perhaps necessary to play with the +/- offset)
         stars = case_when(p < 0.001 ~ "***",
                           p > 0.001 & p <0.01 ~ "**",
                           p > 0.01  & p < 0.05 ~ "*",
                           p > 0.05 & p <0.1 ~ ".",
                           .default = ""),
         stars_location = case_when(Estimate < 0 ~ Estimate - SE - 0.1,
                                    Estimate > 0 ~ Estimate + SE + 0.1)
  ) %>%filter(Variable_labels!="intercept")%>%
  
  ggplot(aes(y = Variable_labels, x = Estimate)) +
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  labs(title = "log(resilience)",
       x = "estimate \u00B1 se") +
  scale_y_discrete(limits = rev) +
  geom_vline(xintercept = 0) +
  geom_col()+
  geom_errorbarh(aes(xmin = Estimate-SE, xmax = Estimate+SE), height = 0.2) +
  geom_text(aes(x = stars_location, label = stars))

#analysis using SPEI6####
# Analysis 1: resistance ----
## Control plot only ----
#using predictors without scaling####

#without interaction of main predictors
resis_unscaled_model<-lmer(log(resistance)~richness+evar+dominant_relative_abund_zero+
                              measurement_scale_cover+richness:spei6_category+dominant_relative_abund_zero:spei6_category+
                             evar:spei6_category+spei6_category+(1|site/experiment/uniqueid)
                            +(1|year), data=plot_ece_control, REML=F)
summary(resis_unscaled_model)
anova(resis_unscaled_model)
simres <- simulateResiduals(resis_unscaled_model)
plot(simres)
check_model(resis_unscaled_model)


#interaction based on hypothesis
resis_unscaled_model1<-lmer(log(resistance)~richness*dominant_relative_abund_zero+evar+
                             measurement_scale_cover+richness:spei6_category+dominant_relative_abund_zero:spei6_category+
                           evar:spei6_category+spei6_category+(1|site/experiment/uniqueid)
                           +(1|year), data=plot_ece_control, REML=F)
summary(resis_unscaled_model1)
anova(resis_unscaled_model1)
simres <- simulateResiduals(resis_unscaled_model1)
plot(simres)
check_model(resis_unscaled_model)
anova(resis_unscaled_model1,resis_unscaled_model) # model without interaction of main predictor better

#model selection without interaction of main predictors
#remove non-significant interaction
resis_unscaled_model2<-lmer(log(resistance)~richness+evar+dominant_relative_abund_zero+
                             measurement_scale_cover+richness:spei6_category+
                             evar:spei6_category+spei6_category+(1|site/experiment/uniqueid)
                           +(1|year), data=plot_ece_control, REML=F)
anova(resis_unscaled_model2)
#compare models
anova(resis_unscaled_model2,resis_unscaled_model)#new model better

#remove other non-significant interactions
resis_unscaled_model3<-lmer(log(resistance)~richness+evar+dominant_relative_abund_zero+
                              measurement_scale_cover+richness:spei6_category+spei6_category+(1|site/experiment/uniqueid)
                            +(1|year), data=plot_ece_control, REML=F)
anova(resis_unscaled_model3)
anova(resis_unscaled_model3,resis_unscaled_model2)#dropped interaction does not sig impact model

#refit additive model with REML
resis_unscaled_model4<-lmer(log(resistance)~richness+evar+dominant_relative_abund_zero+
                              measurement_scale_cover+richness:spei6_category+spei6_category+(1|site/experiment/uniqueid)
                            +(1|year), data=plot_ece_control)
summary(resis_unscaled_model4)
anova(resis_unscaled_model4)
simres <- simulateResiduals(resis_unscaled_model4)
plot(simres)
check_model(resis_unscaled_model4)
#log transforming evenness does not help deal with the low sample at high evenness
#might consider running a model without the high evenness values.


#plot best model
ggpredict(model = resis_unscaled_model4, terms = c("richness", "spei6_category"), back_transform = F ) %>%
  plot(show_data = TRUE)+
  labs(x="richness")
ggpredict(model = resis_unscaled_model4, terms = "evar", back_transform = F) %>%
  plot(show_data = TRUE)+
  labs(x="evenness")
ggpredict(model = resis_unscaled_model4, terms = "richness", back_transform = F) %>%
  plot(show_data = TRUE)+
  labs(x="richness")
ggpredict(model = resis_unscaled_model4, terms = "dominant_relative_abund_zero", back_transform = F) %>%
  plot(show_data = TRUE)+
  labs(x="Relative abundance of dominant species")

ggplot(plot_ece_control, aes(richness, log(resistance)))+
  geom_point()+
  geom_smooth(method = "lm")
ggplot(plot_ece_control, aes(dominant_relative_abund_zero, log(resistance)))+
  geom_point()+
  geom_smooth(method = "lm")
ggplot(plot_ece_control, aes(evar, log(resistance)))+
  geom_point()+
  geom_smooth(method = "lm")

# Plot standardized effect sizes
plot_model(resis_unscaled_model4,
           type = "std",
           rm.terms = "measurement_scale_cover [0.3,0.4,1,10]",
           vline.color = "black",
           sort.est = TRUE,
           ci.lvl = 0.95
) + ylim(-0.05, 0.05) + theme_bw()


#resilience model based on hypotheis
#interaction based on hypothesis
resil_unscaled_model<-lmer(log(resilience)~richness*dominant_relative_abund_zero+evar+
                              measurement_scale_cover+richness:spei6_category+dominant_relative_abund_zero:spei6_category+
                              evar:spei6_category+spei6_category+(1|site/experiment/uniqueid)
                            +(1|year), data=plot_ece_control, REML=F)
summary(resil_unscaled_model)
anova(resil_unscaled_model)
simres <- simulateResiduals(resil_unscaled_model)
plot(simres)
check_model(resil_unscaled_model)


#removing non-significant interactions
resil_unscaled_model1<-lmer(log(resilience)~richness+dominant_relative_abund_zero+evar+
                             measurement_scale_cover+richness:spei6_category+dominant_relative_abund_zero:spei6_category+
                             evar:spei6_category+spei6_category+(1|site/experiment/uniqueid)
                           +(1|year), data=plot_ece_control, REML=F)
summary(resil_unscaled_model1)
anova(resil_unscaled_model,resil_unscaled_model1) #checking if removing interaction was important

#removing non-significant interactions
resil_unscaled_model2<-lmer(log(resilience)~richness+dominant_relative_abund_zero+evar+
                              measurement_scale_cover+richness:spei6_category+
                              evar:spei6_category+spei6_category+(1|site/experiment/uniqueid)
                            +(1|year), data=plot_ece_control, REML=F)
summary(resil_unscaled_model2)
anova(resil_unscaled_model2,resil_unscaled_model1)

#Refit model with REML
resil_unscaled_model3<-lmer(log(resilience)~richness+dominant_relative_abund_zero+evar+
                              measurement_scale_cover+richness:spei6_category+
                              evar:spei6_category+spei6_category+(1|site/experiment/uniqueid)
                            +(1|year), data=plot_ece_control)


summary(resil_unscaled_model3)
anova(resil_unscaled_model3)
mres <- simulateResiduals(resil_unscaled_model3)
plot(simres)
check_model(resil_unscaled_model3)

#plot best model
ggpredict(model = resil_unscaled_model3, terms = c("richness", "spei6_category"), back_transform = F ) %>%
  plot(show_data = TRUE)+
  labs(x="richness")
ggpredict(model = resil_unscaled_model3, terms = c("evar", "spei6_category"), back_transform = F) %>%
  plot(show_data = TRUE)+
  labs(x="evenness")
ggpredict(model = resil_unscaled_model3, terms = "richness", back_transform = F) %>%
  plot(show_data = TRUE)+
  labs(x="richness")
ggpredict(model = resil_unscaled_model3, terms = "evar", back_transform = F) %>%
  plot(show_data = TRUE)+
  labs(x="evenness")

# Plot standardized effect sizes
plot_model(resil_unscaled_model3,
           type = "std",
           rm.terms = "measurement_scale_cover [0.3,0.4,1,10]",
           vline.color = "black",
           sort.est = TRUE,
           ci.lvl = 0.95
) + ylim(-0.05, 0.05) + theme_bw()


#investigating predictors at each site
#select Cedar Creek
plot_ece_control_cdr<-plot_ece_control%>%
  filter(site=="CDR")
#checking measurement scale categories
unique(plot_ece_control_cdr$measurement_scale_cover)
#resisitance
resis_cdr<-lmer(log(resistance)~richness*dominant_relative_abund_zero+evar+
                  measurement_scale_cover+richness:spei6_category+dominant_relative_abund_zero:spei6_category+
                  evar:spei6_category+spei6_category+(1|experiment/uniqueid)
                +(1|year), data=plot_ece_control_cdr, REML=F)
summary(resis_cdr)

#remove non-significant interactions
resis_cdr1<-lmer(log(resistance)~richness+dominant_relative_abund_zero+evar+
                  measurement_scale_cover+richness:spei6_category+dominant_relative_abund_zero:spei6_category+
                  evar:spei6_category+spei6_category+(1|experiment/uniqueid)
                +(1|year), data=plot_ece_control_cdr, REML=F)
summary(resis_cdr1)
anova(resis_cdr,resis_cdr1)

resis_cdr2<-lmer(log(resistance)~richness+dominant_relative_abund_zero+evar+
                   measurement_scale_cover+richness:spei6_category+
                   evar:spei6_category+spei6_category+(1|experiment/uniqueid)
                 +(1|year), data=plot_ece_control_cdr, REML=F)
summary(resis_cdr2)


resis_cdr3<-lmer(log(resistance)~richness+dominant_relative_abund_zero+evar+
                   measurement_scale_cover+richness:spei6_category+
                  +spei6_category+(1|experiment/uniqueid)
                 +(1|year), data=plot_ece_control_cdr, REML=F)
summary(resis_cdr3)
anova(resis_cdr3,resis_cdr2,resis_cdr1,resis_cdr)#model selection

#Refit best model with REML
resis_cdr4<-lmer(log(resistance)~richness+dominant_relative_abund_zero+evar+
                   measurement_scale_cover+richness:spei6_category+
                   +spei6_category+(1|experiment/uniqueid)
                 +(1|year), data=plot_ece_control_cdr)
summary(resis_cdr1)
anova(resis_cdr4)
simres <- simulateResiduals(resis_cdr4)
plot(simres)#might not look good enough
check_model(resis_cdr4)

#plot best model
ggpredict(model = resis_cdr4, terms = c("richness", "spei6_category"), back_transform = F ) %>%
  plot(show_data = TRUE)+
  labs(x="richness")
ggpredict(model = resis_cdr4, terms = "dominant_relative_abund_zero", back_transform = F) %>%
  plot(show_data = TRUE)+
  labs(x="relative abundance of dominant species")
ggpredict(model = resis_cdr4, terms = "spei6_category", back_transform = F) %>%
  plot(show_data = TRUE)+
  labs(x="spei6_category")

#Cedar creek resilience
resil_cdr<-lmer(log(resilience)~richness*dominant_relative_abund_zero+evar+
                  measurement_scale_cover+richness:spei6_category+dominant_relative_abund_zero:spei6_category+
                  evar:spei6_category+spei6_category+(1|experiment/uniqueid)
                +(1|year), data=plot_ece_control_cdr, REML=F)
summary(resil_cdr)

#remove non-significant interactions
resil_cdr1<-lmer(log(resilience)~richness+dominant_relative_abund_zero+evar+
                   measurement_scale_cover+richness:spei6_category+dominant_relative_abund_zero:spei6_category+
                   evar:spei6_category+spei6_category+(1|experiment/uniqueid)
                 +(1|year), data=plot_ece_control_cdr, REML=F)
summary(resil_cdr1)

resil_cdr2<-lmer(log(resilience)~richness+dominant_relative_abund_zero+evar+
                   measurement_scale_cover+richness:spei6_category+
                   evar:spei6_category+spei6_category+(1|experiment/uniqueid)
                 +(1|year), data=plot_ece_control_cdr, REML=F)
summary(resil_cdr2)


resil_cdr3<-lmer(log(resilience)~richness+dominant_relative_abund_zero+evar+
                   measurement_scale_cover+richness:spei6_category+
                   +spei6_category+(1|experiment/uniqueid)
                 +(1|year), data=plot_ece_control_cdr, REML=F)
summary(resil_cdr3)
anova(resil_cdr3,resil_cdr2,resil_cdr1,resil_cdr)#model selection

#Refit best model with REML
resil_cdr4<-lmer(log(resilience)~richness+dominant_relative_abund_zero+evar+
                   measurement_scale_cover+richness:spei6_category+
                   +spei6_category+(1|experiment/uniqueid)
                 +(1|year), data=plot_ece_control_cdr)
summary(resil_cdr1)
anova(resil_cdr4)
simres <- simulateResiduals(resil_cdr4)
plot(simres)#look good enough
check_model(resil_cdr4)

ggpredict(model = resil_cdr4, terms = c("richness", "spei6_category"), back_transform = F ) %>%
  plot(show_data = TRUE)+
  labs(x="richness")
ggpredict(model = resil_cdr4, terms = "dominant_relative_abund_zero", back_transform = F) %>%
  plot(show_data = TRUE)+
  labs(x="relative abundance of dominant species")


#select KBS
plot_ece_control_kbs<-plot_ece_control%>%
  filter(site=="KBS")
#checking categories


#resisitance
resis_kbs<-lmer(log(resistance)~richness*dominant_relative_abund_zero+evar+measurement_scale_cover+
                  richness:spei6_category+dominant_relative_abund_zero:spei6_category+
                  evar:spei6_category+spei6_category+(1|experiment/uniqueid)
                +(1|year), data=plot_ece_control_kbs, REML=F)#it seems there is not enough replicate for each measurement scale cover
summary(resis_kbs)

#remove non-significant interactions
resis_kbs1<-lmer(log(resistance)~richness+dominant_relative_abund_zero+evar+
                   richness:spei6_category+dominant_relative_abund_zero:spei6_category+
                   evar:spei6_category+spei6_category+(1|experiment/uniqueid)
                 +(1|year), data=plot_ece_control_kbs, REML=F)
summary(resis_kbs1)


resis_kbs2<-lmer(log(resistance)~richness+dominant_relative_abund_zero+evar+
                   richness:spei6_category+
                   evar:spei6_category+spei6_category+(1|experiment/uniqueid)
                 +(1|year), data=plot_ece_control_kbs, REML=F)
summary(resis_kbs2)


resis_kbs3<-lmer(log(resistance)~richness+dominant_relative_abund_zero+evar+
                   richness:spei6_category+
                   +spei6_category+(1|experiment/uniqueid)
                 +(1|year), data=plot_ece_control_kbs, REML=F)
summary(resis_kbs3)
anova(resis_kbs3)#model selection
check_model(resis_kbs3)


resis_kbs4<-lmer(log(resistance)~richness+dominant_relative_abund_zero+evar+
                   +spei6_category+(1|experiment/uniqueid)
                 +(1|year), data=plot_ece_control_kbs, REML=F)
summary(resis_kbs4)

resis_kbs5<-lmer(log(resistance)~richness+dominant_relative_abund_zero+
                   +spei6_category+(1|experiment/uniqueid)
                 +(1|year), data=plot_ece_control_kbs, REML=F)
summary(resis_kbs5)

resis_kbs6<-lmer(log(resistance)~richness+dominant_relative_abund_zero+
                   +(1|experiment/uniqueid)
                 +(1|year), data=plot_ece_control_kbs, REML=F)
summary(resis_kbs6)#none of the predictors were significant
check_model(resis_kbs6)

#resilience
resil_kbs<-lmer(log(resilience)~richness*dominant_relative_abund_zero+evar+measurement_scale_cover+
                  richness:spei6_category+dominant_relative_abund_zero:spei6_category+
                  evar:spei6_category+spei6_category+(1|experiment/uniqueid)
                +(1|year), data=plot_ece_control_kbs, REML=F)#it seems there is not enough replicate for each measurement scale cover
summary(resil_kbs)

#remove non-significant interactions
resil_kbs1<-lmer(log(resilience)~richness+dominant_relative_abund_zero+evar+
                   richness:spei6_category+dominant_relative_abund_zero:spei6_category+
                   evar:spei6_category+spei6_category+(1|experiment/uniqueid)
                 +(1|year), data=plot_ece_control_kbs, REML=F)
summary(resil_kbs1)


resil_kbs2<-lmer(log(resilience)~richness+dominant_relative_abund_zero+evar+
                   richness:spei6_category+
                   evar:spei6_category+spei6_category+(1|experiment/uniqueid)
                 +(1|year), data=plot_ece_control_kbs, REML=F)
summary(resil_kbs2)


resil_kbs3<-lmer(log(resilience)~richness+dominant_relative_abund_zero+evar+
                   evar:spei6_category
                   +spei6_category+(1|experiment/uniqueid)
                 +(1|year), data=plot_ece_control_kbs, REML=F)
summary(resil_kbs3)
anova(resil_kbs3)#model selection
check_model(resil_kbs3)

#refit with REML
resil_kbs4<-lmer(log(resilience)~richness+dominant_relative_abund_zero+evar+
                   evar:spei6_category
                 +spei6_category+(1|experiment/uniqueid)
                 +(1|year), data=plot_ece_control_kbs)
summary(resil_kbs4)
anova(resil_kbs4)
check_model(resil_kbs4)

ggpredict(model = resil_kbs4, terms = c("evar", "spei6_category"), back_transform = F ) %>%
  plot(show_data = TRUE)+
  labs(x="richness")

#Konza 
plot_ece_control_knz<-plot_ece_control%>%
  filter(site=="KNZ")
#checking measurement scale categories
unique(plot_ece_control_knz$measurement_scale_cover)
#resisitance
resis_knz<-lmer(log(resistance)~richness*dominant_relative_abund_zero+evar+
                  measurement_scale_cover+richness:spei6_category+dominant_relative_abund_zero:spei6_category+
                  evar:spei6_category+spei6_category+(1|experiment/uniqueid)
                +(1|year), data=plot_ece_control_knz, REML=F)
summary(resis_knz)

#remove non-significant interactions
resis_knz1<-lmer(log(resistance)~richness+dominant_relative_abund_zero+evar+
                   measurement_scale_cover+richness:spei6_category+dominant_relative_abund_zero:spei6_category+
                   evar:spei6_category+spei6_category+(1|experiment/uniqueid)
                 +(1|year), data=plot_ece_control_knz, REML=F)
summary(resis_knz1)


resis_knz2<-lmer(log(resistance)~richness+dominant_relative_abund_zero+evar+
                   measurement_scale_cover+richness:spei6_category+
                   evar:spei6_category+spei6_category+(1|experiment/uniqueid)
                 +(1|year), data=plot_ece_control_knz, REML=F)
summary(resis_knz2)


resis_knz3<-lmer(log(resistance)~richness+dominant_relative_abund_zero+evar+
                   measurement_scale_cover+richness:spei6_category+
                   +spei6_category+(1|experiment/uniqueid)
                 +(1|year), data=plot_ece_control_knz, REML=F)
summary(resis_knz3)

resis_knz4<-lmer(log(resistance)~richness+dominant_relative_abund_zero+evar+
                   measurement_scale_cover+
                   +spei6_category+(1|experiment/uniqueid)
                 +(1|year), data=plot_ece_control_knz, REML=F)
summary(resis_knz4)

resis_knz5<-lmer(log(resistance)~richness+dominant_relative_abund_zero+
                   measurement_scale_cover+
                   +(1|experiment/uniqueid)
                 +(1|year), data=plot_ece_control_knz, REML=F)
summary(resis_knz5)#no predictor significant
check_model(resis_knz5)
anova(resis_knz5,resis_knz4,resis_knz3,resis_knz2,resis_knz1,resis_knz)


#resilience
resil_knz<-lmer(log(resilience)~richness*dominant_relative_abund_zero+evar+
                  measurement_scale_cover+richness:spei6_category+dominant_relative_abund_zero:spei6_category+
                  evar:spei6_category+spei6_category+(1|experiment/uniqueid)
                +(1|year), data=plot_ece_control_knz, REML=F)
summary(resil_knz)

#remove non-significant interactions
resil_knz1<-lmer(log(resilience)~richness+dominant_relative_abund_zero+evar+
                   measurement_scale_cover+richness:spei6_category+dominant_relative_abund_zero:spei6_category+
                   evar:spei6_category+spei6_category+(1|experiment/uniqueid)
                 +(1|year), data=plot_ece_control_knz, REML=F)
summary(resil_knz1)


resil_knz2<-lmer(log(resilience)~richness+dominant_relative_abund_zero+evar+
                   measurement_scale_cover+richness:spei6_category+
                   evar:spei6_category+spei6_category+(1|experiment/uniqueid)
                 +(1|year), data=plot_ece_control_knz, REML=F)
summary(resil_knz2)


resil_knz3<-lmer(log(resilience)~richness+dominant_relative_abund_zero+evar+
                   measurement_scale_cover+richness:spei6_category+
                   +spei6_category+(1|experiment/uniqueid)
                 +(1|year), data=plot_ece_control_knz, REML=F)
summary(resil_knz3)

resil_knz4<-lmer(log(resilience)~richness+dominant_relative_abund_zero+evar+
                   measurement_scale_cover+
                   +spei6_category+(1|experiment/uniqueid)
                 +(1|year), data=plot_ece_control_knz, REML=F)
summary(resil_knz4)
anova(resil_knz4,resil_knz3,resil_knz2,resil_knz1,resil_knz)#model selection
#refit with REML
resil_knz5<-lmer(log(resilience)~richness+dominant_relative_abund_zero+evar+
                   measurement_scale_cover+
                   +spei6_category+(1|experiment/uniqueid)
                 +(1|year), data=plot_ece_control_knz)
summary(resil_knz5)
anova(resil_knz5)
check_model(resil_knz5)

ggpredict(model = resil_knz5, terms = "richness", back_transform = F ) %>%
  plot(show_data = TRUE)+
  labs(x="richness")
ggpredict(model = resil_knz5, terms = "evar", back_transform = F ) %>%
  plot(show_data = TRUE)+
  labs(x="evenness")




#goal: determine an ideal random effect structure
#control only model with all possible main effects and interactions
#using berger parker as dominance metric####
#using scaled predictors
resis_random_model<-lmer(log(resistance) ~ richness_scaled+berger_parker_scaled+evar_scaled+spei6_category+richness_scaled:spei6_category+
                           richness_scaled:berger_parker_scaled+berger_parker_scaled:spei6_category+richness_scaled:berger_parker_scaled:spei6_category+
                           evar_scaled:spei6_category+(1|year) + (1|site/experiment/uniqueid), data = plot_ece_control)
summary(resis_random_model) 
check_model(resis_random_model)#diagnostic visuals
anova(resis_random_model)
#remove site
resis_random_model1<-lmer(log(resistance) ~ richness_scaled+berger_parker_scaled+evar_scaled+spei6_category+richness_scaled:spei6_category+
                            richness_scaled:berger_parker_scaled+berger_parker_scaled:spei6_category+richness_scaled:berger_parker_scaled:spei6_category+
                            evar_scaled:spei6_category+(1|year) + (1|experiment/uniqueid), data = plot_ece_control)
summary(resis_random_model1)
check_model(resis_random_model1)
anova(resis_random_model1)

AIC(resis_random_model1,resis_random_model)#model without site and  has the lowest AIC

#selecting fixed effects but retaining hypothesis
#refit model with "ML"
resis_fixed_model<-lmer(log(resistance) ~ richness_scaled+berger_parker_scaled+evar_scaled+spei6_category+richness_scaled:spei6_category+
                          richness_scaled:berger_parker_scaled+berger_parker_scaled:spei6_category+richness_scaled:berger_parker_scaled:spei6_category+
                          evar_scaled:spei6_category+(1|year) + (1|experiment/uniqueid), data = plot_ece_control, REML = F)
anova(resis_fixed_model)
fixed_model1<-update(resis_fixed_model, .~.-richness_scaled:berger_parker_scaled:spei6_category)
#determine the significance of the dropped term
anova(resis_fixed_model,fixed_model1)#okay to drop term
#drop more terms
fixed_model2<-update(fixed_model1, .~.-berger_parker_scaled:spei6_category)
anova(fixed_model1, fixed_model2)#okay to drop term
anova(fixed_model2)

#more terms to drop
fixed_model3<-update(fixed_model2, .~.-richness_scaled:spei6_category)
anova(fixed_model2, fixed_model3)
anova(fixed_model3)
check_model(fixed_model3)

#model suggests dominance and richness interaction is not significant
#dropping this interaction
fixed_model4<-update(fixed_model3, .~.-berger_parker_scaled:richness_scaled)
anova(fixed_model3,fixed_model4)#okay to drop
anova(fixed_model4)
check_model(fixed_model4)

fixed_model5<-update(fixed_model4, .~.-evar_scaled:spei6_category)
anova(fixed_model4, fixed_model5)
anova(fixed_model5)
check_model(fixed_model5)

#left richness and dominance in the model since it relates to our hypothesis
#Refit model with REML
fixed_model6<-lmer(log(resistance)~richness_scaled+berger_parker_scaled+evar_scaled+spei6_category+
                         (1|year)+(1|experiment/uniqueid), data=plot_ece_control)
anova(fixed_model6)
simres <- simulateResiduals(fixed_model6)
plot(simres)
check_model(fixed_model6)#some spread in variance and normality issue #need covariance structure?
#can't set variance with lme4 so pivot to nlme
library(nlme)
fixed_model7<-lme(log(resistance)~richness_scaled+berger_parker_scaled+evar_scaled+spei6_category,
                    random = ~1|experiment/uniqueid, method="REML", data=plot_ece_control, na.action = na.omit)
#does anyone know a way around setting multiple random effect with lme? can't set year as random effect with another random effect already stated
anova(fixed_model7)

#setting a fixed variance_covariance structure
vfix<-varFixed(~evar_scaled)
fixed_model8<-lme(log(resistance)~richness_scaled+berger_parker_scaled+evar_scaled+spei6_category,
                  random = ~1|experiment/uniqueid, method="REML", data=plot_ece_control, na.action = na.omit,weights= vfix, correlation = corCompSymm(form=~1|experiment/uniqueid))
anova(fixed_model8)
summary(fixed_model8)
check_model(fixed_model8)#not better
anova(fixed_model7,fixed_model8)

AIC(fixed_model7, fixed_model8) #model with the original structure has lower AIC


#runing the model with abundance of dominant species instead of berger parker####
resis_rand_dom_model<-lmer(log(resistance) ~ richness_scaled+dominant_relative_abund_zero_scaled+evar_scaled+spei6_category+richness_scaled:spei6_category+
                           richness_scaled:dominant_relative_abund_zero_scaled+dominant_relative_abund_zero_scaled:spei6_category+richness_scaled:dominant_relative_abund_zero_scaled:spei6_category+
                           evar_scaled:spei6_category+(1|year) + (1|site/experiment/uniqueid), data = plot_ece_control)
summary(resis_rand_dom_model) 
check_model(resis_rand_dom_model)#diagnostic visuals
anova(resis_rand_dom_model)
#remove site
resis_rand_dom_model1<-lmer(log(resistance) ~ richness_scaled+dominant_relative_abund_zero_scaled+evar_scaled+spei6_category+richness_scaled:spei6_category+
                            richness_scaled:dominant_relative_abund_zero_scaled+dominant_relative_abund_zero_scaled:spei6_category+richness_scaled:dominant_relative_abund_zero_scaled:spei6_category+
                            evar_scaled:spei6_category+(1|year) + (1|experiment/uniqueid), data = plot_ece_control)
summary(resis_rand_dom_model1) 
check_model(resis_rand_dom_model1)
anova(resis_rand_dom_model1)
#remove year
resis_rand_dom_model2<-lmer(log(resistance) ~ richness_scaled+dominant_relative_abund_zero_scaled+evar_scaled+spei6_category+richness_scaled:spei6_category+
                              richness_scaled:dominant_relative_abund_zero_scaled+dominant_relative_abund_zero_scaled:spei6_category+richness_scaled:dominant_relative_abund_zero_scaled:spei6_category+
                              evar_scaled:spei6_category+ (1|experiment/uniqueid), data = plot_ece_control)
summary(resis_rand_dom_model2)  
check_model(resis_rand_dom_model2)
anova(resis_rand_dom_model2)

AIC(resis_rand_dom_model,resis_rand_dom_model1,resis_rand_dom_model2)#model with year but without site and  has the lowest AIC

#selecting fixed effects but retaining hypothesis
#refit model with "ML", then proceed with stepwise selection
resis_fixed_dom_model<-lmer(log(resistance) ~ richness_scaled+dominant_relative_abund_zero_scaled+evar_scaled+spei6_category+richness_scaled:spei6_category+
                          richness_scaled:dominant_relative_abund_zero_scaled+dominant_relative_abund_zero_scaled:spei6_category+richness_scaled:dominant_relative_abund_zero_scaled:spei6_category+
                          evar_scaled:spei6_category+(1|year) + (1|experiment/uniqueid), data = plot_ece_control, REML = F)
anova(resis_fixed_dom_model)
fixed_dom_model1<-update(resis_fixed_dom_model, .~.-richness_scaled:dominant_relative_abund_zero_scaled:spei6_category)
anova(fixed_dom_model1, resis_fixed_dom_model)#okay to proceed
anova(fixed_dom_model1)

fixed_dom_model2<-update(fixed_dom_model1, .~.-evar_scaled:spei6_category)
anova(fixed_dom_model1, fixed_dom_model2)
anova(fixed_dom_model2)



fixed_dom_model3<-update(fixed_dom_model2, .~.-richness_scaled:dominant_relative_abund_zero_scaled)
anova(fixed_dom_model3, fixed_dom_model2)
anova(fixed_dom_model3)

fixed_dom_model4<-update(fixed_dom_model3, .~.-dominant_relative_abund_zero_scaled:spei6_category)
anova(fixed_dom_model3, fixed_dom_model4)
anova(fixed_dom_model4)

fixed_dom_model5<-update(fixed_dom_model4, .~.-richness_scaled:spei6_category)
anova(fixed_dom_model5, fixed_dom_model4)
anova(fixed_dom_model5)

fixed_dom_model6<-update(fixed_dom_model5, .~.-dominant_relative_abund_zero_scaled)
anova(fixed_dom_model5, fixed_dom_model6)#fixed_dom_model5 has lower AIC and dominance is part of our initial hypothesis
anova(fixed_dom_model6)

#refit model with REML
fixed_dom_model5_update<-update(fixed_dom_model5, REML=T)
anova(fixed_dom_model5_update)#lower AIC than model using berger parker
summary(fixed_dom_model5_update)
check_model(fixed_dom_model5_update)
simres <- simulateResiduals(fixed_dom_model5_update)
plot(simres)#not that bad, but not the best either-thoughts?



#plot best model
ggpredict(model = fixed_dom_model5_update, terms = c("evar_scaled", "spei6_category"), back_transform = F ) %>%
  plot(show_data = TRUE)+
  labs(x="evenness")
ggpredict(model = fixed_dom_model5_update, terms = "evar_scaled", back_transform = F) %>%
  plot(show_data = TRUE)+
  labs(x="evenness")



#####
#runing resilience model using only control with abundance of dominant species instead of berger parker####
resil_rand_dom_model<-lmer(log(resilience) ~ richness_scaled+dominant_relative_abund_zero_scaled+evar_scaled+spei6_category+richness_scaled:spei6_category+
                             richness_scaled:dominant_relative_abund_zero_scaled+dominant_relative_abund_zero_scaled:spei6_category+richness_scaled:dominant_relative_abund_zero_scaled:spei6_category+
                             evar_scaled:spei6_category+(1|year) + (1|site/experiment/uniqueid), data = plot_ece_control)
summary(resil_rand_dom_model) 
check_model(resil_rand_dom_model)#diagnostic visuals
simres <- simulateResiduals(resil_rand_dom_model)
plot(simres)
anova(resil_rand_dom_model)
#remove site
resil_rand_dom_model1<-lmer(log(resilience) ~ richness_scaled+dominant_relative_abund_zero_scaled+evar_scaled+spei6_category+richness_scaled:spei6_category+
                              richness_scaled:dominant_relative_abund_zero_scaled+dominant_relative_abund_zero_scaled:spei6_category+richness_scaled:dominant_relative_abund_zero_scaled:spei6_category+
                              evar_scaled:spei6_category+(1|year) + (1|experiment), data = plot_ece_control)
summary(resil_rand_dom_model1) 
check_model(resil_rand_dom_model1)
anova(resil_rand_dom_model1)
#remove year
resil_rand_dom_model2<-lmer(log(resilience) ~ richness_scaled+dominant_relative_abund_zero_scaled+evar_scaled+spei6_category+richness_scaled:spei6_category+
                              richness_scaled:dominant_relative_abund_zero_scaled+dominant_relative_abund_zero_scaled:spei6_category+richness_scaled:dominant_relative_abund_zero_scaled:spei6_category+
                              evar_scaled:spei6_category+ (1|experiment), data = plot_ece_control)
summary(resil_rand_dom_model2)  
check_model(resil_rand_dom_model2)
anova(resil_rand_dom_model2)

AIC(resil_rand_dom_model2,resil_rand_dom_model1,resil_rand_dom_model)#model with year but without site and unique id  has the lowest AIC

#selecting fixed effects but retaining hypothesis
#refit model with "ML", then proceed with stepwise selection
resil_fixed_dom_model<-lmer(log(resilience) ~ richness_scaled+dominant_relative_abund_zero_scaled+evar_scaled+spei6_category+richness_scaled:spei6_category+
                                  richness_scaled:dominant_relative_abund_zero_scaled+dominant_relative_abund_zero_scaled:spei6_category+richness_scaled:dominant_relative_abund_zero_scaled:spei6_category+
                                  evar_scaled:spei6_category+(1|year) + (1|experiment), data = plot_ece_control)
summary(resil_fixed_dom_model) 
anova(resil_fixed_dom_model)
resil_fixed_dom_model1<-update(resil_fixed_dom_model, .~.-dominant_relative_abund_zero_scaled:spei6_category)
anova(resil_fixed_dom_model1, resil_fixed_dom_model)#okay to proceed
anova(resil_fixed_dom_model1)

resil_fixed_dom_model2<-update(resil_fixed_dom_model1, .~.-evar_scaled:spei6_category)
anova(resil_fixed_dom_model1, resil_fixed_dom_model2)
anova(resil_fixed_dom_model2)



resil_fixed_dom_model3<-update(resil_fixed_dom_model2, .~.-richness_scaled:spei6_category)
anova(resil_fixed_dom_model3, resil_fixed_dom_model2)
anova(resil_fixed_dom_model3)

resil_fixed_dom_model4<-update(resil_fixed_dom_model3, .~.-richness_scaled:dominant_relative_abund_zero_scaled)
anova(resil_fixed_dom_model3, resil_fixed_dom_model4)
anova(resil_fixed_dom_model4)

resil_fixed_dom_model5<-update(resil_fixed_dom_model4, .~.-richness_scaled:dominant_relative_abund_zero_scaled:spei6_category)
anova(resil_fixed_dom_model5, resil_fixed_dom_model4)
anova(resil_fixed_dom_model5)

resil_fixed_dom_model6<-update(resil_fixed_dom_model5, .~.-spei6_category)
anova(resil_fixed_dom_model6, resil_fixed_dom_model5)#resil_fixed_dom_model6 has lower AIC and dominance is part of our initial hypothesis
anova(resil_fixed_dom_model6)

#refit model with REML
resil_fixed_dom_model5_update<-update(resil_fixed_dom_model5, REML=T)
anova(resil_fixed_dom_model5_update)
summary(resil_fixed_dom_model5_update)
check_model(resil_fixed_dom_model5_update)
simres <- simulateResiduals(resil_fixed_dom_model5_update)
plot(simres)#not that bad, but not the best either-thoughts?



#plot best model
ggpredict(model = resil_fixed_dom_model5_update, terms = c("evar_scaled", "spei6_category"), back_transform = F ) %>%
  plot(show_data = TRUE)+
  labs(x="evenness")
ggpredict(model = resil_fixed_dom_model5_update, terms = "evar_scaled", back_transform = F) %>%
  plot(show_data = TRUE)+
  labs(x="evenness")
ggpredict(model = resil_fixed_dom_model5_update, terms = c("richness_scaled", "spei6_category"), back_transform = F ) %>%
  plot(show_data = TRUE)+
  labs(x="richness")

# Model with all variables (need to add experiment)
resist.control.full <- lmer(log10(resistance) ~ scale(richness)*scale(berger_parker)*site*spei6_category + year + (1|site/uniqueid), data = plot_ece_control, na.action = na.fail) # failed to converge

# Model with all variables (year and experiment as random effects)
resist.control.fuller <- lmer(log10(resistance) ~ scale(richness)*scale(berger_parker)*site*spei6_category + (1|year) + (1|site/experiment/uniqueid), data = plot_ece_control, na.action = na.fail) # failed to converge

# dredge(resist.control.full) # best model only has spei_category lol
# dredge(resist.control.fuller) # best model has richness*spei_category

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

# Drop dominance (but it's in our hypothesis so we might want to include it)
resist.control.rich <- lmer(log10(resistance) ~ scale(richness)*spei6_category + year + (1|site/uniqueid), data = plot_ece_control)
summary(resist.control.rich)
simres <- simulateResiduals(resist.control.rich)
plot(simres) # Very bad!

# Year as random effect
resist.control.rich.ny <- lmer(log10(resistance) ~ scale(richness)*spei6_category + (1|year) + (1|site/uniqueid), data = plot_ece_control)
summary(resist.control.rich.ny)
simres <- simulateResiduals(resist.control.rich.ny)
plot(simres) # Very bad!

# Drop interaction (but it's in our hypothesis so we might want to include it)
resist.control.rich.add <- lmer(log10(resistance) ~ scale(richness) + spei6_category + year + (1|site/uniqueid), data = plot_ece_control)
summary(resist.control.rich.add)
simres <- simulateResiduals(resist.control.rich.add)
plot(simres) # Very bad!

# Drop spei category
resist.control.rich.nc <- lmer(log10(resistance) ~ scale(richness) + year + (1|site/uniqueid), data = plot_ece_control)
summary(resist.control.rich.nc)
simres <- simulateResiduals(resist.control.rich.nc)
plot(simres) # Very bad!

# Only spei category
resist.control.spei <- lmer(log10(resistance) ~ spei6_category + year + (1|site/uniqueid), data = plot_ece_control) # singular fit, site explains zero variance
summary(resist.control.spei)
simres <- simulateResiduals(resist.control.spei)
plot(simres) # Very bad!

# Make a model based on hypotheses
resist.control.hyp <- lmer(log10(resistance) ~ scale(richness)*spei6_category + scale(berger_parker)*scale(richness) + scale(berger_parker)*spei6_category + year + (1|site/uniqueid), data = plot_ece_control)
summary(resist.control.hyp)
simres <- simulateResiduals(resist.control.hyp)
plot(simres) # Very bad!

# Year as random effect (hypothesis)
resist.control.hyp.y <- lmer(log10(resistance) ~ scale(richness)*spei6_category + scale(berger_parker)*scale(richness) + scale(berger_parker)*spei6_category + (1|year) + (1|site/experiment/uniqueid), data = plot_ece_control)
summary(resist.control.hyp.y) 
simres <- simulateResiduals(resist.control.hyp.y)
plot(simres) # Very bad!
vif(resist.control.hyp.y) # berger parker is bad
plot_ece_control %>%
  select(richness, berger_parker, resistance, resilience, spei6_category) %>%
  ggpairs()
D <- cooks.distance(resist.control.hyp.y)
which(D > 0.5) # no influential outliers

# Use richness and BP scaled at the experiment level
# Year as random effect (hypothesis)
resist.control.hyp.scale <- lmer(log10(resistance) ~ richness_scaled*spei6_category + berger_parker_scaled*richness_scaled + berger_parker_scaled*spei6_category + (1|year) + (1|site/experiment/uniqueid), data = plot_ece_control)
summary(resist.control.hyp.scale)
simres <- simulateResiduals(resist.control.hyp.scale)
plot(simres) # Very bad!

AICctab(resist.control.int, resist.control, resist.control.rich, resist.control.rich.add, resist.control.rich.nc, resist.control.spei, resist.control.hyp, resist.control.rich.ny)

# Plot best model
ggpredict(model = resist.control.rich.ny, terms = c("richness", "spei6_category"), back_transform = F) %>%
  plot(show_data = TRUE)

# Plot hypothesis model
ggpredict(model = resist.control.hyp.scale, terms = c("richness_scaled", "spei6_category"), back_transform = F) %>%
  plot(show_data = TRUE)
ggpredict(model = resist.control.hyp.scale, terms = c("berger_parker_scaled", "spei6_category"), back_transform = F) %>%
  plot(show_data = TRUE)
ggpredict(model = resist.control.hyp.scale, terms = c("richness_scaled", "berger_parker_scaled"), back_transform = F) %>%
  plot(show_data = TRUE)
ggpredict(model = resist.control.hyp.scale, terms = c("richness_scaled", "berger_parker_scaled", "spei6_category"), back_transform = F) %>%
  plot(show_data = TRUE, facet = T)



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


# Analysis 2 - resilience ----
## Control plots only ----
# Model with all variables (need to add experiment)
resil.control.full <- lmer(log10(resilience) ~ scale(richness)*scale(berger_parker)*site*spei6_category + year + (1|site/uniqueid), data = plot_ece_control, na.action = na.fail) # singular fit

# Model with all variables (year and experiment random effects)
resil.control.fuller <- lmer(log10(resilience) ~ scale(richness)*scale(berger_parker)*site*spei6_category + (1|year) + (1|site/experiment/uniqueid), data = plot_ece_control, na.action = na.fail) # failed to converge

# dredge(resil.control.full) # best model has richness*spei_cat, second best adds year
# dredge(resil.control.fuller) # best model has site*spei_cat

# Best model as predicted by dredge function (but including year)
resil.control.rich <- lmer(log10(resilience) ~ scale(richness)*spei6_category + year + (1|site/uniqueid), data = plot_ece_control) # singular fit
summary(resil.control.rich)
simres <- simulateResiduals(resil.control.rich)
plot(simres) # Very bad!

# Move year to random effects
resil.control.rich.ny <- lmer(log10(resilience) ~ scale(richness)*spei6_category + (1|year) + (1|site/uniqueid), data = plot_ece_control) # singular fit
summary(resil.control.rich.ny)
simres <- simulateResiduals(resil.control.rich.ny)
plot(simres) # Very bad!

# Additive only model
resil.control.add <- lmer(log10(resilience) ~ richness + berger_parker + site + spei6_category + year + (1|site/uniqueid), data = plot_ece_control) #singular fit

# Additive only model (year as random effect)
resil.control.add.ny <- lmer(log10(resilience) ~ richness + berger_parker + site + spei6_category + (1|year) + (1|site/uniqueid), data = plot_ece_control) #singular fit

# Three way interaction, drop site, scale richness and dominance
resil.control.int <- lmer(log10(resilience) ~ scale(richness)*scale(berger_parker)*spei6_category + year + (1|site/uniqueid), data = plot_ece_control)
summary(resil.control.int)

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

##### SEM #####
library(piecewiseSEM)

names(plot_ece_9_cn)

table ( plot_ece_9_cn$nutrients_grouped )

# get list of experiments that manipulated N
manipulated_n_exp = plot_ece_meta %>% filter (!is.na(nitrogen_amount))  %>% distinct(experiment) %>% select(experiment)

# get subset file for SEM:
plot_ece_meta_selSEM_init = plot_ece_9_cn %>% 
  select (site, experiment, uniqueid,spei9_category,nutrients_grouped ,
                          richness,dominant_relative_abund_zero,resistance,resilience) %>%  # only rows we care about
  filter (  experiment %in%  manipulated_n_exp[,1] & nutrients_grouped != "N" ) 

plot_ece_meta_selSEM = plot_ece_meta_selSEM_init[complete.cases(plot_ece_meta_selSEM_init ) ,]
plot_ece_meta_selSEM$spei9_category = as.factor(plot_ece_meta_selSEM$spei9_category)

table(plot_ece_meta_selSEM$nutrients_grouped)
plot_ece_meta_selSEM$nut_dummy = case_when(plot_ece_meta_selSEM$nutrients_grouped =="no_fertilizer" ~ 0,
                                           plot_ece_meta_selSEM$nutrients_grouped =="Nplus" ~ 1)
dim(plot_ece_meta_selSEM)
unique(plot_ece_meta_selSEM$spei9_category)
table(plot_ece_meta_selSEM$nitrogen_amount)
unique(plot_ece_9_cn$spei9_category)
unique(plot_ece_9_cn$nitrogen_amount)

table (is.na( plot_ece_meta_selSEM ) )


caitlin_resist <- psem(
  lm(log10(resistance) ~ richness + dominant_relative_abund_zero + nut_dummy,
     data = plot_ece_meta_selSEM ),
  lm(richness ~ nut_dummy,
     data = plot_ece_meta_selSEM),
  lm(dominant_relative_abund_zero ~nut_dummy ,
     data = plot_ece_meta_selSEM), 
  
  lm(log10(resilience) ~ richness + dominant_relative_abund_zero + nut_dummy,
     data = plot_ece_meta_selSEM ) ,
  data = plot_ece_meta_selSEM
)

summary(caitlin_resist)
plot(caitlin_resist)
multigroup(caitlin_resist, group = "spei9_category")

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
