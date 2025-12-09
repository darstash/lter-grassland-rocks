# TITLE:        LTER Grassland Rock: Create core analyses 
# AUTHORS:      Ashley Darst, Joshua Ajowele
# COLLABORATORS:  
# DATA INPUT:   Data imported as csv files from shared Google drive L2 folder
# DATA OUTPUT:  Core analyses
# PROJECT:      LTER Grassland Rock
# DATE:         October 2024 , last updated: December 2025

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
library(cowplot)
# Set working directory 
L0_dir <- Sys.getenv("L0DIR")
L1_dir <- Sys.getenv("L1DIR")
L2_dir <- Sys.getenv("L2DIR")
list.files(L2_dir)

# Read in CSV files
plot <- read.csv(file.path(L2_dir, "plot_metrics_SPEI_diversity_L2.csv"))
meta <- read.csv(file.path(L1_dir, "metadata_L1.csv"))
ece_9_norm<-read.csv(file.path(L2_dir, "ece_resist_resil_spei9_norm_L2.csv"))
# Only keep distinct rows in ece
ece_9_norm<-distinct(ece_9_norm)
# Change ex_year to year
ece_9_norm<-ece_9_norm%>%
  rename(year = ex_year,
         resistance_n=resistance,
         resilience_n=resilience)
  

# Merge plot with resistance and resilience
plot_ece_9<-left_join(plot, ece_9_norm)

# Standardize column names
plot_ece_9 <- clean_names(plot_ece_9)


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
plot_ece_9_meta <- left_join(plot_ece_9, meta)
unique(plot_ece_9_meta$experiment)


# Remove NAs for non-extreme years
plot_ece_9_rm_na <- plot_ece_9_meta %>%
  drop_na(resistance_n)
#checking each experiment/study
plot_filter<-plot_ece_9_rm_na%>%
  filter(experiment=="004b_fl")

# Make year a factor
str(plot_ece_9_rm_na)
plot_ece_9_rm_na$year <- as.factor(plot_ece_9_rm_na$year)

#subset to have control and nitrogen (including nitrogen in combination with other nutrients) treated plots
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
         )%>%
  mutate(extreme_after=case_when(after_year_type=="Extreme wet"~"yes",
         after_year_type=="Extreme dry"~"yes",
         .default="no"),
         nitrogen=case_when(nitrogen=="no_fertilizer"~"ambient",
                            nitrogen=="N"~"nutrients"))#when the year after the extreme event is another extreme event

#check treatment groups
unique(plot_ece_9_cn$treatment)
unique(plot_ece_9_cn$nitrogen)
table(plot_ece_9_cn$extreme_after)
table(plot_ece_9_cn$nutrients_grouped)
table(plot_ece_9_cn$nitrogen)


#rename prior year community metrics and create a new dataframe
plot_ece_9_cn_prior<-plot_ece_9_cn%>%
  select(-richness, -evar, -dominant_relative_abund, -dominant_relative_abund_zero,-eq)%>%
  rename(richness=prior_year_rich,
         dominant_relative_abund_zero= prior_year_dom_zero,
         dominant_relative_abund=prior_year_dom,
         evar=prior_year_evar,
         eq=prior_year_eq)

##resistance with control and nitrogen####
resis_norm<-lmer(log(resistance_n)~richness*dominant_relative_abund_zero+nitrogen+richness:nitrogen+evar+
                 dominant_relative_abund_zero:nitrogen+richness:spei9_category+dominant_relative_abund_zero:spei9_category+
                 evar:nitrogen+evar:spei9_category+spei9_category+nitrogen:spei9_category+(1|site/experiment/uniqueid)
               +(1|year), data=plot_ece_9_cn_prior, REML=F)
summary(resis_norm)#singular due to site having zero variance, left it in for consistency since result same with or without
anova(resis_norm)
#same model but with site removed from random effect structure
# resis_norm_site<-lmer(log(resistance_n)~richness*dominant_relative_abund_zero+nitrogen+richness:nitrogen+evar+
#                    dominant_relative_abund_zero:nitrogen+richness:spei9_category+dominant_relative_abund_zero:spei9_category+
#                    evar:nitrogen+evar:spei9_category+spei9_category+nitrogen:spei9_category+(1|experiment/uniqueid)
#                  +(1|year), data=plot_ece_9_cn_prior, REML=F)
# summary(resis_norm_site)
#model update
resis_norm1<-lmer(log(resistance_n)~richness*dominant_relative_abund_zero+nitrogen+richness:nitrogen+evar+
                   dominant_relative_abund_zero:nitrogen+richness:spei9_category+dominant_relative_abund_zero:spei9_category+
                   evar:spei9_category+spei9_category+nitrogen:spei9_category+(1|site/experiment/uniqueid)
                 +(1|year), data=plot_ece_9_cn_prior, REML=F)
anova(resis_norm1)
#model update
resis_norm2<-lmer(log(resistance_n)~richness*dominant_relative_abund_zero+nitrogen+richness:nitrogen+evar+
                    +richness:spei9_category+dominant_relative_abund_zero:spei9_category+
                    evar:spei9_category+spei9_category+nitrogen:spei9_category+(1|site/experiment/uniqueid)
                  +(1|year), data=plot_ece_9_cn_prior, REML=F)
#model comparison with loglik
anova(resis_norm2, resis_norm1)
#model update
resis_norm3<-lmer(log(resistance_n)~richness*dominant_relative_abund_zero+nitrogen+richness:nitrogen+evar+
                    +richness:spei9_category+
                    evar:spei9_category+spei9_category+nitrogen:spei9_category+(1|site/experiment/uniqueid)
                  +(1|year), data=plot_ece_9_cn_prior, REML=F)
#model comparison
anova(resis_norm3,resis_norm2)
#model update
resis_norm4<-lmer(log(resistance_n)~richness+dominant_relative_abund_zero+nitrogen+richness:nitrogen+evar+
                    +richness:spei9_category+
                    evar:spei9_category+spei9_category+nitrogen:spei9_category+(1|site/experiment/uniqueid)
                  +(1|year), data=plot_ece_9_cn_prior, REML=F)
#model comparison
anova(resis_norm4,resis_norm3)
#model update
resis_norm5<-lmer(log(resistance_n)~richness+dominant_relative_abund_zero+nitrogen+richness:nitrogen+evar+
                    +richness:spei9_category+
                    evar:spei9_category+spei9_category+(1|site/experiment/uniqueid)
                  +(1|year), data=plot_ece_9_cn_prior, REML=F)
#model comparison
anova(resis_norm5, resis_norm4)
#model update
resis_norm6<-lmer(log(resistance_n)~richness+dominant_relative_abund_zero+nitrogen+richness:nitrogen+evar+
                    +richness:spei9_category+
                    spei9_category+(1|site/experiment/uniqueid)
                  +(1|year), data=plot_ece_9_cn_prior, REML=F)
#model comparison
anova(resis_norm6,resis_norm5)
#model update
resis_norm7<-lmer(log(resistance_n)~richness+dominant_relative_abund_zero+nitrogen+evar+
                    +richness:spei9_category+
                    spei9_category+(1|site/experiment/uniqueid)
                  +(1|year), data=plot_ece_9_cn_prior, REML=F)
#model comparison
anova(resis_norm7, resis_norm6)
#model selection with log likelihood ratio
anova(resis_norm7,resis_norm6,resis_norm5,resis_norm4,resis_norm3,resis_norm2,resis_norm1,resis_norm)
#refit best model with REML
resis_norm8<-lmer(log(resistance_n)~richness+dominant_relative_abund_zero+nitrogen+evar+
                    +richness:spei9_category+
                    spei9_category+(1|site/experiment/uniqueid)
                  +(1|year), data=plot_ece_9_cn_prior)
anova(resis_norm8)
summary(resis_norm8)
#model diagnostics
simres <- simulateResiduals(resis_norm8)
plot(simres)
resis8<-plot(check_model(resis_norm8))
resis8[[5]]
###figure for full resistance model####
resis_full_std <- update(resis_norm8, 
                             data = plot_ece_9_cn_prior %>% 
                               mutate(richness = scale(richness),
                                      evar = scale(evar),
                                      dominant_relative_abund_zero = scale(dominant_relative_abund_zero),
                               ))

resis_full_estim <-coef(summary(resis_full_std)) %>%
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
                                      "nitrogennutrients",
                                      "spei9_categoryExtreme wet",
                                      "richness:spei9_categoryExtreme wet"
                           )),
         # cosmetics in the names
         Variable_labels = factor(case_when(Variable %in% "(Intercept)" ~ "intercept",
                                            Variable %in% "evar" ~ "evenness",
                                            Variable %in% "dominant_relative_abund_zero" ~ "dominant",
                                            Variable %in% "nitrogennutrients" ~ "nutrients",
                                            Variable %in% "spei9_categoryExtreme wet" ~ "extreme wet",
                                            Variable %in% "richness:spei9_categoryExtreme wet" ~ "richness:extreme wet",
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
  ggplot(aes(y = Variable_labels, x = Estimate))+
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  labs(title = "ln(resistance)",
       x = "estimate \u00B1 se") +
  scale_y_discrete(limits = rev) +
  geom_vline(xintercept = 0,linetype=2) +
  geom_point(col = "#A020F0", size =3)+
  geom_errorbarh(aes(xmin = Estimate-SE, xmax = Estimate+SE), height = 0.2,col = "#A020F0") +
  geom_text(aes(x = stars_location, label = stars))

#view figure
print(resis_full_estim)
#standardized coefficient for result table
summary(resis_full_std)


###examine wet and dry separately based on the combined model####
plot_ece_9_cn_prior_wet<-plot_ece_9_cn_prior%>%
  filter(spei9_category=="Extreme wet")
plot_ece_9_cn_prior_dry<-plot_ece_9_cn_prior%>%
  filter(spei9_category=="Extreme dry")

#model for extreme wet
resis_norm_wet<-lmer(log(resistance_n)~richness+dominant_relative_abund_zero+nitrogen+evar+
                    (1|site/experiment/uniqueid)
                  +(1|year), data=plot_ece_9_cn_prior_wet)
anova(resis_norm_wet)#warning is due to site variance being low. similar output without convergence warning when site is removed
summary(resis_norm_wet)
simres <- simulateResiduals(resis_norm_wet)
plot(simres)
check_model(resis_norm_wet)

#model for extreme dry
resis_norm_dry<-lmer(log(resistance_n)~richness+dominant_relative_abund_zero+nitrogen+evar+
                       (1|site/experiment/uniqueid)
                     +(1|year), data=plot_ece_9_cn_prior_dry)
anova(resis_norm_dry)
summary(resis_norm_dry)
simres <- simulateResiduals(resis_norm_dry)
plot(simres)
check_model(resis_norm_dry)

#refit with EQ as evenness metric
#resis_wet_EQ<-lmer(log(resistance_n)~richness+dominant_relative_abund_zero+nitrogen+eq+
#                     (1|site/experiment/uniqueid)
#                   +(1|year), data=plot_ece_9_cn_prior_wet)
#anova(resis_wet_EQ)#did not change model result
#resis_dry_EQ<-lmer(log(resistance_n)~richness+dominant_relative_abund_zero+nitrogen+eq+
#                       (1|site/experiment/uniqueid)
#                     +(1|year), data=plot_ece_9_cn_prior_dry)
#anova(resis_dry_EQ)
#####create figure for resistance extreme dry and extreme wet####
resis_norm_wet_std <- update(resis_norm_wet, 
                        data = plot_ece_9_cn_prior_wet %>% 
                          mutate(richness = scale(richness),
                                 evar = scale(evar),
                                 dominant_relative_abund_zero = scale(dominant_relative_abund_zero),
                          ))

resis_norm_estim_wet <-coef(summary(resis_norm_wet_std)) %>%
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
                                      "nitrogennutrients"
                           )),
         # cosmetics in the names
         Variable_labels = factor(case_when(Variable %in% "(Intercept)" ~ "intercept",
                                            Variable %in% "evar" ~ "evenness",
                                            Variable %in% "dominant_relative_abund_zero" ~ "dominant",
                                            Variable %in% "nitrogennutrients" ~ "nutrients",
                                            .default = Variable)),
         Variable_labels = fct_reorder(Variable_labels, as.numeric(Variable)),
         # get significances and their plotting location (perhaps necessary to play with the +/- offset)
         stars = case_when(p < 0.001 ~ "***",
                           p > 0.001 & p <0.01 ~ "**",
                           p > 0.01  & p < 0.05 ~ "*",
                           p > 0.05 & p <0.1 ~ ".",
                           .default = ""),
         stars_location = case_when(Estimate < 0 ~ Estimate - SE - 0.01,
                                    Estimate > 0 ~ Estimate + SE + 0.01)
  ) %>%
  filter(Variable!="(Intercept)")%>%
  mutate(event="wet")%>%
  
  ggplot(aes(y = Variable_labels, x = Estimate, fill=event))+
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  labs(title = "ln(resistance)",
       x = "estimate \u00B1 se") +
  scale_y_discrete(limits = rev) +
  geom_vline(xintercept = 0, linetype=2) +
  geom_point(col = "#377EB8", size =3)+
  geom_errorbarh(aes(xmin = Estimate-SE, xmax = Estimate+SE), height = 0.2,col = "#377EB8") +
  geom_text(aes(x = stars_location, label = stars))
#standardized estimates
summary(resis_norm_wet_std)
#view
print(resis_norm_estim_wet)
  
#dry
resis_norm_dry_std <- update(resis_norm_dry, 
                             data = plot_ece_9_cn_prior_dry %>% 
                               mutate(richness = scale(richness),
                                      evar = scale(evar),
                                      dominant_relative_abund_zero = scale(dominant_relative_abund_zero),
                               ))

resis_norm_estim_dry <-coef(summary(resis_norm_dry_std)) %>%
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
                                      "nitrogennutrients"
                           )),
         # cosmetics in the names
         Variable_labels = factor(case_when(Variable %in% "(Intercept)" ~ "intercept",
                                            Variable %in% "evar" ~ "evenness",
                                            Variable %in% "dominant_relative_abund_zero" ~ "dominant",
                                            Variable %in% "nitrogennutrients" ~ "nutrients",
                                            .default = Variable)),
         Variable_labels = fct_reorder(Variable_labels, as.numeric(Variable)),
         # get significances and their plotting location (perhaps necessary to play with the +/- offset)
         stars = case_when(p < 0.001 ~ "***",
                           p > 0.001 & p <0.01 ~ "**",
                           p > 0.01  & p < 0.05 ~ "*",
                           p > 0.05 & p <0.1 ~ ".",
                           .default = ""),
         stars_location = case_when(Estimate < 0 ~ Estimate - SE - 0.01,
                                    Estimate > 0 ~ Estimate + SE + 0.01)
  ) %>%
  filter(Variable!="(Intercept)")%>%
  mutate(event="dry")%>%
  
  ggplot(aes(y = Variable_labels, x = Estimate, fill=event))+
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  labs(title = "ln(resistance)",
       x = "estimate \u00B1 se") +
  scale_y_discrete(limits = rev) +
  geom_vline(xintercept = 0, linetype=2) +
  geom_point(col = "#E41A1C", size =3)+
  geom_errorbarh(aes(xmin = Estimate-SE, xmax = Estimate+SE), height = 0.2,col = "#E41A1C") +
  geom_text(aes(x = stars_location, label = stars))
#standardized estimates 
summary(resis_norm_dry_std)
#view
print(resis_norm_estim_dry)

####Sensitivity analysis for resistance####
# Sensitivity analysis dry (leave-one-site out) 
sites <- unique(plot_ece_9_cn_prior_dry$site)
results_list <- list()

for (s in sites) {
  cat("Leaving out site:", s, "\n")
  # remove one site
  df_subset <- plot_ece_9_cn_prior_dry %>% filter(site != s)
  # refit model
  mod <- lmer(log(resistance_n)~richness+dominant_relative_abund_zero+nitrogen+evar+
                (1|site/experiment/uniqueid)
              +(1|year), data = df_subset)
  
  # store the model summary (you can change this)
  results_list[[as.character(s)]] <- summary(mod)
}

# Extract fixed effects from the list of summaries
coef_df <- map_df(
  names(results_list),
  ~ {
    sm <- results_list[[.x]]
    fe <- as.data.frame(sm$coefficients)
    fe$term <- rownames(fe)
    fe$site_left_out <- .x
    fe
  }
)
coef_df<-coef_df%>%
  mutate(term=case_when(term=="dominant_relative_abund_zero"~"dominant",
                        term=="nitrogennutrients"~"nutrients",
                        term=="evar"~"evenness",
                        .default=term))
resist_dry_sens<-ggplot(coef_df, aes(x = site_left_out, y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = Estimate - `Std. Error`,
                    ymax = Estimate + `Std. Error`),
                width = 0.2) +
  facet_wrap(~ term, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(
    title = "sensitivity analysis of resistance to dry extreme events",
    x = "Site left out",
    y = "ln(Resistance ± SE)"
  )

# Sensitivity analysis resistance wet (leave-one-site out) 
sites <- unique(plot_ece_9_cn_prior_wet$site)
results_list <- list()

for (s in sites) {
  cat("Leaving out site:", s, "\n")
  # remove one site
  df_subset <- plot_ece_9_cn_prior_wet %>% filter(site != s)
  # refit model
  mod <- lmer(log(resistance_n)~richness+dominant_relative_abund_zero+nitrogen+evar+
                (1|site/experiment/uniqueid)
              +(1|year), data = df_subset)
  
  # store the model summary (you can change this)
  results_list[[as.character(s)]] <- summary(mod)
}

# Extract fixed effects from the list of summaries
coef_df <- map_df(
  names(results_list),
  ~ {
    sm <- results_list[[.x]]
    fe <- as.data.frame(sm$coefficients)
    fe$term <- rownames(fe)
    fe$site_left_out <- .x
    fe
  }
)
coef_df<-coef_df%>%
  mutate(term=case_when(term=="dominant_relative_abund_zero"~"dominant",
                        term=="nitrogennutrients"~"nutrients",
                        term=="evar"~"evenness",
                        .default=term))

resist_wet_sens<-ggplot(coef_df, aes(x = site_left_out, y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = Estimate - `Std. Error`,
                    ymax = Estimate + `Std. Error`),
                width = 0.2) +
  facet_wrap(~ term, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(
    title = "Sensitivity analysis of resistance to wet extreme events",
    x = "Site left out",
    y = "ln(Resistance ± SE)"
  )

######sensitivity anlaysis by year####
#resistance to dry
years <- unique(plot_ece_9_cn_prior_dry$year)
results_list <- list()

for (s in years) {
  cat("Leaving out year:", s, "\n")
  # remove one site
  df_subset <- plot_ece_9_cn_prior_dry %>% filter(year != s)
  # refit model
  mod <- lmer(log(resistance_n)~richness+dominant_relative_abund_zero+nitrogen+evar+
                (1|site/experiment/uniqueid)
              +(1|year), data = df_subset)
  
  # store the model summary (you can change this)
  results_list[[as.character(s)]] <- summary(mod)
}

# Extract fixed effects from the list of summaries
coef_df <- map_df(
  names(results_list),
  ~ {
    sm <- results_list[[.x]]
    fe <- as.data.frame(sm$coefficients)
    fe$term <- rownames(fe)
    fe$year_left_out <- .x
    fe
  }
)
coef_df<-coef_df%>%
  mutate(term=case_when(term=="dominant_relative_abund_zero"~"dominant",
                        term=="nitrogennutrients"~"nutrients",
                        term=="evar"~"evenness",
                        .default=term))
year_resist_dry_sens<-ggplot(coef_df, aes(x = year_left_out, y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = Estimate - `Std. Error`,
                    ymax = Estimate + `Std. Error`),
                width = 0.2) +
  facet_wrap(~ term, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(
    title = "sensitivity analysis of resistance to dry extreme events",
    x = "year left out",
    y = "ln(Resistance ± SE)"
  )
#resistance to extreme wet
years <- unique(plot_ece_9_cn_prior_wet$year)
results_list <- list()

for (s in years) {
  cat("Leaving out year:", s, "\n")
  # remove one site
  df_subset <- plot_ece_9_cn_prior_wet %>% filter(year != s)
  # refit model
  mod <- lmer(log(resistance_n)~richness+dominant_relative_abund_zero+nitrogen+evar+
                (1|site/experiment/uniqueid)
              +(1|year), data = df_subset)
  
  # store the model summary (you can change this)
  results_list[[as.character(s)]] <- summary(mod)
}

# Extract fixed effects from the list of summaries
coef_df <- map_df(
  names(results_list),
  ~ {
    sm <- results_list[[.x]]
    fe <- as.data.frame(sm$coefficients)
    fe$term <- rownames(fe)
    fe$year_left_out <- .x
    fe
  }
)
coef_df<-coef_df%>%
  mutate(term=case_when(term=="dominant_relative_abund_zero"~"dominant",
                        term=="nitrogennutrients"~"nutrients",
                        term=="evar"~"evenness",
                        .default=term))
year_resist_wet_sens<-ggplot(coef_df, aes(x = year_left_out, y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = Estimate - `Std. Error`,
                    ymax = Estimate + `Std. Error`),
                width = 0.2) +
  facet_wrap(~ term, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(
    title = "sensitivity analysis of resistance to wet extreme events",
    x = "year left out",
    y = "ln(Resistance ± SE)"
  )


#remove resilience values where extreme event occured after an extreme event####
plot_ece_9_cn_prior_rm<-plot_ece_9_cn_prior%>%
  filter(extreme_after=="no")

###resilience analysis after removing resilience in years with extreme events after an extreme event####
resil_norm<-lmer(log(resilience_n)~richness*dominant_relative_abund_zero+nitrogen+richness:nitrogen+evar+
                   dominant_relative_abund_zero:nitrogen+richness:spei9_category+dominant_relative_abund_zero:spei9_category+
                   evar:nitrogen+evar:spei9_category+spei9_category+nitrogen:spei9_category+(1|site/experiment/uniqueid)
                 +(1|year), data=plot_ece_9_cn_prior_rm, REML=F)
anova(resil_norm)
#model update
resil_norm1<-lmer(log(resilience_n)~richness*dominant_relative_abund_zero+nitrogen+evar+
                   dominant_relative_abund_zero:nitrogen+richness:spei9_category+dominant_relative_abund_zero:spei9_category+
                   evar:nitrogen+evar:spei9_category+spei9_category+nitrogen:spei9_category+(1|site/experiment/uniqueid)
                 +(1|year), data=plot_ece_9_cn_prior_rm, REML=F)
#model comparison with loglik
anova(resil_norm1, resil_norm)
#model update
resil_norm2<-lmer(log(resilience_n)~richness*dominant_relative_abund_zero+nitrogen+evar+
                    dominant_relative_abund_zero:nitrogen+richness:spei9_category+dominant_relative_abund_zero:spei9_category+
                    evar:nitrogen+spei9_category+nitrogen:spei9_category+(1|site/experiment/uniqueid)
                  +(1|year), data=plot_ece_9_cn_prior_rm, REML=F)
#mode, comparison
anova(resil_norm2, resil_norm1)
#model update
resil_norm3<-lmer(log(resilience_n)~richness+dominant_relative_abund_zero+nitrogen+evar+
                    dominant_relative_abund_zero:nitrogen+richness:spei9_category+dominant_relative_abund_zero:spei9_category+
                    evar:nitrogen+spei9_category+nitrogen:spei9_category+(1|site/experiment/uniqueid)
                  +(1|year), data=plot_ece_9_cn_prior_rm, REML=F)
#model comparison
anova(resil_norm3, resil_norm2)

#model update
resil_norm4<-lmer(log(resilience_n)~richness+dominant_relative_abund_zero+nitrogen+evar+
                    dominant_relative_abund_zero:nitrogen+richness:spei9_category+
                    evar:nitrogen+spei9_category+nitrogen:spei9_category+(1|site/experiment/uniqueid)
                  +(1|year), data=plot_ece_9_cn_prior_rm, REML=F)
#model comparison
anova(resil_norm4, resil_norm3)

#model update
resil_norm5<-lmer(log(resilience_n)~richness+dominant_relative_abund_zero+nitrogen+evar+
                    dominant_relative_abund_zero:nitrogen+richness:spei9_category+
                    evar:nitrogen+spei9_category+(1|site/experiment/uniqueid)
                  +(1|year), data=plot_ece_9_cn_prior_rm, REML=F)
#model comparison
anova(resil_norm5, resil_norm4)
#model selection
anova(resil_norm5,resil_norm4,resil_norm3,resil_norm2,resil_norm1,resil_norm)

#refit best model with REML
resil_norm6<-lmer(log(resilience_n)~richness+dominant_relative_abund_zero+nitrogen+evar+
                    dominant_relative_abund_zero:nitrogen+richness:spei9_category+
                    evar:nitrogen+spei9_category+(1|site/experiment/uniqueid)
                  +(1|year), data=plot_ece_9_cn_prior_rm)
anova(resil_norm6)
summary(resil_norm6)
simres <- simulateResiduals(resil_norm6)
plot(simres)

####figure for full resilience model####
resil_full_std <- update(resil_norm6, 
                             data = plot_ece_9_cn_prior_rm %>% 
                               mutate(richness = scale(richness),
                                      evar = scale(evar),
                                      dominant_relative_abund_zero = scale(dominant_relative_abund_zero),
                               ))

resil_full_estim <-coef(summary(resil_full_std)) %>%
  data.frame()%>%
  rownames_to_column("Variable") %>%
  rename("SE" = "Std..Error",
         "p" = "Pr...t..") %>%
  # order things neatly
  mutate(Variable = factor(Variable,
                           levels = c("(Intercept)",
                                      "richness",
                                      "dominant_relative_abund_zero",
                                      "evar",
                                      "nitrogennutrients",
                                      "spei9_categoryExtreme wet",
                                      "richness:spei9_categoryExtreme wet",
                                      "dominant_relative_abund_zero:nitrogennutrients",
                                      "nitrogennutrients:evar"
                           )),
         # cosmetics in the names
         Variable_labels = factor(case_when(Variable %in% "(Intercept)" ~ "intercept",
                                            Variable %in% "evar" ~ "evenness",
                                            Variable %in% "dominant_relative_abund_zero" ~ "dominant",
                                            Variable %in% "nitrogennutrients" ~ "nutrients",
                                            Variable %in% "spei9_categoryExtreme wet" ~ "extreme wet",
                                            Variable %in% "richness:spei9_categoryExtreme wet" ~ "richness:extreme wet",
                                            Variable %in% "dominant_relative_abund_zero:nitrogennutrients" ~ "dominant:nutrients",
                                            Variable %in% "nitrogennutrients:evar" ~ "evenness:nutrients",
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
  
  ggplot(aes(y = Variable_labels, x = Estimate))+
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  labs(title = "ln(resilience)",
       x = "estimate \u00B1 se") +
  scale_y_discrete(limits = rev) +
  geom_vline(xintercept = 0, linetype=2) +
  geom_point(col = "#A020F0", size =3)+
  geom_errorbarh(aes(xmin = Estimate-SE, xmax = Estimate+SE), height = 0.2,col = "#A020F0") +
  geom_text(aes(x = stars_location, label = stars))
#standardized estimates
summary(resil_full_std)
#view
print(resil_full_estim)
#####combine full model figure####
resis_full_estim+resil_full_estim& plot_annotation(tag_levels = 'A')
####split into wet and dry####
plot_ece_9_cn_prior_rm_wet<-plot_ece_9_cn_prior_rm%>%
  filter(spei9_category=="Extreme wet")
plot_ece_9_cn_prior_rm_dry<-plot_ece_9_cn_prior_rm%>%
  filter(spei9_category=="Extreme dry")

#model for extreme wet
resil_norm_wet<-lmer(log(resilience_n)~richness+dominant_relative_abund_zero+nitrogen+evar+
                    dominant_relative_abund_zero:nitrogen+
                    evar:nitrogen+(1|site/experiment/uniqueid)
                  +(1|year), data=plot_ece_9_cn_prior_rm_wet)
summary(resil_norm_wet)
simres <- simulateResiduals(resil_norm_wet)
plot(simres)
check_model(resil_norm_wet)
#model for extreme dry
resil_norm_dry<-lmer(log(resilience_n)~richness+dominant_relative_abund_zero+nitrogen+evar+
                       dominant_relative_abund_zero:nitrogen+
                       evar:nitrogen+(1|site/experiment/uniqueid)
                     +(1|year), data=plot_ece_9_cn_prior_rm_dry)
summary(resil_norm_dry)
simres <- simulateResiduals(resil_norm_dry)
plot(simres)
check_model(resil_norm_dry)

#refit with EQ as evenness
#resil_wet_eq<-lmer(log(resilience_n)~richness+dominant_relative_abund_zero+nitrogen+eq+
#                       dominant_relative_abund_zero:nitrogen+
#                       eq:nitrogen+(1|site/experiment/uniqueid)
#                     +(1|year), data=plot_ece_9_cn_prior_rm_wet)
#anova(resil_wet_eq)
#resil_dry_eq<-lmer(log(resilience_n)~richness+dominant_relative_abund_zero+nitrogen+eq+
#                       dominant_relative_abund_zero:nitrogen+
#                       eq:nitrogen+(1|site/experiment/uniqueid)
#                     +(1|year), data=plot_ece_9_cn_prior_rm_dry)
#summary(resil_dry_eq)
#largely similar
####create figure for resiience wet and dry####
resil_norm_wet_std <- update(resil_norm_wet, 
                             data = plot_ece_9_cn_prior_rm_wet %>% 
                               mutate(richness = scale(richness),
                                      evar = scale(evar),
                                      dominant_relative_abund_zero = scale(dominant_relative_abund_zero),
                               ))

resil_norm_estim_wet <-coef(summary(resil_norm_wet_std)) %>%
  data.frame()%>%
  rownames_to_column("Variable") %>%
  rename("SE" = "Std..Error",
         "p" = "Pr...t..") %>%
  # order things neatly
  mutate(Variable = factor(Variable,
                           levels = c("(Intercept)",
                                      "richness",
                                      "dominant_relative_abund_zero",
                                      "evar",
                                      "nitrogennutrients",
                                      "dominant_relative_abund_zero:nitrogennutrients",
                                      "nitrogennutrients:evar"
                           )),
         # cosmetics in the names
         Variable_labels = factor(case_when(Variable %in% "(Intercept)" ~ "intercept",
                                            Variable %in% "evar" ~ "evenness",
                                            Variable %in% "dominant_relative_abund_zero" ~ "dominant",
                                            Variable %in% "nitrogennutrients" ~ "nutrients",
                                            Variable %in% "dominant_relative_abund_zero:nitrogennutrients" ~ "dominant:nutrients",
                                            Variable %in% "nitrogennutrients:evar" ~ "evenness:nutrients",
                                            .default = Variable)),
         Variable_labels = fct_reorder(Variable_labels, as.numeric(Variable)),
         # get significances and their plotting location (perhaps necessary to play with the +/- offset)
         stars = case_when(p < 0.001 ~ "***",
                           p > 0.001 & p <0.01 ~ "**",
                           p > 0.01  & p < 0.05 ~ "*",
                           p > 0.05 & p <0.1 ~ ".",
                           .default = ""),
         stars_location = case_when(Estimate < 0 ~ Estimate - SE - 0.03,
                                    Estimate > 0 ~ Estimate + SE + 0.03)
  ) %>%
  filter(Variable!="(Intercept)")%>%
  
  ggplot(aes(y = Variable_labels, x = Estimate))+
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  labs(title = "ln(resilience)",
       x = "estimate \u00B1 se") +
  scale_y_discrete(limits = rev) +
  geom_vline(xintercept = 0, linetype=2) +
  geom_point(col = "#377EB8", size =3)+
  geom_errorbarh(aes(xmin = Estimate-SE, xmax = Estimate+SE), height = 0.2,col = "#377EB8") +
  geom_text(aes(x = stars_location, label = stars))
#standardized coefficient
summary(resil_norm_wet_std)
#view
print(resil_norm_estim_wet)
#dry
resil_norm_dry_std <- update(resil_norm_dry, 
                             data = plot_ece_9_cn_prior_rm_dry %>% 
                               mutate(richness = scale(richness),
                                      evar = scale(evar),
                                      dominant_relative_abund_zero = scale(dominant_relative_abund_zero),
                               ))

resil_norm_estim_dry <-coef(summary(resil_norm_dry_std)) %>%
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
                                      "nitrogennutrients",
                                      "dominant_relative_abund_zero:nitrogennutrients",
                                      "nitrogennutrients:evar"
                           )),
         # cosmetics in the names
         Variable_labels = factor(case_when(Variable %in% "(Intercept)" ~ "intercept",
                                            Variable %in% "evar" ~ "evenness",
                                            Variable %in% "dominant_relative_abund_zero" ~ "dominant",
                                            Variable %in% "nitrogennutrients" ~ "nutrients",
                                            Variable %in% "dominant_relative_abund_zero:nitrogennutrients" ~ "dominant:nutrients",
                                            Variable %in% "nitrogennutrients:evar" ~ "evenness:nutrients",
                                            .default = Variable)),
         Variable_labels = fct_reorder(Variable_labels, as.numeric(Variable)),
         # get significances and their plotting location (perhaps necessary to play with the +/- offset)
         stars = case_when(p < 0.001 ~ "***",
                           p > 0.001 & p <0.01 ~ "**",
                           p > 0.01  & p < 0.05 ~ "*",
                           p > 0.05 & p <0.1 ~ ".",
                           .default = ""),
         stars_location = case_when(Estimate < 0 ~ Estimate - SE - 0.03,
                                    Estimate > 0 ~ Estimate + SE + 0.03)
  ) %>%
  filter(Variable!="(Intercept)")%>%
  
  ggplot(aes(y = Variable_labels, x = Estimate))+
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  labs(title = "ln(resilience)",
       x = "estimate \u00B1 se") +
  scale_y_discrete(limits = rev) +
  geom_vline(xintercept = 0, linetype=2) +
  geom_point(col = "#E41A1C", size =3)+
  geom_errorbarh(aes(xmin = Estimate-SE, xmax = Estimate+SE), height = 0.2,col = "#E41A1C") +
  geom_text(aes(x = stars_location, label = stars))
#standardized estimates 
summary(resil_norm_dry_std)
#view
print(resil_norm_estim_dry)
##combine figure for resistance and resilience wet and dry####
#using patchwork
resis_norm_estim_dry + resis_norm_estim_wet+ plot_layout(guides = "collect")+resil_norm_estim_dry+resil_norm_estim_wet& plot_annotation(tag_levels = 'A')&theme(legend.position = "bottom")
#&scale_fill_discrete(limits =c(resis_norm_estim_dry$event_type, resis_norm_estim_wet$event_type))

####Interaction figures####
ggpredict(model = resil_norm_dry_std, terms = c("evar", "nitrogen"), back_transform = F) %>%
  plot(show_data = F)+
  labs(title=NULL, x="Evenness", y="ln(Resilience)")+
  scale_fill_manual(values=c("#377EB8","#E41A1C"))+
  scale_color_manual(values=c("#377EB8","#E41A1C"),labels = c("control", "nutrients"))+
  
  theme(legend.title=element_blank())#+
  theme_classic()
ggpredict(model = resil_norm_wet_std, terms =  c("dominant_relative_abund_zero", "nitrogen"), back_transform = F) %>%
  plot(show_data = F,colors= c("#377EB8","#E41A1C"))
ggpredic
?labs()
####sensitivity analyses#####
# Sensitivity analysis resilience wet (leave-one-site out) 
sites <- unique(plot_ece_9_cn_prior_rm_wet$site)
results_list <- list()

for (s in sites) {
  cat("Leaving out site:", s, "\n")
  # remove one site
  df_subset <- plot_ece_9_cn_prior_rm_wet %>% filter(site != s)
  # refit model
  mod <- lmer(log(resilience_n)~richness+dominant_relative_abund_zero+nitrogen+evar+
                dominant_relative_abund_zero:nitrogen+
                evar:nitrogen+(1|site/experiment/uniqueid)
              +(1|year), data = df_subset)
  
  # store the model summary (you can change this)
  results_list[[as.character(s)]] <- summary(mod)
}

# Extract fixed effects from the list of summaries
coef_df <- map_df(
  names(results_list),
  ~ {
    sm <- results_list[[.x]]
    fe <- as.data.frame(sm$coefficients)
    fe$term <- rownames(fe)
    fe$site_left_out <- .x
    fe
  }
)
#rename terms appropriately
coef_df<-coef_df%>%
  mutate(term=case_when(term=="dominant_relative_abund_zero"~"dominant",
                        term=="nitrogennutrients"~"nutrients",
                        term=="evar"~"evenness",
                        term %in% "nitrogennutrients" ~ "nutrients",
                        term %in% "dominant_relative_abund_zero:nitrogennutrients" ~ "dominant:nutrients",
                        term %in% "nitrogennutrients:evar" ~ "evenness:nutrients",
                        .default=term))
resil_wet_sens<-ggplot(coef_df, aes(x = site_left_out, y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = Estimate - `Std. Error`,
                    ymax = Estimate + `Std. Error`),
                width = 0.2) +
  facet_wrap(~ term, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(
    title = "sensitivity analysis of resilience to wet extreme events",
    x = "Site left out",
    y = "ln(Resilience ± SE)"
  )


# Sensitivity analysis resilience dry (leave-one-site out) 
sites <- unique(plot_ece_9_cn_prior_rm_dry$site)
results_list <- list()
#did not work-probably because not all site had all the required interactions
for (s in sites) {
  cat("Leaving out site:", s, "\n")
  # remove one site
  df_subset <- plot_ece_9_cn_prior_rm_dry %>% filter(site != s)
  # refit model
  mod <- lmer(log(resilience_n)~richness+dominant_relative_abund_zero+nitrogen+evar+
                dominant_relative_abund_zero:nitrogen+
                evar:nitrogen+(1|site/experiment/uniqueid)
              +(1|year), data = df_subset)
  
  # store the model summary (you can change this)
  results_list[[as.character(s)]] <- summary(mod)
}

# Extract fixed effects from the list of summaries
coef_df <- map_df(
  names(results_list),
  ~ {
    sm <- results_list[[.x]]
    fe <- as.data.frame(sm$coefficients)
    fe$term <- rownames(fe)
    fe$site_left_out <- .x
    fe
  }
)
ggplot(coef_df, aes(x = site_left_out, y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = Estimate - `Std. Error`,
                    ymax = Estimate + `Std. Error`),
                width = 0.2) +
  facet_wrap(~ term, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(
    title = "Leave-one-site-out sensitivity analysis",
    x = "Site left out",
    y = "In(Resilience_dry ± SE)"
  )
######sensitivity analysis by year####
# Sensitivity analysis resilience wet (leave-one-year out) 
years <- unique(plot_ece_9_cn_prior_rm_wet$year)
results_list <- list()

for (s in years) {
  cat("Leaving out year:", s, "\n")
  # remove one site
  df_subset <- plot_ece_9_cn_prior_rm_wet %>% filter(year != s)
  # refit model
  mod <- lmer(log(resilience_n)~richness+dominant_relative_abund_zero+nitrogen+evar+
                dominant_relative_abund_zero:nitrogen+
                evar:nitrogen+(1|site/experiment/uniqueid)
              +(1|year), data = df_subset)
  
  # store the model summary (you can change this)
  results_list[[as.character(s)]] <- summary(mod)
}

# Extract fixed effects from the list of summaries
coef_df <- map_df(
  names(results_list),
  ~ {
    sm <- results_list[[.x]]
    fe <- as.data.frame(sm$coefficients)
    fe$term <- rownames(fe)
    fe$year_left_out <- .x
    fe
  }
)
#rename terms appropriately
coef_df<-coef_df%>%
  mutate(term=case_when(term=="dominant_relative_abund_zero"~"dominant",
                        term=="nitrogennutrients"~"nutrients",
                        term=="evar"~"evenness",
                        term %in% "nitrogennutrients" ~ "nutrients",
                        term %in% "dominant_relative_abund_zero:nitrogennutrients" ~ "dominant:nutrients",
                        term %in% "nitrogennutrients:evar" ~ "evenness:nutrients",
                        .default=term))
year_resil_wet_sens<-ggplot(coef_df, aes(x = year_left_out, y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = Estimate - `Std. Error`,
                    ymax = Estimate + `Std. Error`),
                width = 0.2) +
  facet_wrap(~ term, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(
    title = "sensitivity analysis of resilience to wet extreme events",
    x = "year left out",
    y = "ln(Resilience ± SE)"
  )

# Sensitivity analysis resilience dry (leave-one-year out) 
years <- unique(plot_ece_9_cn_prior_rm_dry$year)
results_list <- list()

for (s in years) {
  cat("Leaving out year:", s, "\n")
  # remove one site
  df_subset <- plot_ece_9_cn_prior_rm_dry %>% filter(year != s)
  # refit model
  mod <- lmer(log(resilience_n)~richness+dominant_relative_abund_zero+nitrogen+evar+
                dominant_relative_abund_zero:nitrogen+
                evar:nitrogen+(1|site/experiment/uniqueid)
              +(1|year), data = df_subset)

  # store the model summary (you can change this)
  results_list[[as.character(s)]] <- summary(mod)
}

# Extract fixed effects from the list of summaries
coef_df <- map_df(
  names(results_list),
  ~ {
    sm <- results_list[[.x]]
    fe <- as.data.frame(sm$coefficients)
    fe$term <- rownames(fe)
    fe$year_left_out <- .x
    fe
  }
)
#rename terms appropriately
coef_df<-coef_df%>%
  mutate(term=case_when(term=="dominant_relative_abund_zero"~"dominant",
                        term=="nitrogennutrients"~"nutrients",
                        term=="evar"~"evenness",
                        term %in% "nitrogennutrients" ~ "nutrients",
                        term %in% "dominant_relative_abund_zero:nitrogennutrients" ~ "dominant:nutrients",
                        term %in% "nitrogennutrients:evar" ~ "evenness:nutrients",
                        .default=term))
year_resil_dry_sens<-ggplot(coef_df, aes(x = year_left_out, y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = Estimate - `Std. Error`,
                    ymax = Estimate + `Std. Error`),
                width = 0.2) +
  facet_wrap(~ term, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(
    title = "sensitivity analysis of resilience to dry extreme events",
    x = "year left out",
    y = "ln(Resilience ± SE)"
  )
#####combine sensitivity figures for sites####
resist_dry_sens / resist_wet_sens+ plot_layout(guides = "collect")+resil_wet_sens& plot_annotation(tag_levels = 'A')&theme(legend.position = "bottom")
#&scale_fil

#####combine sensitivity fig for years
year_resist_dry_sens / year_resist_wet_sens+ plot_layout(guides = "collect")+year_resil_wet_sens& plot_annotation(tag_levels = 'A')&theme(legend.position = "bottom")
#&scale_fil

#calculating Mean annual temp and precipitation for each site####
precip_temp<-meta%>%
  group_by(site)%>%
  summarise(meantemp=(mean(meantemp, na.rm=T)),
            meanprecip=(mean(annualprecip, na.rm=T)))




