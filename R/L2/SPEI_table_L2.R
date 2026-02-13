# TITLE:        LTER Grassland Rock: SPEI exploration and LRR analyses 
# AUTHORS:      Ashley Darst
# COLLABORATORS:  
# DATA INPUT:   Data imported as csv files from shared Google drive L2 (plot_metrics_SPEI) folder
# DATA OUTPUT:  SPEI figures, analyses to choose SPEI duration, and LRR analyses
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
library(patchwork)
library(ggeffects)
library(ggsignif)
library(purrr)

# Set working directory 
L0_dir <- Sys.getenv("L0DIR")
L1_dir <- Sys.getenv("L1DIR")
L2_dir <- Sys.getenv("L2DIR")
list.files(L2_dir)

# Read in CSV files ----
plot_spei <- read.csv(file.path(L2_dir, "plot_metrics_SPEI_diversity_L2.csv"))
metadata <- read.csv(file.path(L1_dir, "metadata_L1.csv"))
# ece <- read.csv(file.path(L2_dir, "ece_resist_resil_spei9_L2.csv"))

# Make sure SPEI category is a factor
str(plot_spei)
plot_spei$spei12_category <- factor(plot_spei$spei12_category, levels = c("Extreme dry", "Moderate dry", "Normal", "Moderate wet", "Extreme wet"))
plot_spei$spei9_category <- factor(plot_spei$spei9_category, levels = c("Extreme dry", "Moderate dry", "Normal", "Moderate wet", "Extreme wet"))

# Summarize SPEI by year
spei_summary <- plot_spei %>%
  group_by(year, spei12_category, site) %>%
  summarize(spei12 = mean(spei12),
            spei9 = mean(spei9),
            spei6 = mean(spei6),
            spei3 = mean(spei3))
spei9_summary <- plot_spei %>%
  drop_na(spei9_category) %>%
  group_by(year, spei9_category, site) %>%
  summarize(spei9 = mean(spei9))

# Merge metadata with SPEI to get treatment column
metadata_filter <- metadata %>%
  drop_na(uniqueid)
plot_filter <- plot_spei %>%
  drop_na(uniqueid) %>%
  filter(higher_order_organization != "RaMPs_Block1",
         higher_order_organization != "RaMPs_Block2",
         higher_order_organization != "RaMPs_Block3")
plot_trt <- left_join(plot_filter, metadata_filter)


# Look at SPEI plots for each site in the range the data exist
cdr_spei <- spei9_summary %>%
  filter(site == "CDR") %>%
  ggplot(aes (x = year, y = spei9)) +
  geom_line(alpha = 0.5) + 
  geom_point( aes(color = spei9_category), size = 1.25) + theme_bw() +
  scale_x_continuous(breaks = seq(1980, 2025, by = 10)) +
  scale_y_continuous(limits = c(-2.2, 2.2)) +
  annotate("text", label = "CDR",
            x = 1986, y = 2, size = 5, colour = "black") + 
  scale_color_manual(values = c("#F5191CFF", "#E78200FF", "#EAC728FF", "#81BB95FF", "#3B99B1FF"))+
  geom_hline(yintercept = -1.28, col = "#F5191CFF", linetype = "dashed")+
  geom_hline(yintercept = 1.28, col = "#3B99B1FF", linetype = "dashed") +
  labs(x = "Year", y = "SPEI-9", color = "Event Type")

kbs_spei <- spei9_summary %>%
  filter(site == "KBS") %>%
  ggplot(aes (x = year, y = spei9)) +
  geom_line(alpha = 0.5) + 
  geom_point( aes(color = spei9_category), size = 1.25) + theme_bw() +
  scale_x_continuous(breaks = seq(1990, 2025, by = 10)) +
  scale_y_continuous(limits = c(-2.2, 2.2)) +
  annotate("text", label = "KBS",
            x = 1993, y = 2, size = 5, colour = "black") + 
  scale_color_manual(values = c("#F5191CFF", "#E78200FF", "#EAC728FF", "#81BB95FF", "#3B99B1FF"))+
  geom_hline(yintercept = -1.28, col = "#F5191CFF", linetype = "dashed")+
  geom_hline(yintercept = 1.28, col = "#3B99B1FF", linetype = "dashed") +
  labs(x = "Year", y = "SPEI-9", color = "Event Type")

knz_spei <- spei9_summary %>%
  filter(site == "KNZ") %>%
  ggplot(aes(x = year, y = spei9)) +
  geom_line(alpha = 0.5) + 
  geom_point( aes(color = spei9_category), size = 1.25) + theme_bw() +
  scale_x_continuous(breaks = seq(1990, 2025, by = 10)) +
  scale_y_continuous(limits = c(-2.2, 2.2)) +
  annotate("text", label = "KNZ",
            x = 1988, y = 2, size = 5, colour = "black") + 
  scale_color_manual(values = c("#F5191CFF", "#E78200FF", "#EAC728FF", "#81BB95FF", "#3B99B1FF"))+
  geom_hline(yintercept = -1.28, col = "#F5191CFF", linetype = "dashed")+
  geom_hline(yintercept = 1.28, col = "#3B99B1FF", linetype = "dashed") +
  labs(x = "Year", y = "SPEI-9", color = "Event Type")

# FIGURE S2 ----
cdr_spei + kbs_spei + knz_spei + plot_layout(guides = 'collect', axes = 'collect_y') & plot_annotation(tag_levels = 'A') 

# Choosing SPEI ----
# What duration of SPEI best explains plot biomass?

# Separate out control plots
plot_control <- plot_trt %>%
  filter(treatment == "control")

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
control_spei3.lm2.ln.ml <- lmer(log1p(plot_biomass) ~ spei3 + I(spei3^2) + (1|site/experiment/uniqueid) + (1|year), data = plot_control, REML = F)
summary(control_spei3.lm2.ln) # Quadratic term significant
control_spei3.lm.ln <- lmer(log1p(plot_biomass) ~ spei3 + (1|site/experiment/uniqueid) + (1|year), data = plot_control)
control_spei3.lm.ln.ml <- lmer(log1p(plot_biomass) ~ spei3 + (1|site/experiment/uniqueid) + (1|year), data = plot_control, REML = F)

simres <- simulateResiduals(control_spei3.lm2.ln)
plot(simres) # Better

control_spei6.lm2.ln <- lmer(log1p(plot_biomass) ~ spei6 + I(spei6^2) + (1|site/experiment/uniqueid)
                             + (1|year), data = plot_control)
control_spei6.lm2.ln.ml <- lmer(log1p(plot_biomass) ~ spei6 + I(spei6^2) + (1|site/experiment/uniqueid)
                                + (1|year), data = plot_control, REML = F)
summary(control_spei6.lm2.ln) # Quadratic term not significant
control_spei6.lm.ln <- lmer(log1p(plot_biomass) ~ spei6 + (1|site/experiment/uniqueid)
                            + (1|year), data = plot_control)
control_spei6.lm.ln.ml <- lmer(log1p(plot_biomass) ~ spei6 + (1|site/experiment/uniqueid)
                               + (1|year), data = plot_control, REML = F)
simres <- simulateResiduals(control_spei6.lm2.ln)
plot(simres) # Better

control_spei9.lm2.ln <- lmer(log1p(plot_biomass) ~ spei9 + I(spei9^2) + (1|site/experiment/uniqueid) + (1|year), data = plot_control)
control_spei9.lm2.ln.ml <- lmer(log1p(plot_biomass) ~ spei9 + I(spei9^2) + (1|site/experiment/uniqueid) + (1|year), data = plot_control, REML = F)
summary(control_spei9.lm2.ln) # Quadratic term marginally significant
control_spei9.lm.ln <- lmer(log1p(plot_biomass) ~ spei9 + (1|site/experiment/uniqueid) + (1|year), data = plot_control)
control_spei9.lm.ln.ml <- lmer(log1p(plot_biomass) ~ spei9 + (1|site/experiment/uniqueid) + (1|year), data = plot_control, REML = F)
summary(control_spei9.lm.ln)
simres <- simulateResiduals(control_spei9.lm2.ln)
plot(simres) # Better
simres <- simulateResiduals(control_spei9.lm.ln)
plot(simres) # Okay

control_spei12.lm2.ln <- lmer(log1p(plot_biomass) ~ spei12 + I(spei12^2) + (1|site/experiment/uniqueid) + (1|year), data = plot_control, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) # Failed to converge with defaults
control_spei12.lm2.ln.ml <- lmer(log1p(plot_biomass) ~ spei12 + I(spei12^2) + (1|site/experiment/uniqueid) + (1|year), data = plot_control, REML = F)
summary(control_spei12.lm2.ln)
control_spei12.lm.ln <- lmer(log1p(plot_biomass) ~ spei12 + (1|site/experiment/uniqueid) + (1|year), data = plot_control, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) # Failed to converge with defaults
control_spei12.lm.ln.ml <- lmer(log1p(plot_biomass) ~ spei12 + (1|site/experiment/uniqueid) + (1|year), data = plot_control, REML = F)
simres <- simulateResiduals(control_spei12.lm2.ln)
plot(simres) # Better

# Need to use ML to compare AIC
AICctab(control_spei3.lm2.ln.ml, control_spei3.lm.ln.ml, control_spei6.lm2.ln.ml, control_spei6.lm.ln.ml, control_spei9.lm2.ln.ml, control_spei9.lm.ln.ml, control_spei12.lm2.ln.ml, control_spei12.lm.ln.ml)

# Compare model R2 like Robinson (using REML)
r.squaredGLMM(control_spei3.lm2.ln)
r.squaredGLMM(control_spei3.lm.ln)
r.squaredGLMM(control_spei6.lm2.ln)
r.squaredGLMM(control_spei6.lm.ln)
r.squaredGLMM(control_spei9.lm2.ln) # best R2m
r.squaredGLMM(control_spei9.lm.ln)  # second best R2m
r.squaredGLMM(control_spei12.lm2.ln)
r.squaredGLMM(control_spei12.lm.ln)

# FIGURE S1 ----
# Plot the SPEI9 model
plot_model(
  control_spei9.lm2.ln,
  type = "pred",
  terms="spei9  [all]",
  show.data = TRUE,
  title = "",
  axis.title = c("SPEI-9", expression("Aboveground biomass (g/m"^2*")"))) +
  theme_bw() +
  scale_y_continuous(trans = "log1p", breaks = c(0, 10, 100, 1000, 10000))

# LRR ----
# Add experiment column to all plots
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

# Subset nitrogen fertilized plots and control plots
plot_n_sub <- plot_trt %>%
  filter(nutrients_added == "NPK+" | nutrients_added == "N" | nutrients_added == "NP" | nutrients_added == "NPK" | nutrients_added == "NK" | nutrients_added == "NPK+Fence" | nutrients_added == "none" | nutrients_added == "no_fertilizer")

# Add column with nitrogen yes/no
plot_n_sub <- plot_n_sub %>%
  mutate(nitrogen = dplyr::recode(nutrients_added, "NPK+" = "N", "NP" = "N", "NPK" = "N", "NK" = "N", "NPK+Fence" = "N", "none" = "no_fertilizer"))

# Filter out treatments
plot_n_sub <- plot_n_sub %>%
  filter(disturbance != "disturbed") %>%
  filter(treatment != "irrigated") %>%
  filter(treatment != "altered") %>%
  filter(grazing != "grazed") %>%
  filter(treatment != "insecticide") #remove plots with additional treatment

# Make year factor for random effect
plot_n_sub$year <- as.factor(plot_n_sub$year)

# LRR Nitrogen plots

# Subset nitrogen only
plot_n <- plot_trt %>%
  filter(nutrients_added == "NPK+" | nutrients_added == "N" | nutrients_added == "NP" | nutrients_added == "NPK" | nutrients_added == "NK" | nutrients_added == "NPK+Fence")

# Make histogram of nitrogen amounts (total observations)
plot_n %>%
  ggplot(aes(nitrogen_amount)) +
  geom_histogram(color = "#000000", fill = "darkgreen") +
  theme_bw()

# Make histogram of nitrogen amounts (plot level)
plot_n %>%
  distinct(uniqueid, .keep_all = TRUE) %>%
  ggplot(aes(nitrogen_amount)) +
  geom_histogram(color = "#000000", fill = "darkgreen") +
  theme_bw()


## ANPP ----
plot_sub_lrr9 <- plot_n_sub %>%
  group_by(uniqueid) %>%
  arrange(uniqueid, year) %>%
  mutate(prior_year_biomass = lag(plot_biomass),
         prior_year_type = lag(spei9_category)) %>%
  ungroup()
normal_biomass9sub_norm <- plot_sub_lrr9  %>%
  group_by(uniqueid) %>%
  filter(spei9_category == "Normal") %>% # ONLY normal year biomass
  summarize(average_normal_biomass = mean(plot_biomass, na.rm=T))
plot_sub_lrr9_norm <- left_join(plot_sub_lrr9, normal_biomass9sub_norm)

lrr9sub_norm <- plot_sub_lrr9_norm %>%
  # filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
  filter(spei9_category == "Extreme wet" | spei9_category == "Extreme dry") %>%
  mutate(LRR = log(plot_biomass / average_normal_biomass))  # Calculate the log response ratio

# Take absolute value of SPEI9
lrr9sub_norm <- lrr9sub_norm %>%
  mutate(abs_spei9 = abs(spei9))

# Simple model
lrr.lm9sub2_norm <- lmer(LRR ~ spei9_category*nitrogen + (1|site/experiment/uniqueid) + (1|year), data = lrr9sub_norm)
summary(lrr.lm9sub2_norm)
emm_biomass <- emmeans(lrr.lm9sub2_norm, pairwise ~ spei9_category * nitrogen, infer = T)
summary(emm_biomass)
simres <- simulateResiduals(lrr.lm9sub2_norm)
plot(simres)

# Sensitivity analysis (leave-one-site out) # created with ChatGPT
sites <- unique(lrr9sub_norm$site)
results_list <- list()

for (s in sites) {
  cat("Leaving out site:", s, "\n")
  # remove one site
  df_subset <- lrr9sub_norm %>% filter(site != s)
  # refit model
  mod <- lmer(
    LRR ~ spei9_category * nitrogen +
      (1 | site/experiment/uniqueid) +
      (1 | year),
    data = df_subset
  )
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
coef_df <- coef_df %>%
  mutate(term = case_when(
    term == "(Intercept)" ~ "(Intercept)",
    term == "spei9_categoryExtreme wet" ~ "event type",
    term == "nitrogenno_fertilizer" ~ "nutrients",
    term == "spei9_categoryExtreme wet:nitrogenno_fertilizer" ~ "event type:nutrients"
  )) %>%
  mutate(across(term, ~factor(., levels=c("(Intercept)","event type","nutrients", "event type:nutrients"))))

fig_sens_anpp_site <- ggplot(coef_df, aes(x = site_left_out, y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = Estimate - `Std. Error`,
                    ymax = Estimate + `Std. Error`),
                width = 0.2) +
  facet_wrap(~ term, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(
    title = "Leave-one-site-out biomass",
    x = "Site left out",
    y = "Fixed effect estimate ± SE"
  )

# Leave one year out
years <- unique(lrr9sub_norm$year)
results_list <- list()

for (s in years) {
  cat("Leaving out year:", s, "\n")
  # remove one year
  df_subset <- lrr9sub_norm %>% filter(year != s)
  # refit model
  mod <- lmer(
    LRR ~ spei9_category * nitrogen +
      (1 | site/experiment/uniqueid) +
      (1 | year),
    data = df_subset
  )
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
coef_df <- coef_df %>%
  mutate(term = case_when(
    term == "(Intercept)" ~ "(Intercept)",
    term == "spei9_categoryExtreme wet" ~ "event type",
    term == "nitrogenno_fertilizer" ~ "nutrients",
    term == "spei9_categoryExtreme wet:nitrogenno_fertilizer" ~ "event type:nutrients"
  )) %>%
  mutate(across(term, ~factor(., levels=c("(Intercept)","event type","nutrients", "event type:nutrients"))))

fig_sens_anpp_year <- ggplot(coef_df, aes(x = site_left_out, y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = Estimate - `Std. Error`,
                    ymax = Estimate + `Std. Error`),
                width = 0.2) +
  facet_wrap(~ term, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(
    title = "Leave-one-year-out biomass",
    x = "Year left out",
    y = "Fixed effect estimate ± SE"
  )

plot_model(
  lrr.lm9sub2_norm,
  type = "pred",
  terms= c("spei9_category", "nitrogen"),
  show.data = F,
  title = "") + 
  theme_bw() +
  geom_hline(yintercept=0, linetype='dashed', col = 'black') +
  labs(x = "Event type", y = "ANPP LRR") +
  theme(legend.position="none")

# LRR plot for ANPP
cols <- c("Extreme dry" = "#E41A1C", "Extreme wet" = "#377EB8")
lrr9sub_norm %>%
  ggplot(aes(x = spei9_category, y = LRR, col = spei9_category)) +
  stat_summary(fun.data = mean_cl_boot, position = position_dodge(0.2), aes(shape = nitrogen)) + 
  geom_hline(yintercept=0, linetype='dashed', col = 'black')+
  theme_bw() +
  scale_color_manual(values = cols)
anpp_plot_norm <- lrr9sub_norm %>%
  mutate(nitrogen = relevel(as.factor(nitrogen), 'no_fertilizer', 'nitrogen')) %>%
  ggplot(aes(x = spei9_category, y = LRR, col = spei9_category)) +
  stat_summary(fun.data = mean_cl_boot, position = position_dodge(0.2), aes(shape = nitrogen)) + 
  geom_hline(yintercept=0, linetype='dashed', col = 'black')+
  theme_bw() +
  labs(x = "Event type", y = "ANPP LRR") +
  theme(legend.position="none") +
  scale_color_manual(values = cols)

# Using ggpredict
lrr.lm9sub2_norm.df <- predict_response(lrr.lm9sub2_norm, terms = c("spei9_category", "nitrogen"))

# Add significance column
sig <- c("TN" = 17, "FN" = 2, "TC" = 16, "FC" = 1)
sig_anpp <- c("TN", "FC", "TN", "FC")
lrr.lm9sub2_norm.df$sig <- sig_anpp

anpp_plot_norm_pred <- lrr.lm9sub2_norm.df %>%
  mutate(group = relevel(as.factor(group), 'no_fertilizer', 'nitrogen')) %>%
  drop_na() %>%
  ggplot(aes(x, predicted, color = x)) +
  geom_point(aes(shape = sig), position = position_dodge(0.2), size = 3) +
  geom_errorbar(aes(ymin = (predicted - std.error), ymax = (predicted + std.error), group = group), width = 0.2, position = position_dodge(0.2)) +
  geom_hline(yintercept=0, linetype='dashed', col = 'black')+
  theme_bw() +
  labs(x = "Event type", y = "Aboveground biomass LRR") +
  scale_color_manual(values = cols) +
  scale_shape_manual(values = sig) +
  theme(legend.position="none") +
  geom_signif(map_signif_level = F, tip_length = 0, y_position = 0.25, xmin = 0.95, xmax = 1.95, annotations = "***", color = "black") +
  geom_signif(map_signif_level = F, tip_length = 0, y_position = 0.45, xmin = 1.05, xmax = 2.05, annotations = "***", color = "black") +
  annotate(geom = "text", x=1, y = 0.05, label="***", color="black") +
  annotate(geom = "text", x=2, y = 0.4, label="**", color="black") +
  coord_cartesian(ylim = c(NA, 0.48))

# Three way interaction (for supplement)
lrr9sub_norm$nitrogen <- fct_recode(lrr9sub_norm$nitrogen, "nutrients" = "N", "control" = "no_fertilizer")
lrr.lm9sub_norm3 <- lmer(LRR ~ abs_spei9*spei9_category*nitrogen + (1|site/experiment/uniqueid) + (1|year), data = lrr9sub_norm)
summary(lrr.lm9sub_norm3)
emtrends(lrr.lm9sub_norm3, pairwise ~ spei9_category * nitrogen, var = "abs_spei9", infer = T)
simres <- simulateResiduals(lrr.lm9sub_norm3)
plot(simres)

anpp_plot_norm3 <- plot_model(
  lrr.lm9sub_norm3,
  type = "pred",
  terms= c("abs_spei9", "spei9_category", "nitrogen"),
  show.data = T,
  title = "",
  axis.title = c("|SPEI-9|", "Aboveground biomass LRR"),
  legend.title = "Event type") +
  theme_bw() +
  geom_hline(yintercept=0, linetype='dashed', col = 'black')


## Richness ----
plot_sub_lrr9_rich <- plot_n_sub %>%
  group_by(uniqueid) %>%
  arrange(uniqueid, year) %>%
  mutate(prior_year_richness = lag(Richness),
         prior_year_type = lag(spei9_category)) %>%
  ungroup()
normal_biomass9sub_rich_norm <- plot_sub_lrr9_rich  %>%
  group_by(uniqueid) %>%
  filter(spei9_category == "Normal") %>%
  summarize(average_normal_rich = mean(Richness, na.rm=T))
plot_sub_lrr9_rich_norm <- left_join(plot_sub_lrr9_rich, normal_biomass9sub_rich_norm)

lrr9sub_rich_norm <- plot_sub_lrr9_rich_norm %>%
  # filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
  filter(spei9_category == "Extreme wet" | spei9_category == "Extreme dry") %>%
  mutate(LRR = log(Richness / average_normal_rich))  # Calculate the log response ratio

# Take absolute value of SPEI9
lrr9sub_rich_norm <- lrr9sub_rich_norm %>%
  mutate(abs_spei9 = abs(spei9))

lrr.lm9sub_rich_norm <- lmer(LRR ~ spei9_category*nitrogen + (1|site/experiment/uniqueid) + (1|year), data = lrr9sub_rich_norm)
summary(lrr.lm9sub_rich_norm)
emm_rich <- emmeans(lrr.lm9sub_rich_norm, pairwise ~ spei9_category * nitrogen, infer = T)
summary(emm_rich)
simres <- simulateResiduals(lrr.lm9sub_rich_norm)
plot(simres)

# Sensitivity analysis (leave-one-site out) # created with ChatGPT
sites <- unique(lrr9sub_rich_norm$site)
results_list <- list()

for (s in sites) {
  cat("Leaving out site:", s, "\n")
  # remove one site
  df_subset <- lrr9sub_rich_norm %>% filter(site != s)
  # refit model
  mod <- lmer(
    LRR ~ spei9_category * nitrogen +
      (1 | site/experiment/uniqueid) +
      (1 | year),
    data = df_subset
  )
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
coef_df <- coef_df %>%
  mutate(term = case_when(
    term == "(Intercept)" ~ "(Intercept)",
    term == "spei9_categoryExtreme wet" ~ "event type",
    term == "nitrogenno_fertilizer" ~ "nutrients",
    term == "spei9_categoryExtreme wet:nitrogenno_fertilizer" ~ "event type:nutrients"
  )) %>%
  mutate(across(term, ~factor(., levels=c("(Intercept)","event type","nutrients", "event type:nutrients"))))

fig_sens_rich_site <- ggplot(coef_df, aes(x = site_left_out, y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = Estimate - `Std. Error`,
                    ymax = Estimate + `Std. Error`),
                width = 0.2) +
  facet_wrap(~ term, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(
    title = "Leave-one-site-out richness",
    x = "Site left out",
    y = "Fixed effect estimate ± SE"
  )

# Leave one year out
years <- unique(lrr9sub_rich_norm$year)
results_list <- list()

for (s in years) {
  cat("Leaving out year:", s, "\n")
  # remove one year
  df_subset <- lrr9sub_rich_norm %>% filter(year != s)
  # refit model
  mod <- lmer(
    LRR ~ spei9_category * nitrogen +
      (1 | site/experiment/uniqueid) +
      (1 | year),
    data = df_subset
  )
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
coef_df <- coef_df %>%
  mutate(term = case_when(
    term == "(Intercept)" ~ "(Intercept)",
    term == "spei9_categoryExtreme wet" ~ "event type",
    term == "nitrogenno_fertilizer" ~ "nutrients",
    term == "spei9_categoryExtreme wet:nitrogenno_fertilizer" ~ "event type:nutrients"
  )) %>%
  mutate(across(term, ~factor(., levels=c("(Intercept)","event type","nutrients", "event type:nutrients"))))

fig_sens_rich_year <- ggplot(coef_df, aes(x = site_left_out, y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = Estimate - `Std. Error`,
                    ymax = Estimate + `Std. Error`),
                width = 0.2) +
  facet_wrap(~ term, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(
    title = "Leave-one-year-out richness",
    x = "Year left out",
    y = "Fixed effect estimate ± SE"
  )

# LRR plot for richness 
lrr9sub_rich_norm %>%
  ggplot(aes(x = spei9_category, y = LRR, col = nitrogen)) +
  stat_summary(fun.data = mean_cl_boot, position = position_dodge(0.2)) + 
  geom_hline(yintercept=0, linetype='dashed', col = 'black')+
  theme_bw()
rich_plot_norm <- lrr9sub_rich_norm %>%
  mutate(nitrogen = relevel(as.factor(nitrogen), 'no_fertilizer', 'nitrogen')) %>%
  ggplot(aes(x = spei9_category, y = LRR, col = spei9_category)) +
  stat_summary(fun.data = mean_cl_boot, position = position_dodge(0.2), aes(shape = nitrogen)) + 
  geom_hline(yintercept=0, linetype='dashed', col = 'black')+
  theme_bw() +
  labs(x = "Event type", y = "Richness LRR") +
  theme(legend.position="none") +
  scale_color_manual(values = cols)

plot_model(
  lrr.lm9sub_rich_norm,
  type = "pred",
  terms= c("spei9_category", "nitrogen"),
  show.data = F
) + theme_bw()

# Using ggpredict
lrr9sub_rich_norm.df <- predict_response(lrr.lm9sub_rich_norm, terms = c("spei9_category", "nitrogen"))

#Add significance column
sig_rich <- c("TN", "FC", "FN", "FC")
lrr9sub_rich_norm.df$sig <- sig_rich

rich_plot_norm_pred <- lrr9sub_rich_norm.df %>%
  mutate(group = relevel(as.factor(group), 'no_fertilizer', 'nitrogen')) %>%
  drop_na() %>%
  ggplot(aes(x, predicted, color = x)) +
  geom_point(aes(shape = sig), position = position_dodge(0.2), size = 3) +
  geom_errorbar(aes(ymin = (predicted - std.error), ymax = (predicted + std.error), group = group), width = 0.2, position = position_dodge(0.2)) +
  geom_hline(yintercept=0, linetype='dashed', col = 'black')+
  theme_bw() +
  labs(x = "Event type", y = "Species richness LRR") +
  scale_color_manual(values = cols) +
  scale_shape_manual(values = sig) +
  theme(legend.position="none") +
  annotate(geom = "text", x=1, y = 0.1, label= "***", color="black") +
  annotate(geom = "text", x=2, y = 0.05, label= "+", color="black")

# Three way interaction (for supplement)
lrr9sub_rich_norm$nitrogen <- fct_recode(lrr9sub_rich_norm$nitrogen, "nutrients" = "N", "control" = "no_fertilizer")
lrr.lm9sub_rich_norm3 <- lmer(LRR ~ abs_spei9*spei9_category*nitrogen + (1|site/experiment/uniqueid) + (1|year), data = lrr9sub_rich_norm)
summary(lrr.lm9sub_rich_norm3)
emtrends(lrr.lm9sub_rich_norm3, pairwise ~ spei9_category * nitrogen, var = "abs_spei9", infer = T)
simres <- simulateResiduals(lrr.lm9sub_rich_norm3)
plot(simres)

rich_plot_norm3 <- plot_model(
  lrr.lm9sub_rich_norm3,
  type = "pred",
  terms= c("abs_spei9", "spei9_category", "nitrogen"),
  show.data = T,
  title = "",
  axis.title = c("|SPEI-9|", "Species richness LRR"),
  legend.title = "Event type") + 
  theme_bw() +
  geom_hline(yintercept=0, linetype='dashed', col = 'black')


## Dominance ----
# Zeros were not used in these analyses, similar to when treated as NA
plot_sub_lrr9_dom <- plot_n_sub %>%
  group_by(uniqueid) %>%
  arrange(uniqueid, year) %>%
  mutate(prior_year_dom = lag(dominant_relative_abund),
         prior_year_type = lag(spei9_category)) %>%
  ungroup()
normal_biomass9sub_dom_norm <- plot_sub_lrr9_dom  %>%
  group_by(uniqueid) %>%
  filter(spei9_category == "Normal") %>%
  summarize(average_normal_dom = mean(dominant_relative_abund, na.rm=T))
plot_sub_lrr9_dom_norm <- left_join(plot_sub_lrr9_dom, normal_biomass9sub_dom_norm)

lrr9sub_dom_norm <- plot_sub_lrr9_dom_norm %>%
  # filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
  filter(spei9_category == "Extreme wet" | spei9_category == "Extreme dry") %>%
  mutate(LRR = log(dominant_relative_abund / average_normal_dom))  # Calculate the log response ratio

# Take absolute value of SPEI9
lrr9sub_dom_norm <- lrr9sub_dom_norm %>%
  mutate(abs_spei9 = abs(spei9))

lrr.lm9sub_dom_norm <- lmer(LRR ~ spei9_category*nitrogen + (1|site/experiment/uniqueid) + (1|year), data = lrr9sub_dom_norm)
summary(lrr.lm9sub_dom_norm)
emm_dom <- emmeans(lrr.lm9sub_dom_norm, pairwise ~ spei9_category * nitrogen, infer = T)
summary(emm_dom)
simres <- simulateResiduals(lrr.lm9sub_dom_norm)
plot(simres)

# Sensitivity analysis (leave-one-site out) # created with ChatGPT
sites <- unique(lrr9sub_dom_norm$site)
results_list <- list()

for (s in sites) {
  cat("Leaving out site:", s, "\n")
  # remove one site
  df_subset <- lrr9sub_dom_norm %>% filter(site != s)
  # refit model
  mod <- lmer(
    LRR ~ spei9_category * nitrogen +
      (1 | site/experiment/uniqueid) +
      (1 | year),
    data = df_subset
  )
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
coef_df <- coef_df %>%
  mutate(term = case_when(
    term == "(Intercept)" ~ "(Intercept)",
    term == "spei9_categoryExtreme wet" ~ "event type",
    term == "nitrogenno_fertilizer" ~ "nutrients",
    term == "spei9_categoryExtreme wet:nitrogenno_fertilizer" ~ "event type:nutrients"
  )) %>%
  mutate(across(term, ~factor(., levels=c("(Intercept)","event type","nutrients", "event type:nutrients"))))

fig_sens_dom_site <- ggplot(coef_df, aes(x = site_left_out, y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = Estimate - `Std. Error`,
                    ymax = Estimate + `Std. Error`),
                width = 0.2) +
  facet_wrap(~ term, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(
    title = "Leave-one-site-out dominance",
    x = "Site left out",
    y = "Fixed effect estimate ± SE"
  )

# Leave one year out
years <- unique(lrr9sub_dom_norm$year)
results_list <- list()

for (s in years) {
  cat("Leaving out year:", s, "\n")
  # remove one year
  df_subset <- lrr9sub_dom_norm %>% filter(year != s)
  # refit model
  mod <- lmer(
    LRR ~ spei9_category * nitrogen +
      (1 | site/experiment/uniqueid) +
      (1 | year),
    data = df_subset
  )
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
coef_df <- coef_df %>%
  mutate(term = case_when(
    term == "(Intercept)" ~ "(Intercept)",
    term == "spei9_categoryExtreme wet" ~ "event type",
    term == "nitrogenno_fertilizer" ~ "nutrients",
    term == "spei9_categoryExtreme wet:nitrogenno_fertilizer" ~ "event type:nutrients"
  )) %>%
  mutate(across(term, ~factor(., levels=c("(Intercept)","event type","nutrients", "event type:nutrients"))))

fig_sens_dom_year <- ggplot(coef_df, aes(x = site_left_out, y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = Estimate - `Std. Error`,
                    ymax = Estimate + `Std. Error`),
                width = 0.2) +
  facet_wrap(~ term, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(
    title = "Leave-one-year-out dominance",
    x = "Year left out",
    y = "Fixed effect estimate ± SE"
  )

# LRR plot for dominance normal only
lrr9sub_dom_norm %>%
  ggplot(aes(x = spei9_category, y = LRR, col = nitrogen)) +
  stat_summary(fun.data = mean_cl_boot, position = position_dodge(0.2)) + 
  geom_hline(yintercept=0, linetype='dashed', col = 'black')+
  theme_bw()
dom_plot_norm <- lrr9sub_dom_norm %>%
  mutate(nitrogen = relevel(as.factor(nitrogen), 'no_fertilizer', 'nitrogen')) %>%
  ggplot(aes(x = spei9_category, y = LRR, col = spei9_category)) +
  stat_summary(fun.data = mean_cl_boot, position = position_dodge(0.2), aes(shape = nitrogen)) + 
  geom_hline(yintercept=0, linetype='dashed', col = 'black')+
  theme_bw() +
  labs(x = "Event type", y = "Dominance LRR") +
  theme(legend.position="none") +
  scale_color_manual(values = cols)

plot_model(
  lrr.lm9sub_dom_norm,
  type = "pred",
  terms= c("spei9_category", "nitrogen"),
  show.data = F
) + theme_bw()

# Using ggpredict
lrr.lm9sub_dom_norm.df <- predict_response(lrr.lm9sub_dom_norm, terms = c("spei9_category", "nitrogen"))

#Add significance column
sig_dom <- c("TN", "FC", "FN", "TC")
lrr.lm9sub_dom_norm.df$sig <- sig_dom

dom_plot_norm_pred <- lrr.lm9sub_dom_norm.df %>%
  mutate(group = relevel(as.factor(group), 'no_fertilizer', 'nitrogen')) %>%
  mutate(sig = relevel(as.factor(sig), 'TC', 'TN', 'FN', 'FC',)) %>%
  drop_na() %>%
  ggplot(aes(x, predicted, color = x)) +
  geom_point(aes(shape = sig), position = position_dodge(0.2), size = 3) +
  geom_errorbar(aes(ymin = (predicted - std.error), ymax = (predicted + std.error), group = group), width = 0.2, position = position_dodge(0.2)) +
  geom_hline(yintercept=0, linetype='dashed', col = 'black')+
  theme_bw() +
  labs(x = "Event type", y = "Dominance LRR") +
  scale_color_manual(values = cols) +
  scale_shape_manual(values = sig) +
  theme(legend.position="none")

# Three way interaction (for supplement)
lrr9sub_dom_norm$nitrogen <- fct_recode(lrr9sub_dom_norm$nitrogen, "nutrients" = "N", "control" = "no_fertilizer")
lrr.lm9sub_dom_norm3 <- lmer(LRR ~ abs_spei9*spei9_category*nitrogen + (1|site/experiment/uniqueid) + (1|year), data = lrr9sub_dom_norm)
summary(lrr.lm9sub_dom_norm3)
emtrends(lrr.lm9sub_dom_norm3, pairwise ~ spei9_category * nitrogen, var = "abs_spei9", infer = T)
simres <- simulateResiduals(lrr.lm9sub_dom_norm3)
plot(simres)

dom_plot_norm3 <- plot_model(
  lrr.lm9sub_dom_norm3,
  type = "pred",
  terms= c("abs_spei9", "spei9_category", "nitrogen"),
  show.data = T,
  title = "",
  axis.title = c("|SPEI-9|", "Dominance LRR"),
  legend.title = "Event type") +
  theme_bw() +
  geom_hline(yintercept=0, linetype='dashed', col = 'black')


## Evenness ----
plot_sub_lrr9_ev <- plot_n_sub %>%
  group_by(uniqueid) %>%
  arrange(uniqueid, year) %>%
  mutate(prior_year_ev = lag(Evar),
         prior_year_type = lag(spei9_category)) %>%
  ungroup()
normal_biomass9sub_ev_norm <- plot_sub_lrr9_ev  %>%
  group_by(uniqueid) %>%
  filter(spei9_category == "Normal") %>%
  summarize(average_normal_ev = mean(Evar, na.rm=T))
plot_sub_lrr9_ev_norm <- left_join(plot_sub_lrr9_ev, normal_biomass9sub_ev_norm)

lrr9sub_ev_norm <- plot_sub_lrr9_ev_norm %>%
  # filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
  filter(spei9_category == "Extreme wet" | spei9_category == "Extreme dry") %>%
  mutate(LRR = log(Evar / average_normal_ev))  # Calculate the log response ratio

# Take absolute value of SPEI9
lrr9sub_ev_norm <- lrr9sub_ev_norm %>%
  mutate(abs_spei9 = abs(spei9))

lrr.lm9sub_ev_norm <- lmer(LRR ~ spei9_category*nitrogen + (1|site/experiment/uniqueid) + (1|year), data = lrr9sub_ev_norm)
summary(lrr.lm9sub_ev_norm)
emm_ev <- emmeans(lrr.lm9sub_ev_norm, pairwise ~ spei9_category * nitrogen, infer = T)
summary(emm_ev)
simres <- simulateResiduals(lrr.lm9sub_ev_norm)
plot(simres)

# Sensitivity analysis (leave-one-site out) # created with ChatGPT
sites <- unique(lrr9sub_ev_norm$site)
results_list <- list()

for (s in sites) {
  cat("Leaving out site:", s, "\n")
  # remove one site
  df_subset <- lrr9sub_ev_norm %>% filter(site != s)
  # refit model
  mod <- lmer(
    LRR ~ spei9_category * nitrogen +
      (1 | site/experiment/uniqueid) +
      (1 | year),
    data = df_subset
  )
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
coef_df <- coef_df %>%
  mutate(term = case_when(
    term == "(Intercept)" ~ "(Intercept)",
    term == "spei9_categoryExtreme wet" ~ "event type",
    term == "nitrogenno_fertilizer" ~ "nutrients",
    term == "spei9_categoryExtreme wet:nitrogenno_fertilizer" ~ "event type:nutrients"
  )) %>%
  mutate(across(term, ~factor(., levels=c("(Intercept)","event type","nutrients", "event type:nutrients"))))

fig_sens_ev_site <- ggplot(coef_df, aes(x = site_left_out, y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = Estimate - `Std. Error`,
                    ymax = Estimate + `Std. Error`),
                width = 0.2) +
  facet_wrap(~ term, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(
    title = "Leave-one-site-out evenness",
    x = "Site left out",
    y = "Fixed effect estimate ± SE"
  )

# Leave one year out
years <- unique(lrr9sub_ev_norm$year)
results_list <- list()

for (s in years) {
  cat("Leaving out year:", s, "\n")
  # remove one year
  df_subset <- lrr9sub_ev_norm %>% filter(year != s)
  # refit model
  mod <- lmer(
    LRR ~ spei9_category * nitrogen +
      (1 | site/experiment/uniqueid) +
      (1 | year),
    data = df_subset
  )
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
coef_df <- coef_df %>%
  mutate(term = case_when(
    term == "(Intercept)" ~ "(Intercept)",
    term == "spei9_categoryExtreme wet" ~ "event type",
    term == "nitrogenno_fertilizer" ~ "nutrients",
    term == "spei9_categoryExtreme wet:nitrogenno_fertilizer" ~ "event type:nutrients"
  )) %>%
  mutate(across(term, ~factor(., levels=c("(Intercept)","event type","nutrients", "event type:nutrients"))))

fig_sens_ev_year <- ggplot(coef_df, aes(x = site_left_out, y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = Estimate - `Std. Error`,
                    ymax = Estimate + `Std. Error`),
                width = 0.2) +
  facet_wrap(~ term, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(
    title = "Leave-one-year-out evenness",
    x = "Year left out",
    y = "Fixed effect estimate ± SE"
  )

# LRR plot for evenness normal only
lrr9sub_ev_norm %>%
  ggplot(aes(x = spei9_category, y = LRR, col = nitrogen)) +
  stat_summary(fun.data = mean_cl_boot, position = position_dodge(0.2)) + 
  geom_hline(yintercept=0, linetype='dashed', col = 'black')+
  theme_bw()
ev_plot_norm <- lrr9sub_ev_norm %>%
  mutate(nitrogen = relevel(as.factor(nitrogen), 'no_fertilizer', 'nitrogen')) %>%
  ggplot(aes(x = spei9_category, y = LRR, col = spei9_category)) +
  stat_summary(fun.data = mean_cl_boot, position = position_dodge(0.2), aes(shape = nitrogen)) + 
  geom_hline(yintercept=0, linetype='dashed', col = 'black')+
  theme_bw() +
  labs(x = "Event type", y = "Evenness LRR", shape = "Fertilizer") +
  scale_color_manual(values = cols, guide = "none") +
  scale_shape_manual(values = c(19, 17), labels = c('no nutrients', 'nutrients')) +
  theme(legend.title=element_blank())

plot_model(
  lrr.lm9sub_ev_norm,
  type = "pred",
  terms= c("spei9_category", "nitrogen"),
  show.data = F
) + theme_bw()

# Using ggpredict
lrr.lm9sub_ev_norm.df <- predict_response(lrr.lm9sub_ev_norm, terms = c("spei9_category", "nitrogen"))

#Add significance column
sig_ev <- c("TN", "FC", "FN", "TC")
lrr.lm9sub_ev_norm.df$sig <- sig_ev

ev_plot_norm_pred <- lrr.lm9sub_ev_norm.df %>%
  mutate(group = relevel(as.factor(group), 'no_fertilizer', 'nitrogen')) %>%
  drop_na() %>%
  ggplot(aes(x, predicted, color = x, shape = group)) +
  geom_point(position = position_dodge(0.2), size = 3) +
  geom_errorbar(aes(ymin = (predicted - std.error), ymax = (predicted + std.error), group = group), width = 0.2, position = position_dodge(0.2)) +
  geom_hline(yintercept=0, linetype='dashed', col = 'black')+
  theme_bw() +
  labs(x = "Event type", y = "Evenness LRR") +
  scale_color_manual(values = cols, guide = "none") +
  scale_shape_manual(values = c(1, 2), labels = c('control', 'nutrients')) +
  theme(legend.title=element_blank()) +
  annotate(geom = "text", x=2, y=0.13, label= "**", color="black") +
  geom_signif(map_signif_level = F, tip_length = 0, y_position = 0.44, xmin = 0.95, xmax = 1.95, annotations = "*", color = "black") +
  geom_signif(map_signif_level = F, tip_length = 0, y_position = 0.5, xmin = 1.05, xmax = 2.05, annotations = "***", color = "black") +
  coord_cartesian(ylim = c(NA, 0.52))

# Three way interaction (for supplement)
lrr9sub_ev_norm$nitrogen <- fct_recode(lrr9sub_ev_norm$nitrogen, "nutrients" = "N", "control" = "no_fertilizer")
lrr.lm9sub_ev_norm3 <- lmer(LRR ~ abs_spei9*spei9_category*nitrogen + (1|site/experiment/uniqueid) + (1|year), data = lrr9sub_ev_norm)
summary(lrr.lm9sub_ev_norm3)
emtrends(lrr.lm9sub_ev_norm3, pairwise ~ spei9_category * nitrogen, var = "abs_spei9", infer = T)
simres <- simulateResiduals(lrr.lm9sub_ev_norm3)
plot(simres)

ev_plot_norm3 <- plot_model(
  lrr.lm9sub_ev_norm3,
  type = "pred",
  terms= c("abs_spei9", "spei9_category", "nitrogen"),
  show.data = T,
  title = "",
  axis.title = c("|SPEI-9|", "Evenness LRR"),
  legend.title = "Event type") +
  theme_bw() +
  geom_hline(yintercept=0, linetype='dashed', col = 'black')


## Final LRR plot NORMAL
anpp_plot_norm + rich_plot_norm + dom_plot_norm + ev_plot_norm + plot_layout(guides = 'collect') & plot_annotation(tag_levels = 'A')

## FIGURE 3 ----
anpp_plot_norm_pred + rich_plot_norm_pred + dom_plot_norm_pred + ev_plot_norm_pred + plot_layout(guides = 'collect') & plot_annotation(tag_levels = 'A')

## FIGURE S3 ----
anpp_plot_norm3 + rich_plot_norm3 + dom_plot_norm3 + ev_plot_norm3 + plot_layout(guides = 'collect') & plot_annotation(tag_levels = 'A')

## FIGURE Supp ----
fig_sens_anpp_site + fig_sens_rich_site + fig_sens_dom_site + fig_sens_ev_site & plot_annotation(tag_levels = 'A')

fig_sens_anpp_year + fig_sens_rich_year + fig_sens_dom_year + fig_sens_ev_year & plot_annotation(tag_levels = 'A')

# ARCHIVED CODE ----
# # Calculate the number of extreme wet years
# CDR_ew <- length(which(spei_summary$site == "CDR" & spei_summary$spei12_category == "Extreme wet"))
# KBS_ew <- length(which(spei_summary$site == "KBS" & spei_summary$spei12_category == "Extreme wet"))
# KNZ_ew <- length(which(spei_summary$site == "KNZ" & spei_summary$spei12_category == "Extreme wet"))
# extreme_wet <- c(CDR_ew, KBS_ew, KNZ_ew)
# sites <- c("CDR", "KBS", "KNZ")
# 
# # Calculate the number of moderate wet years
# CDR_mw <- length(which(spei_summary$site == "CDR" & spei_summary$spei12_category == "Moderate wet"))
# KBS_mw <- length(which(spei_summary$site == "KBS" & spei_summary$spei12_category == "Moderate wet"))
# KNZ_mw <- length(which(spei_summary$site == "KNZ" & spei_summary$spei12_category == "Moderate wet"))
# moderate_wet <- c(CDR_mw, KBS_mw, KNZ_mw)
# 
# # Calculate the number of extreme dry years
# CDR_ed <- length(which(spei_summary$site == "CDR" & spei_summary$spei12_category == "Extreme dry"))
# KBS_ed <- length(which(spei_summary$site == "KBS" & spei_summary$spei12_category == "Extreme dry"))
# KNZ_ed <- length(which(spei_summary$site == "KNZ" & spei_summary$spei12_category == "Extreme dry"))
# extreme_dry <- c(CDR_ed, KBS_ed, KNZ_ed)
# 
# # Calculate the number of moderate dry years
# CDR_md <- length(which(spei_summary$site == "CDR" & spei_summary$spei12_category == "Moderate dry"))
# KBS_md <- length(which(spei_summary$site == "KBS" & spei_summary$spei12_category == "Moderate dry"))
# KNZ_md <- length(which(spei_summary$site == "KNZ" & spei_summary$spei12_category == "Moderate dry"))
# moderate_dry <- c(CDR_md, KBS_md, KNZ_md)
# 
# # Calculate the number of extreme wet then extreme dry years
# spei_summary_CDR <- spei_summary %>%
#   filter(site == "CDR")
# CDR_ewd <- sum(spei_summary_CDR$spei12_category[-nrow(spei_summary_CDR)] == "Extreme wet" & spei_summary_CDR$spei12_category[-1] == "Extreme dry") # This was written by ChatGPT
# spei_summary_KBS <- spei_summary %>%
#   filter(site == "KBS")
# KBS_ewd <- sum(spei_summary_KBS$spei12_category[-nrow(spei_summary_KBS)] == "Extreme wet" & spei_summary_KBS$spei12_category[-1] == "Extreme dry")
# spei_summary_KNZ <- spei_summary %>%
#   filter(site == "KNZ")
# KNZ_ewd <- sum(spei_summary_KNZ$spei12_category[-nrow(spei_summary_KNZ)] == "Extreme wet" & spei_summary_KNZ$spei12_category[-1] == "Extreme dry")
# 
# extreme_wet_dry <- c(CDR_ewd, KBS_ewd, KNZ_ewd)
# 
# # Calculate the number of extreme dry then extreme wet years
# CDR_edw <- sum(spei_summary_CDR$spei12_category[-nrow(spei_summary_CDR)] == "Extreme dry" & spei_summary_CDR$spei12_category[-1] == "Extreme wet")
# KBS_edw <- sum(spei_summary_KBS$spei12_category[-nrow(spei_summary_KBS)] == "Extreme dry" & spei_summary_KBS$spei12_category[-1] == "Extreme wet")
# KNZ_edw <- sum(spei_summary_KNZ$spei12_category[-nrow(spei_summary_KNZ)] == "Extreme dry" & spei_summary_KNZ$spei12_category[-1] == "Extreme wet")
# 
# extreme_dry_wet <- c(CDR_edw, KBS_edw, KNZ_edw)
# 
# # Calculate the number of extreme wet then extreme wet years
# CDR_eww <- sum(spei_summary_CDR$spei12_category[-nrow(spei_summary_CDR)] == "Extreme wet" & spei_summary_CDR$spei12_category[-1] == "Extreme wet")
# KBS_eww <- sum(spei_summary_KBS$spei12_category[-nrow(spei_summary_KBS)] == "Extreme wet" & spei_summary_KBS$spei12_category[-1] == "Extreme wet")
# KNZ_eww <- sum(spei_summary_KNZ$spei12_category[-nrow(spei_summary_KNZ)] == "Extreme wet" & spei_summary_KNZ$spei12_category[-1] == "Extreme wet")
# 
# extreme_wet_wet <- c(CDR_eww, KBS_eww, KNZ_eww)
# 
# # Calculate the number of extreme dry then extreme dry years
# CDR_edd <- sum(spei_summary_CDR$spei12_category[-nrow(spei_summary_CDR)] == "Extreme dry" & spei_summary_CDR$spei12_category[-1] == "Extreme dry")
# KBS_edd <- sum(spei_summary_KBS$spei12_category[-nrow(spei_summary_KBS)] == "Extreme dry" & spei_summary_KBS$spei12_category[-1] == "Extreme dry")
# KNZ_edd <- sum(spei_summary_KNZ$spei12_category[-nrow(spei_summary_KNZ)] == "Extreme dry" & spei_summary_KNZ$spei12_category[-1] == "Extreme dry")
# 
# extreme_dry_dry <- c(CDR_edd, KBS_edd, KNZ_edd)
# 
# spei_table <- data.frame(sites, extreme_wet, moderate_wet, extreme_dry, moderate_dry, extreme_wet_dry, extreme_dry_wet, extreme_wet_wet, extreme_dry_dry)
# 
# # Not enough extreme events. Let's lump all moderate and extreme events into just wet and dry categories
# spei_lump <- spei_summary %>%
#   mutate(spei12_category = ifelse(as.character(spei12_category) == "Extreme dry", "Dry", as.character(spei12_category))) %>%
#   mutate(spei12_category = ifelse(as.character(spei12_category) == "Moderate dry", "Dry", as.character(spei12_category))) %>%
#   mutate(spei12_category = ifelse(as.character(spei12_category) == "Extreme wet", "Wet", as.character(spei12_category))) %>%
#   mutate(spei12_category = ifelse(as.character(spei12_category) == "Moderate wet", "Wet", as.character(spei12_category)))
# 
# # Calculate the number of moderate/extreme wet years
# CDR_w <- length(which(spei_lump$site == "CDR" & spei_lump$spei12_category == "Wet"))
# KBS_w <- length(which(spei_lump$site == "KBS" & spei_lump$spei12_category == "Wet"))
# KNZ_w <- length(which(spei_lump$site == "KNZ" & spei_lump$spei12_category == "Wet"))
# wet <- c(CDR_w, KBS_w, KNZ_w)
# 
# # Calculate the number of moderate/extreme dry years
# CDR_d <- length(which(spei_lump$site == "CDR" & spei_lump$spei12_category == "Dry"))
# KBS_d <- length(which(spei_lump$site == "KBS" & spei_lump$spei12_category == "Dry"))
# KNZ_d <- length(which(spei_lump$site == "KNZ" & spei_lump$spei12_category == "Dry"))
# dry <- c(CDR_d, KBS_d, KNZ_d)
# 
# # Calculate the number of wet then dry years
# spei_lump_CDR <- spei_lump %>%
#   filter(site == "CDR")
# CDR_wd <- sum(spei_lump_CDR$spei12_category[-nrow(spei_lump_CDR)] == "Wet" & spei_lump_CDR$spei12_category[-1] == "Dry")
# spei_lump_KBS <- spei_lump %>%
#   filter(site == "KBS")
# KBS_wd <- sum(spei_lump_KBS$spei12_category[-nrow(spei_lump_KBS)] == "Wet" & spei_lump_KBS$spei12_category[-1] == "Dry")
# spei_lump_KNZ <- spei_lump %>%
#   filter(site == "KNZ")
# KNZ_wd <- sum(spei_lump_KNZ$spei12_category[-nrow(spei_lump_KNZ)] == "Wet" & spei_lump_KNZ$spei12_category[-1] == "Dry")
# 
# wet_dry <- c(CDR_wd, KBS_wd, KNZ_wd)
# 
# # Calculate the number of wet then dry years
# spei_lump_CDR <- spei_lump %>%
#   filter(site == "CDR")
# CDR_wd <- sum(spei_lump_CDR$spei12_category[-nrow(spei_lump_CDR)] == "Wet" & spei_lump_CDR$spei12_category[-1] == "Dry")
# spei_lump_KBS <- spei_lump %>%
#   filter(site == "KBS")
# KBS_wd <- sum(spei_lump_KBS$spei12_category[-nrow(spei_lump_KBS)] == "Wet" & spei_lump_KBS$spei12_category[-1] == "Dry")
# spei_lump_KNZ <- spei_lump %>%
#   filter(site == "KNZ")
# KNZ_wd <- sum(spei_lump_KNZ$spei12_category[-nrow(spei_lump_KNZ)] == "Wet" & spei_lump_KNZ$spei12_category[-1] == "Dry")
# 
# wet_dry <- c(CDR_wd, KBS_wd, KNZ_wd)
# 
# # Calculate the number of dry then wet years
# CDR_dw <- sum(spei_lump_CDR$spei12_category[-nrow(spei_lump_CDR)] == "Dry" & spei_lump_CDR$spei12_category[-1] == "Wet")
# KBS_dw <- sum(spei_lump_KBS$spei12_category[-nrow(spei_lump_KBS)] == "Dry" & spei_lump_KBS$spei12_category[-1] == "Wet")
# KNZ_dw <- sum(spei_lump_KNZ$spei12_category[-nrow(spei_lump_KNZ)] == "Dry" & spei_lump_KNZ$spei12_category[-1] == "Wet")
# 
# dry_wet <- c(CDR_dw, KBS_dw, KNZ_dw)
# 
# # Calculate the number of wet then wet years
# CDR_ww <- sum(spei_lump_CDR$spei12_category[-nrow(spei_lump_CDR)] == "Wet" & spei_lump_CDR$spei12_category[-1] == "Wet")
# KBS_ww <- sum(spei_lump_KBS$spei12_category[-nrow(spei_lump_KBS)] == "Wet" & spei_lump_KBS$spei12_category[-1] == "Wet")
# KNZ_ww <- sum(spei_lump_KNZ$spei12_category[-nrow(spei_lump_KNZ)] == "Wet" & spei_lump_KNZ$spei12_category[-1] == "Wet")
# 
# wet_wet <- c(CDR_ww, KBS_ww, KNZ_ww)
# 
# # Calculate the number of dry then dry years
# CDR_dd <- sum(spei_lump_CDR$spei12_category[-nrow(spei_lump_CDR)] == "Dry" & spei_lump_CDR$spei12_category[-1] == "Dry")
# KBS_dd <- sum(spei_lump_KBS$spei12_category[-nrow(spei_lump_KBS)] == "Dry" & spei_lump_KBS$spei12_category[-1] == "Dry")
# KNZ_dd <- sum(spei_lump_KNZ$spei12_category[-nrow(spei_lump_KNZ)] == "Dry" & spei_lump_KNZ$spei12_category[-1] == "Dry")
# 
# dry_dry <- c(CDR_dd, KBS_dd, KNZ_dd)
# 
# spei_table_lump <- data.frame(sites, wet, dry, wet_dry, dry_wet, wet_wet, dry_dry)

# ## All biomass
# # Look at plant biomass averages by event category
# plot_spei %>%
#   ggplot(aes(x = spei12_category, y = plot_biomass)) + 
#   geom_sina() +
#   stat_summary(fun.data = mean_cl_boot, color = "red")
# 
# # Are there any significant differences?
# biomass.lm <- lmer(plot_biomass ~ spei12_category + (1|site) + (1|site:plot), data = plot_spei)
# anova(biomass.lm)
# emmeans(biomass.lm, list(pairwise ~ spei12_category), adjust = "tukey")
# plot_model(
#   biomass.lm,
#   type = "pred",
#   terms = c("spei12_category")
# )
# 
# # Look at plant biomass vs spei12
# plot_spei %>%
#   ggplot(aes(x = spei12, y = plot_biomass)) + 
#   geom_point()
# 
# plot_spei %>%
#   ggplot(aes(x = spei12, y = plot_biomass)) + 
#   geom_point() +
#   facet_wrap(~site)
# 
# # Is there a quadratic relationship?
# biomass_spei12.lm <- lmer(plot_biomass ~ spei12 + I(spei12^2) + (1|site) + (1|site:plot), data = plot_spei)
# summary(biomass_spei12.lm) # quadratic term significant
# biomass_spei12.lm2 <- lmer(plot_biomass ~ spei12 + (1|site) + (1|site:plot), data = plot_spei)
# 
# 
# plot_model(
#   biomass_spei12.lm,
#   type = "pred",
#   terms="spei12[all]",
#   show.data = TRUE
# )
# 
# # Why 1 missing values for plot biomass
# sum(is.na(plot_spei$plot_biomass))
# 
# 
# # Compare spei3, 6, 9, 12 models
# # They all appear to have potential quadratic terms # Quadratic makes sense conceptually
# plot_spei %>%
#   ggplot(aes(x = spei3, y = plot_biomass)) + 
#   geom_point()
# 
# plot_spei %>%
#   ggplot(aes(x = spei6, y = plot_biomass)) + 
#   geom_point()
# 
# plot_spei %>%
#   ggplot(aes(x = spei9, y = plot_biomass)) + 
#   geom_point()
# 
# biomass_spei9.lm <- lmer(plot_biomass ~ spei9 + I(spei9^2) + (1|site) + (1|site:plot), data = plot_spei)
# summary(biomass_spei9.lm) # Quadratic term significant
# biomass_spei9.lm2 <- lmer(plot_biomass ~ spei9 + (1|site) + (1|site:plot), data = plot_spei)
# 
# 
# biomass_spei6.lm <- lmer(plot_biomass ~ spei6 + I(spei6^2) + (1|site) + (1|site:plot), data = plot_spei)
# summary(biomass_spei6.lm)  # Quadratic term significant
# biomass_spei6.lm2 <- lmer(plot_biomass ~ spei6 + (1|site) + (1|site:plot), data = plot_spei)
# 
# 
# biomass_spei3.lm <- lmer(plot_biomass ~ spei3 + I(spei3^2) + (1|site) + (1|site:plot), data = plot_spei)
# summary(biomass_spei3.lm) # Quadratic term significant (is this just because of sample size?)
# biomass_spei3.lm2 <- lmer(plot_biomass ~ spei3 + (1|site) + (1|site:plot), data = plot_spei)
# 
# 
# # Compare models using AICc like Robinson paper
# AICctab(biomass_spei3.lm, biomass_spei3.lm2, biomass_spei6.lm, biomass_spei6.lm2, biomass_spei9.lm, biomass_spei9.lm2, biomass_spei12.lm, biomass_spei12.lm2) # spei9 lowest
# 
# # Compare model R2 like Robinson
# r.squaredGLMM(biomass_spei3.lm)
# r.squaredGLMM(biomass_spei3.lm2)
# r.squaredGLMM(biomass_spei6.lm)
# r.squaredGLMM(biomass_spei6.lm2)
# r.squaredGLMM(biomass_spei9.lm)
# r.squaredGLMM(biomass_spei9.lm2)
# r.squaredGLMM(biomass_spei12.lm)
# r.squaredGLMM(biomass_spei12.lm2)
# 
# # Plot spei 9 moel
# plot_model(
#   biomass_spei9.lm,
#   type = "pred",
#   terms="spei9[all]",
#   show.data = TRUE
# )
# 
# # Look at spei 9 model diagnostic plots
# simulationOutput <- simulateResiduals(fittedModel = biomass_spei9.lm, plot = F)
# plot(simulationOutput) # something going on here, but could be due to sample size
# 
# # Comparing spei by site
# # CDR
# biomass_cdr_spei12.lm <- lmer(plot_biomass ~ spei12 + I(spei12^2) + (1|plot), data = plot_spei, subset = site == "CDR")
# biomass_cdr_spei9.lm <- lmer(plot_biomass ~ spei9 + I(spei9^2) + (1|plot), data = plot_spei, subset = site == "CDR")
# biomass_cdr_spei6.lm <- lmer(plot_biomass ~ spei6 + I(spei6^2) + (1|plot), data = plot_spei, subset = site == "CDR")
# biomass_cdr_spei3.lm <- lmer(plot_biomass ~ spei3 + I(spei3^2) + (1|plot), data = plot_spei, subset = site == "CDR")
# 
# AICctab(biomass_cdr_spei12.lm, biomass_cdr_spei9.lm, biomass_cdr_spei6.lm, biomass_cdr_spei3.lm) # spei6 slightly better than 9
# 
# # KBS
# biomass_kbs_spei12.lm <- lmer(plot_biomass ~ spei12 + I(spei12^2) + (1|plot), data = plot_spei, subset = site == "KBS")
# biomass_kbs_spei9.lm <- lmer(plot_biomass ~ spei9 + I(spei9^2) + (1|plot), data = plot_spei, subset = site == "KBS")
# biomass_kbs_spei6.lm <- lmer(plot_biomass ~ spei6 + I(spei6^2) + (1|plot), data = plot_spei, subset = site == "KBS")
# biomass_kbs_spei3.lm <- lmer(plot_biomass ~ spei3 + I(spei3^2) + (1|plot), data = plot_spei, subset = site == "KBS")
# 
# AICctab(biomass_kbs_spei12.lm, biomass_kbs_spei9.lm, biomass_kbs_spei6.lm, biomass_kbs_spei3.lm) # spei9 better
# 
# # KNZ
# biomass_knz_spei12.lm <- lmer(plot_biomass ~ spei12 + I(spei12^2) + (1|plot), data = plot_spei, subset = site == "KNZ")
# biomass_knz_spei9.lm <- lmer(plot_biomass ~ spei9 + I(spei9^2) + (1|plot), data = plot_spei, subset = site == "KNZ")
# biomass_knz_spei6.lm <- lmer(plot_biomass ~ spei6 + I(spei6^2) + (1|plot), data = plot_spei, subset = site == "KNZ")
# biomass_knz_spei3.lm <- lmer(plot_biomass ~ spei3 + I(spei3^2) + (1|plot), data = plot_spei, subset = site == "KNZ")
# 
# AICctab(biomass_knz_spei12.lm, biomass_knz_spei9.lm, biomass_knz_spei6.lm, biomass_knz_spei3.lm) # spei9 better

# # Control biomass
# plot_control <- plot_trt %>%
#   filter(treatment == "control")
# 
# # Look at plant biomass averages by event category
# plot_control %>%
#   ggplot(aes(x = spei12_category, y = plot_biomass)) + 
#   geom_sina(alpha = 0.2) +
#   stat_summary(fun.data = mean_cl_boot, color = "red")
# 
# # Look at plant biomass vs spei12
# plot_control %>%
#   ggplot(aes(x = spei12, y = plot_biomass)) + 
#   geom_point(alpha=0.2) +
#   geom_smooth() +
#   ggpubr::stat_cor()
# 
# plot_control %>%
#   ggplot(aes(x = spei12, y = plot_biomass)) + 
#   geom_point(alpha  = 0.2) +
#   facet_wrap(~site)
# 
# # Model comparison with spei12, 9, 6, 3 for ONLY control plots
# control_spei3.lm2 <- lmer(plot_biomass ~ spei3 + I(spei3^2) + (1|site) + (1|site:plot), data = plot_control)
# summary(control_spei3.lm2) # Quadratic term significant
# control_spei3.lm <- lmer(plot_biomass ~ spei3 + (1|site) + (1|site:plot), data = plot_control)
# simres <- simulateResiduals(control_spei3.lm2)
# plot(simres)
# 
# control_spei6.lm2 <- lmer(plot_biomass ~ spei6 + I(spei6^2) + (1|site) + (1|site:plot), data = plot_control)
# summary(control_spei6.lm2) # Quadratic term significant
# control_spei6.lm <- lmer(plot_biomass ~ spei6 + (1|site) + (1|site:plot), data = plot_control)
# 
# control_spei9.lm2 <- lmer(plot_biomass ~ spei9 + I(spei9^2) + (1|site) + (1|site:plot), data = plot_control)
# summary(control_spei9.lm2) # Quadratic term significant
# control_spei9.lm <- lmer(plot_biomass ~ spei9 + (1|site) + (1|site:plot), data = plot_control)
# 
# control_spei12.lm2 <- lmer(plot_biomass ~ spei12 + I(spei12^2) + (1|site) + (1|site:plot), data = plot_control)
# summary(control_spei12.lm2) # Quadratic term significant
# control_spei12.lm <- lmer(plot_biomass ~ spei12 + (1|site) + (1|site:plot), data = plot_control)
# 
# AICctab(control_spei3.lm2, control_spei3.lm, control_spei6.lm2, control_spei6.lm, control_spei9.lm2, control_spei9.lm, control_spei12.lm2, control_spei12.lm)
# 
# # SPEI 6 with a quadratic term best model when all sites in one model
# 
# # Compare model R2 like Robinson
# r.squaredGLMM(control_spei3.lm2)
# r.squaredGLMM(control_spei3.lm)
# r.squaredGLMM(control_spei6.lm2) # best R2m
# r.squaredGLMM(control_spei6.lm)
# r.squaredGLMM(control_spei9.lm2)
# r.squaredGLMM(control_spei9.lm)
# r.squaredGLMM(control_spei12.lm2)
# r.squaredGLMM(control_spei12.lm)
# 
# # Plot the quadratic SPEI6 model
# plot_model(
#   control_spei6.lm2,
#   type = "pred",
#   terms="spei6[all]",
#   show.data = TRUE
# )
# 


# # Try logging biomass and adding more random effects
# # Model comparison with spei12, 9, 6, 3 for ONLY control plots
# nc_spei3.lm2.ln <- lmer(log1p(plot_biomass) ~ spei3 + nitrogen + I(spei3^2) + (1|site/experiment/uniqueid) + (1|year), data = plot_n_sub)
# summary(nc_spei3.lm2.ln) # Quadratic term significant
# nc_spei3.lm.ln <- lmer(log1p(plot_biomass) ~ spei3 + nitrogen + (1|site/experiment/uniqueid) + (1|year), data = plot_n_sub)
# simres <- simulateResiduals(nc_spei3.lm.ln)
# plot(simres) # Better 
# 
# nc_spei6.lm2.ln <- lmer(log1p(plot_biomass) ~ spei6 + I(spei6^2) + nitrogen + (1|site/experiment/uniqueid) + (1|year), data = plot_n_sub, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) # Failed to converge with defaults
# summary(nc_spei6.lm2.ln) # Quadratic term not significant
# nc_spei6.lm.ln <- lmer(log1p(plot_biomass) ~ spei6 + nitrogen + (1|site/experiment/uniqueid) + (1|year), data = plot_n_sub)
# simres <- simulateResiduals(nc_spei6.lm2.ln)
# plot(simres) # Better 
# 
# nc_spei9.lm2.ln <- lmer(log1p(plot_biomass) ~ spei9 + I(spei9^2) + nitrogen + (1|site/experiment/uniqueid) + (1|year), data = plot_n_sub)
# summary(nc_spei9.lm2.ln) # Quadratic term marginally significant
# nc_spei9.lm.ln <- lmer(log1p(plot_biomass) ~ spei9 + nitrogen + (1|site/experiment/uniqueid) + (1|year), data = plot_n_sub)
# summary(nc_spei9.lm.ln)
# simres <- simulateResiduals(nc_spei9.lm2.ln)
# plot(simres) # Better 
# simres <- simulateResiduals(nc_spei9.lm.ln)
# plot(simres) # Okay 
# 
# nc_spei12.lm2.ln <- lmer(log1p(plot_biomass) ~ spei12 + I(spei12^2) + nitrogen + (1|site/experiment/uniqueid) + (1|year), data = plot_n_sub)
# summary(nc_spei12.lm2.ln)
# nc_spei12.lm.ln <- lmer(log1p(plot_biomass) ~ spei12 + nitrogen + (1|site/experiment/uniqueid) + (1|year), data = plot_n_sub) 
# simres <- simulateResiduals(nc_spei12.lm2.ln)
# plot(simres) # Better 
# 
# AICctab(nc_spei3.lm2.ln, nc_spei3.lm.ln, nc_spei6.lm2.ln, nc_spei6.lm.ln, nc_spei9.lm2.ln, nc_spei9.lm.ln, nc_spei12.lm2.ln, nc_spei12.lm.ln)
# 
# # Compare model R2 like Robinson
# r.squaredGLMM(nc_spei3.lm2.ln)
# r.squaredGLMM(nc_spei3.lm.ln)
# r.squaredGLMM(nc_spei6.lm2.ln)
# r.squaredGLMM(nc_spei6.lm.ln)
# r.squaredGLMM(nc_spei9.lm2.ln)
# r.squaredGLMM(nc_spei9.lm.ln)
# r.squaredGLMM(nc_spei12.lm2.ln) # best R2m
# r.squaredGLMM(nc_spei12.lm.ln)  # second best R2m
# 
# # Plot the SPEI9 model
# plot_model(
#   nc_spei9.lm.ln,
#   type = "pred",
#   terms= c("spei9", "nitrogen"),
#   show.data = TRUE
# )
# 
# # Log response of SPEI
# # Using SPEI3
# plot_control_lrr3 <- plot_control %>%
#   group_by(uniqueid) %>%
#   arrange(uniqueid, year) %>%
#   mutate(prior_year_biomass = lag(plot_biomass),
#          prior_year_type = lag(spei3_category)) %>%
#   ungroup()
# normal_biomass3 <- plot_control_lrr3  %>%
#   group_by(uniqueid) %>%
#   filter(spei3_category == "Normal" | spei3_category == "Moderate wet" | spei3_category == "Moderate dry") %>%
#   summarize(average_normal_biomass = mean(plot_biomass, na.rm=T))
# plot_control_lrr3 <- left_join(plot_control_lrr3, normal_biomass3)
# 
# lrr3 <- plot_control_lrr3 %>%
#   filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
#   filter(spei3_category == "Extreme wet" | spei3_category == "Extreme dry") %>%
#   mutate(LRR = log(plot_biomass / average_normal_biomass))  # Calculate the log response ratio
# 
# # Take absolute value of SPEI3
# lrr3 <- lrr3 %>%
#   mutate(abs_spei3 = abs(spei3))
# 
# lrr.lm3 <- lmer(LRR ~ abs_spei3*spei3_category + (1|site/experiment/uniqueid) + (1|year), data = lrr3)
# summary(lrr.lm3)
# simres <- simulateResiduals(lrr.lm3)
# plot(simres) # Weird
# 
# # Plot the LRR model
# plot_model(
#   lrr.lm3,
#   type = "pred",
#   terms= c("abs_spei3", "spei3_category"),
#   show.data = TRUE
# )
# 
# # Using SPEI6
# plot_control_lrr <- plot_control %>%
#   group_by(uniqueid) %>%
#   arrange(uniqueid, year) %>%
#   mutate(prior_year_biomass = lag(plot_biomass),
#          prior_year_type = lag(spei6_category)) %>%
#   ungroup()
# normal_biomass <- plot_control_lrr  %>%
#   group_by(uniqueid) %>%
#   filter(spei6_category == "Normal" | spei6_category == "Moderate wet" | spei6_category == "Moderate dry") %>%
#   summarize(average_normal_biomass = mean(plot_biomass, na.rm=T))
# plot_control_lrr <- left_join(plot_control_lrr, normal_biomass)
# 
# lrr <- plot_control_lrr %>%
#   filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
#   filter(spei6_category == "Extreme wet" | spei6_category == "Extreme dry") %>%
#   mutate(LRR = log(plot_biomass / average_normal_biomass))  # Calculate the log response ratio
# 
# # Take absolute value of SPEI6
# lrr <- lrr %>%
#   mutate(abs_spei6 = abs(spei6))
# 
# lrr.lm <- lmer(LRR ~ abs_spei6*spei6_category + (1|site/experiment/uniqueid) + (1|year), data = lrr)
# summary(lrr.lm)
# simres <- simulateResiduals(lrr.lm)
# plot(simres) # Weird
# 
# # Plot the LRR model
# plot_model(
#   lrr.lm,
#   type = "pred",
#   terms= c("abs_spei6", "spei6_category"),
#   show.data = TRUE
# )
# 
# # Using SPEI9
# plot_control_lrr9 <- plot_control %>%
#   group_by(uniqueid) %>%
#   arrange(uniqueid, year) %>%
#   mutate(prior_year_biomass = lag(plot_biomass),
#          prior_year_type = lag(spei9_category)) %>%
#   ungroup()
# normal_biomass9 <- plot_control_lrr9  %>%
#   group_by(uniqueid) %>%
#   filter(spei9_category == "Normal" | spei9_category == "Moderate wet" | spei9_category == "Moderate dry") %>%
#   summarize(average_normal_biomass = mean(plot_biomass, na.rm=T))
# plot_control_lrr9 <- left_join(plot_control_lrr9, normal_biomass9)
# 
# lrr9 <- plot_control_lrr9 %>%
#   filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
#   filter(spei9_category == "Extreme wet" | spei9_category == "Extreme dry") %>%
#   mutate(LRR = log(plot_biomass / average_normal_biomass))  # Calculate the log response ratio
# 
# # Take absolute value of SPEI9
# lrr9 <- lrr9 %>%
#   mutate(abs_spei9 = abs(spei9))
# 
# lrr.lm9 <- lmer(LRR ~ abs_spei9*spei9_category + (1|site/experiment/uniqueid) + (1|year), data = lrr9)
# summary(lrr.lm9)
# Anova(lrr.lm9, type = "III")
# simres <- simulateResiduals(lrr.lm9)
# plot(simres)
# 
# # Plot the LRR model
# plot_model(
#   lrr.lm9,
#   type = "pred",
#   terms= c("abs_spei9", "spei9_category"),
#   show.data = TRUE
# )
# 
# # Using SPEI12
# plot_control_lrr12 <- plot_control %>%
#   group_by(uniqueid) %>%
#   arrange(uniqueid, year) %>%
#   mutate(prior_year_biomass = lag(plot_biomass),
#          prior_year_type = lag(spei12_category)) %>%
#   ungroup()
# normal_biomass12 <- plot_control_lrr12  %>%
#   group_by(uniqueid) %>%
#   filter(spei12_category == "Normal" | spei12_category == "Moderate wet" | spei12_category == "Moderate dry") %>%
#   summarize(average_normal_biomass = mean(plot_biomass, na.rm=T))
# plot_control_lrr12 <- left_join(plot_control_lrr12, normal_biomass12)
# 
# lrr12 <- plot_control_lrr12 %>%
#   filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
#   filter(spei12_category == "Extreme wet" | spei12_category == "Extreme dry") %>%
#   mutate(LRR = log(plot_biomass / average_normal_biomass))  # Calculate the log response ratio
# 
# # Take absolute value of SPEI6
# lrr12 <- lrr12 %>%
#   mutate(abs_spei12 = abs(spei12))
# 
# lrr.lm12 <- lmer(LRR ~ abs_spei12*spei12_category + (1|site/experiment/uniqueid) + (1|year), data = lrr12)
# summary(lrr.lm12)
# simres <- simulateResiduals(lrr.lm12)
# plot(simres) # Weird
# 
# # Plot the LRR model
# plot_model(
#   lrr.lm12,
#   type = "pred",
#   terms= c("abs_spei12", "spei12_category"),
#   show.data = TRUE
# )
# Using SPEI3
# plot_n_lrr3 <- plot_n %>%
#   group_by(uniqueid) %>%
#   arrange(uniqueid, year) %>%
#   mutate(prior_year_biomass = lag(plot_biomass),
#          prior_year_type = lag(spei3_category)) %>%
#   ungroup()
# normal_biomass3n <- plot_n_lrr3  %>%
#   group_by(uniqueid) %>%
#   filter(spei3_category == "Normal" | spei3_category == "Moderate wet" | spei3_category == "Moderate dry") %>%
#   summarize(average_normal_biomass = mean(plot_biomass, na.rm=T))
# plot_n_lrr3 <- left_join(plot_n_lrr3, normal_biomass3n)
# 
# lrr3n <- plot_n_lrr3 %>%
#   filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
#   filter(spei3_category == "Extreme wet" | spei3_category == "Extreme dry") %>%
#   mutate(LRR = log(plot_biomass / average_normal_biomass))  # Calculate the log response ratio
# 
# # Take absolute value of SPEI3
# lrr3n <- lrr3n %>%
#   mutate(abs_spei3 = abs(spei3))
# 
# lrr.lm3n <- lmer(LRR ~ abs_spei3*spei3_category + (1|site/experiment/uniqueid) + (1|year), data = lrr3n)
# summary(lrr.lm3n)
# simres <- simulateResiduals(lrr.lm3n)
# plot(simres) # Weird
# 
# # Plot the LRR model
# plot_model(
#   lrr.lm3n,
#   type = "pred",
#   terms= c("abs_spei3", "spei3_category"),
#   show.data = TRUE
# )
# 
# # Using SPEI6
# plot_n_lrr6 <- plot_n %>%
#   group_by(uniqueid) %>%
#   arrange(uniqueid, year) %>%
#   mutate(prior_year_biomass = lag(plot_biomass),
#          prior_year_type = lag(spei6_category)) %>%
#   ungroup()
# normal_biomass6n <- plot_n_lrr6  %>%
#   group_by(uniqueid) %>%
#   filter(spei6_category == "Normal" | spei6_category == "Moderate wet" | spei6_category == "Moderate dry") %>%
#   summarize(average_normal_biomass = mean(plot_biomass, na.rm=T))
# plot_n_lrr6 <- left_join(plot_n_lrr6, normal_biomass6n)
# 
# lrr6n <- plot_n_lrr6 %>%
#   filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
#   filter(spei6_category == "Extreme wet" | spei6_category == "Extreme dry") %>%
#   mutate(LRR = log(plot_biomass / average_normal_biomass))  # Calculate the log response ratio
# 
# # Take absolute value of SPEI6
# lrr6n <- lrr6n %>%
#   mutate(abs_spei6 = abs(spei6))
# 
# # Remove zero value in LRR
# lrr6n <- lrr6n %>% filter(plot_biomass != 0)
# 
# lrr.lm6n <- lmer(LRR ~ abs_spei6*spei6_category + (1|site/experiment/uniqueid) + (1|year), data = lrr6n)
# summary(lrr.lm6n)
# simres <- simulateResiduals(lrr.lm6n)
# plot(simres) # Weird
# 
# # Plot the LRR model
# plot_model(
#   lrr.lm6n,
#   type = "pred",
#   terms= c("abs_spei6", "spei6_category"),
#   show.data = TRUE
# )
# 
# # Using SPEI9
# plot_n_lrr9 <- plot_n %>%
#   group_by(uniqueid) %>%
#   arrange(uniqueid, year) %>%
#   mutate(prior_year_biomass = lag(plot_biomass),
#          prior_year_type = lag(spei9_category)) %>%
#   ungroup()
# normal_biomass9n <- plot_n_lrr9  %>%
#   group_by(uniqueid) %>%
#   filter(spei9_category == "Normal" | spei9_category == "Moderate wet" | spei9_category == "Moderate dry") %>%
#   summarize(average_normal_biomass = mean(plot_biomass, na.rm=T))
# plot_n_lrr9 <- left_join(plot_n_lrr9, normal_biomass9n)
# 
# lrr9n <- plot_n_lrr9 %>%
#   filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
#   filter(spei9_category == "Extreme wet" | spei9_category == "Extreme dry") %>%
#   mutate(LRR = log(plot_biomass / average_normal_biomass))  # Calculate the log response ratio
# 
# # Take absolute value of SPEI9
# lrr9n <- lrr9n %>%
#   mutate(abs_spei9 = abs(spei9))
# 
# # Remove zero value in LRR
# lrr9n <- lrr9n %>% filter(plot_biomass != 0)
# 
# lrr.lm9n <- lmer(LRR ~ abs_spei9*spei9_category + (1|site/experiment/uniqueid) + (1|year), data = lrr9n)
# summary(lrr.lm9n)
# simres <- simulateResiduals(lrr.lm9n)
# plot(simres) # Weird
# 
# # Plot the LRR model
# plot_model(
#   lrr.lm9n,
#   type = "pred",
#   terms= c("abs_spei9", "spei9_category"),
#   show.data = TRUE
# )
# 
# # Using SPEI12
# plot_n_lrr12 <- plot_n %>%
#   group_by(uniqueid) %>%
#   arrange(uniqueid, year) %>%
#   mutate(prior_year_biomass = lag(plot_biomass),
#          prior_year_type = lag(spei12_category)) %>%
#   ungroup()
# normal_biomass12n <- plot_n_lrr12  %>%
#   group_by(uniqueid) %>%
#   filter(spei12_category == "Normal" | spei12_category == "Moderate wet" | spei12_category == "Moderate dry") %>%
#   summarize(average_normal_biomass = mean(plot_biomass, na.rm=T))
# plot_n_lrr12 <- left_join(plot_n_lrr12, normal_biomass12n)
# 
# lrr12n <- plot_n_lrr12 %>%
#   filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
#   filter(spei12_category == "Extreme wet" | spei12_category == "Extreme dry") %>%
#   mutate(LRR = log(plot_biomass / average_normal_biomass))  # Calculate the log response ratio
# 
# # Take absolute value of SPEI12
# lrr12n <- lrr12n %>%
#   mutate(abs_spei12 = abs(spei12))
# 
# # Remove zero value in LRR
# lrr12n <- lrr12n %>% filter(plot_biomass != 0)
# 
# lrr.lm12n <- lmer(LRR ~ abs_spei12*spei12_category + (1|site/experiment/uniqueid) + (1|year), data = lrr12n, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) # Failed to converge with defaults
# summary(lrr.lm12n)
# simres <- simulateResiduals(lrr.lm12n)
# plot(simres) # Weird
# 
# # Plot the LRR model
# plot_model(
#   lrr.lm12n,
#   type = "pred",
#   terms= c("abs_spei12", "spei12_category"),
#   show.data = TRUE
# )
# 
# 
# # LRR with nitrogen yes/no as a predictor
# plot_nc <- plot_trt %>%
#   filter(nutrients_added == "NPK+" | nutrients_added == "N" | nutrients_added == "NP" | nutrients_added == "NPK" | nutrients_added == "NK" | nutrients_added == "NPK+Fence" | nutrients_added == "none" | nutrients_added == "no_fertilizer")
# 
# # Add column with nitrogen yes/no
# plot_nc <- plot_nc %>%
#   mutate(nitrogen = dplyr::recode(nutrients_added, "NPK+" = "N", "NP" = "N", "NPK" = "N", "NK" = "N", "NPK+Fence" = "N", "none" = "no_fertilizer"))
# 
# # Using SPEI3
# plot_nc_lrr3 <- plot_nc %>%
#   group_by(uniqueid) %>%
#   arrange(uniqueid, year) %>%
#   mutate(prior_year_biomass = lag(plot_biomass),
#          prior_year_type = lag(spei3_category)) %>%
#   ungroup()
# normal_biomass3nc <- plot_nc_lrr3  %>%
#   group_by(uniqueid) %>%
#   filter(spei3_category == "Normal" | spei3_category == "Moderate wet" | spei3_category == "Moderate dry") %>%
#   summarize(average_normal_biomass = mean(plot_biomass, na.rm=T))
# plot_nc_lrr3 <- left_join(plot_nc_lrr3, normal_biomass3nc)
# 
# lrr3nc <- plot_nc_lrr3 %>%
#   filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
#   filter(spei3_category == "Extreme wet" | spei3_category == "Extreme dry") %>%
#   mutate(LRR = log(plot_biomass / average_normal_biomass))  # Calculate the log response ratio
# 
# # Take absolute value of SPEI3
# lrr3nc <- lrr3nc %>%
#   mutate(abs_spei3 = abs(spei3))
# 
# lrr.lm3nc <- lmer(LRR ~ abs_spei3*spei3_category*nitrogen + (1|site/experiment/uniqueid) + (1|year), data = lrr3nc)
# summary(lrr.lm3nc)
# Anova(lrr.lm3nc, type = "III")
# simres <- simulateResiduals(lrr.lm3nc)
# plot(simres)
# 
# # Plot the LRR model
# plot_model(
#   lrr.lm3nc,
#   type = "pred",
#   terms= c("abs_spei3", "spei3_category", "nitrogen"),
#   show.data = TRUE
# )
# 
# # Using SPEI6
# plot_nc_lrr6 <- plot_nc %>%
#   group_by(uniqueid) %>%
#   arrange(uniqueid, year) %>%
#   mutate(prior_year_biomass = lag(plot_biomass),
#          prior_year_type = lag(spei6_category)) %>%
#   ungroup()
# normal_biomass6nc <- plot_nc_lrr6  %>%
#   group_by(uniqueid) %>%
#   filter(spei6_category == "Normal" | spei6_category == "Moderate wet" | spei6_category == "Moderate dry") %>%
#   summarize(average_normal_biomass = mean(plot_biomass, na.rm=T))
# plot_nc_lrr6 <- left_join(plot_nc_lrr6, normal_biomass6nc)
# 
# lrr6nc <- plot_nc_lrr6 %>%
#   filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
#   filter(spei6_category == "Extreme wet" | spei6_category == "Extreme dry") %>%
#   mutate(LRR = log(plot_biomass / average_normal_biomass))  # Calculate the log response ratio
# 
# # Take absolute value of SPEI6
# lrr6nc <- lrr6nc %>%
#   mutate(abs_spei6 = abs(spei6))
# 
# # Remove zero value in LRR
# lrr6nc <- lrr6nc %>% filter(plot_biomass != 0)
# 
# lrr.lm6nc <- lmer(LRR ~ abs_spei6*spei6_category*nitrogen + (1|site/experiment/uniqueid) + (1|year), data = lrr6nc)
# summary(lrr.lm6nc)
# Anova(lrr.lm6nc, type = "III")
# simres <- simulateResiduals(lrr.lm6nc)
# plot(simres)
# 
# # Plot the LRR model
# plot_model(
#   lrr.lm6nc,
#   type = "pred",
#   terms= c("abs_spei6", "spei6_category", "nitrogen"),
#   show.data = TRUE
# )
# 
# # Using SPEI9
# plot_nc_lrr9 <- plot_nc %>%
#   group_by(uniqueid) %>%
#   arrange(uniqueid, year) %>%
#   mutate(prior_year_biomass = lag(plot_biomass),
#          prior_year_type = lag(spei9_category)) %>%
#   ungroup()
# normal_biomass9nc <- plot_nc_lrr9  %>%
#   group_by(uniqueid) %>%
#   filter(spei9_category == "Normal" | spei9_category == "Moderate wet" | spei9_category == "Moderate dry") %>%
#   summarize(average_normal_biomass = mean(plot_biomass, na.rm=T))
# plot_nc_lrr9 <- left_join(plot_nc_lrr9, normal_biomass9nc)
# 
# lrr9nc <- plot_nc_lrr9 %>%
#   filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
#   filter(spei9_category == "Extreme wet" | spei9_category == "Extreme dry") %>%
#   mutate(LRR = log(plot_biomass / average_normal_biomass))  # Calculate the log response ratio
# 
# # Take absolute value of SPEI9
# lrr9nc <- lrr9nc %>%
#   mutate(abs_spei9 = abs(spei9))
# 
# # Remove zero value in LRR
# lrr9nc <- lrr9nc %>% filter(plot_biomass != 0)
# 
# lrr.lm9nc <- lmer(LRR ~ abs_spei9*spei9_category*nitrogen + (1|site/experiment/uniqueid) + (1|year), data = lrr9nc)
# summary(lrr.lm9nc)
# Anova(lrr.lm9nc, type = "III")
# simres <- simulateResiduals(lrr.lm9nc)
# plot(simres)
# 
# # Plot the LRR model
# plot_model(
#   lrr.lm9nc,
#   type = "pred",
#   terms= c("abs_spei9", "spei9_category", "nitrogen"),
#   show.data = TRUE
# )
# 
# # Using SPEI12
# plot_nc_lrr12 <- plot_nc %>%
#   group_by(uniqueid) %>%
#   arrange(uniqueid, year) %>%
#   mutate(prior_year_biomass = lag(plot_biomass),
#          prior_year_type = lag(spei12_category)) %>%
#   ungroup()
# normal_biomass12nc <- plot_nc_lrr12  %>%
#   group_by(uniqueid) %>%
#   filter(spei12_category == "Normal" | spei12_category == "Moderate wet" | spei12_category == "Moderate dry") %>%
#   summarize(average_normal_biomass = mean(plot_biomass, na.rm=T))
# plot_nc_lrr12 <- left_join(plot_nc_lrr12, normal_biomass12nc)
# 
# lrr12nc <- plot_nc_lrr12 %>%
#   filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
#   filter(spei12_category == "Extreme wet" | spei12_category == "Extreme dry") %>%
#   mutate(LRR = log(plot_biomass / average_normal_biomass))  # Calculate the log response ratio
# 
# # Take absolute value of SPEI12
# lrr12nc <- lrr12nc %>%
#   mutate(abs_spei12 = abs(spei12))
# 
# # Remove zero value in LRR
# lrr12nc <- lrr12nc %>% filter(plot_biomass != 0)
# 
# lrr.lm12nc <- lmer(LRR ~ abs_spei12*spei12_category*nitrogen + (1|site/experiment/uniqueid) + (1|year), data = lrr12nc)
# summary(lrr.lm12nc)
# Anova(lrr.lm12nc, type = "III")
# simres <- simulateResiduals(lrr.lm12nc)
# plot(simres)
# 
# # Plot the LRR model
# plot_model(
#   lrr.lm12nc,
#   type = "pred",
#   terms= c("abs_spei12", "spei12_category", "nitrogen"),
#   show.data = TRUE
# )
# 
# # LRR Nitrogen subset (no treatments)
# # Using SPEI3
# plot_sub_lrr3 <- plot_n_sub %>%
#   group_by(uniqueid) %>%
#   arrange(uniqueid, year) %>%
#   mutate(prior_year_biomass = lag(plot_biomass),
#          prior_year_type = lag(spei3_category)) %>%
#   ungroup()
# normal_biomass3sub <- plot_sub_lrr3  %>%
#   group_by(uniqueid) %>%
#   filter(spei3_category == "Normal" | spei3_category == "Moderate wet" | spei3_category == "Moderate dry") %>%
#   summarize(average_normal_biomass = mean(plot_biomass, na.rm=T))
# plot_sub_lrr3 <- left_join(plot_sub_lrr3, normal_biomass3sub)
# 
# lrr3sub <- plot_sub_lrr3 %>%
#   filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
#   filter(spei3_category == "Extreme wet" | spei3_category == "Extreme dry") %>%
#   mutate(LRR = log(plot_biomass / average_normal_biomass))  # Calculate the log response ratio
# 
# # Take absolute value of SPEI3
# lrr3sub <- lrr3sub %>%
#   mutate(abs_spei3 = abs(spei3))
# 
# lrr.lm3sub <- lmer(LRR ~ abs_spei3*spei3_category*nitrogen + (1|site/experiment/uniqueid) + (1|year), data = lrr3sub)
# summary(lrr.lm3sub)
# Anova(lrr.lm3sub, type = "III")
# simres <- simulateResiduals(lrr.lm3sub)
# plot(simres)
# 
# # Plot the LRR model
# plot_model(
#   lrr.lm3sub,
#   type = "pred",
#   terms= c("abs_spei3", "spei3_category", "nitrogen"),
#   show.data = TRUE
# )
# 
# # Maowei's way for SPEI3
# lrr.lm3sub2 <- lmer(LRR ~ spei3_category*nitrogen + (1|site/experiment/uniqueid) + (1|year), data = lrr3sub)
# summary(lrr.lm3sub2)
# simres <- simulateResiduals(lrr.lm3sub2)
# plot(simres)
# 
# # Plot the LRR model
# plot_model(
#   lrr.lm3sub2,
#   type = "pred",
#   terms= c("spei3_category", "nitrogen"),
#   show.data = TRUE
# )
# 
# 
# # Using SPEI6
# plot_sub_lrr6 <- plot_n_sub %>%
#   group_by(uniqueid) %>%
#   arrange(uniqueid, year) %>%
#   mutate(prior_year_biomass = lag(plot_biomass),
#          prior_year_type = lag(spei6_category)) %>%
#   ungroup()
# normal_biomass6sub <- plot_sub_lrr6  %>%
#   group_by(uniqueid) %>%
#   filter(spei6_category == "Normal" | spei6_category == "Moderate wet" | spei6_category == "Moderate dry") %>%
#   summarize(average_normal_biomass = mean(plot_biomass, na.rm=T))
# plot_sub_lrr6 <- left_join(plot_sub_lrr6, normal_biomass6sub)
# 
# lrr6sub <- plot_sub_lrr6 %>%
#   filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
#   filter(spei6_category == "Extreme wet" | spei6_category == "Extreme dry") %>%
#   mutate(LRR = log(plot_biomass / average_normal_biomass))  # Calculate the log response ratio
# 
# # Take absolute value of SPEI6
# lrr6sub <- lrr6sub %>%
#   mutate(abs_spei6 = abs(spei6))
# 
# lrr.lm6sub <- lmer(LRR ~ abs_spei6*spei6_category*nitrogen + (1|site/experiment/uniqueid) + (1|year), data = lrr6sub, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) # Failed to converge with defaults)
# summary(lrr.lm6sub)
# Anova(lrr.lm6sub, type = "III")
# simres <- simulateResiduals(lrr.lm6sub)
# plot(simres)
# 
# # Plot the LRR model
# plot_model(
#   lrr.lm6sub,
#   type = "pred",
#   terms= c("abs_spei6", "spei6_category", "nitrogen"),
#   show.data = TRUE
# )
# 
# # Using SPEI9
# plot_sub_lrr9 <- plot_n_sub %>%
#   group_by(uniqueid) %>%
#   arrange(uniqueid, year) %>%
#   mutate(prior_year_biomass = lag(plot_biomass),
#          prior_year_type = lag(spei9_category)) %>%
#   ungroup()
# normal_biomass9sub <- plot_sub_lrr9  %>%
#   group_by(uniqueid) %>%
#   filter(spei9_category == "Normal" | spei9_category == "Moderate wet" | spei9_category == "Moderate dry") %>%
#   summarize(average_normal_biomass = mean(plot_biomass, na.rm=T))
# plot_sub_lrr9 <- left_join(plot_sub_lrr9, normal_biomass9sub)
# 
# lrr9sub <- plot_sub_lrr9 %>%
#   # filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
#   filter(spei9_category == "Extreme wet" | spei9_category == "Extreme dry") %>%
#   mutate(LRR = log(plot_biomass / average_normal_biomass))  # Calculate the log response ratio
# 
# # Take absolute value of SPEI9
# lrr9sub <- lrr9sub %>%
#   mutate(abs_spei9 = abs(spei9))
# 
# lrr.lm9sub <- lmer(LRR ~ abs_spei9*spei9_category*nitrogen + (1|site/experiment/uniqueid) + (1|year), data = lrr9sub)
# summary(lrr.lm9sub)
# Anova(lrr.lm9sub, type = "III")
# simres <- simulateResiduals(lrr.lm9sub)
# plot(simres)
# 
# # Plot the LRR model
# plot_model(
#   lrr.lm9sub,
#   type = "pred",
#   terms= c("abs_spei9", "spei9_category", "nitrogen"),
#   show.data = TRUE
# )
# 
# # Simple model
# lrr.lm9sub2 <- lmer(LRR ~ spei9_category*nitrogen + (1|site/experiment/uniqueid) + (1|year), data = lrr9sub)
# summary(lrr.lm9sub2)
# emm_biomass <- emmeans(lrr.lm9sub2, pairwise ~ spei9_category * nitrogen, infer = T)
# summary(emm_biomass)
# simres <- simulateResiduals(lrr.lm9sub2)
# plot(simres)
# 
# plot_model(
#   lrr.lm9sub2,
#   type = "pred",
#   terms= c("spei9_category", "nitrogen"),
#   show.data = F
# ) + theme_bw()
# 
# ## Final LRR plot for ANPP
# cols <- c("Extreme dry" = "#E41A1C", "Extreme wet" = "#377EB8")
# lrr9sub %>%
#   ggplot(aes(x = spei9_category, y = LRR, col = spei9_category)) +
#   stat_summary(fun.data = mean_cl_boot, position = position_dodge(0.2), aes(shape = nitrogen)) + 
#   geom_hline(yintercept=0, linetype='dashed', col = 'black')+
#   theme_bw() +
#   scale_color_manual(values = cols)
# anpp_plot <- lrr9sub %>%
#   mutate(nitrogen = relevel(as.factor(nitrogen), 'no_fertilizer', 'nitrogen')) %>%
#   ggplot(aes(x = spei9_category, y = LRR, col = spei9_category)) +
#   stat_summary(fun.data = mean_cl_boot, position = position_dodge(0.2), aes(shape = nitrogen)) + 
#   geom_hline(yintercept=0, linetype='dashed', col = 'black')+
#   theme_bw() +
#   labs(x = "Event type", y = "ANPP LRR") +
#   theme(legend.position="none") +
#   scale_color_manual(values = cols)
# 
# # LRR for richness
# plot_sub_lrr9_rich <- plot_n_sub %>%
#   group_by(uniqueid) %>%
#   arrange(uniqueid, year) %>%
#   mutate(prior_year_richness = lag(Richness),
#          prior_year_type = lag(spei9_category)) %>%
#   ungroup()
# normal_biomass9sub_rich <- plot_sub_lrr9_rich  %>%
#   group_by(uniqueid) %>%
#   filter(spei9_category == "Normal" | spei9_category == "Moderate wet" | spei9_category == "Moderate dry") %>%
#   summarize(average_normal_rich = mean(Richness, na.rm=T))
# plot_sub_lrr9_rich <- left_join(plot_sub_lrr9_rich, normal_biomass9sub_rich)
# 
# lrr9sub_rich <- plot_sub_lrr9_rich %>%
#   # filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
#   filter(spei9_category == "Extreme wet" | spei9_category == "Extreme dry") %>%
#   mutate(LRR = log(Richness / average_normal_rich))  # Calculate the log response ratio
# 
# lrr.lm9sub_rich <- lmer(LRR ~ spei9_category*nitrogen + (1|site/experiment/uniqueid) + (1|year), data = lrr9sub_rich)
# summary(lrr.lm9sub_rich)
# emm_rich <- emmeans(lrr.lm9sub_rich, pairwise ~ spei9_category * nitrogen, infer = T)
# summary(emm_rich)
# simres <- simulateResiduals(lrr.lm9sub_rich)
# plot(simres)
# 
# ## Final LRR plot for richness
# lrr9sub_rich %>%
#   ggplot(aes(x = spei9_category, y = LRR, col = nitrogen)) +
#   stat_summary(fun.data = mean_cl_boot, position = position_dodge(0.2)) + 
#   geom_hline(yintercept=0, linetype='dashed', col = 'black')+
#   theme_bw()
# rich_plot <- lrr9sub_rich %>%
#   mutate(nitrogen = relevel(as.factor(nitrogen), 'no_fertilizer', 'nitrogen')) %>%
#   ggplot(aes(x = spei9_category, y = LRR, col = spei9_category)) +
#   stat_summary(fun.data = mean_cl_boot, position = position_dodge(0.2), aes(shape = nitrogen)) + 
#   geom_hline(yintercept=0, linetype='dashed', col = 'black')+
#   theme_bw() +
#   labs(x = "Event type", y = "Richness LRR") +
#   theme(legend.position="none") +
#   scale_color_manual(values = cols)
# 
# plot_model(
#   lrr.lm9sub_rich,
#   type = "pred",
#   terms= c("spei9_category", "nitrogen"),
#   show.data = F
# ) + theme_bw()
# 
# # LRR for dominance (NO ZEROS!!!)
# plot_sub_lrr9_dom <- plot_n_sub %>%
#   group_by(uniqueid) %>%
#   arrange(uniqueid, year) %>%
#   mutate(prior_year_dom = lag(dominant_relative_abund),
#          prior_year_type = lag(spei9_category)) %>%
#   ungroup()
# normal_biomass9sub_dom <- plot_sub_lrr9_dom  %>%
#   group_by(uniqueid) %>%
#   filter(spei9_category == "Normal" | spei9_category == "Moderate wet" | spei9_category == "Moderate dry") %>%
#   summarize(average_normal_dom = mean(dominant_relative_abund, na.rm=T))
# plot_sub_lrr9_dom <- left_join(plot_sub_lrr9_dom, normal_biomass9sub_dom)
# 
# lrr9sub_dom <- plot_sub_lrr9_dom %>%
#   # filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
#   filter(spei9_category == "Extreme wet" | spei9_category == "Extreme dry") %>%
#   mutate(LRR = log(dominant_relative_abund / average_normal_dom))  # Calculate the log response ratio
# 
# lrr.lm9sub_dom <- lmer(LRR ~ spei9_category*nitrogen + (1|site/experiment/uniqueid) + (1|year), data = lrr9sub_dom)
# summary(lrr.lm9sub_dom)
# emm_dom <- emmeans(lrr.lm9sub_dom, pairwise ~ spei9_category * nitrogen, infer = T)
# summary(emm_dom)
# simres <- simulateResiduals(lrr.lm9sub_dom)
# plot(simres)
# 
# ## Final LRR plot for dominance
# lrr9sub_dom %>%
#   ggplot(aes(x = spei9_category, y = LRR, col = nitrogen)) +
#   stat_summary(fun.data = mean_cl_boot, position = position_dodge(0.2)) + 
#   geom_hline(yintercept=0, linetype='dashed', col = 'black')+
#   theme_bw()
# dom_plot <- lrr9sub_dom %>%
#   mutate(nitrogen = relevel(as.factor(nitrogen), 'no_fertilizer', 'nitrogen')) %>%
#   ggplot(aes(x = spei9_category, y = LRR, col = spei9_category)) +
#   stat_summary(fun.data = mean_cl_boot, position = position_dodge(0.2), aes(shape = nitrogen)) + 
#   geom_hline(yintercept=0, linetype='dashed', col = 'black')+
#   theme_bw() +
#   labs(x = "Event type", y = "Dominance LRR", shape = "Fertilizer") +
#   scale_color_manual(values = cols, guide = "none") +
#   scale_shape_manual(values = c(19, 17), labels = c('no nutrients', 'nutrients')) +
#   theme(legend.title=element_blank())
# 
# plot_model(
#   lrr.lm9sub_dom,
#   type = "pred",
#   terms= c("spei9_category", "nitrogen"),
#   show.data = F
# ) + theme_bw()
# 
# # LRR for dominance (Top 6 species) 
# # Andropogon gerardii/Andropogon gerardi (typo), Elymus repens, Poa pratensis, Schizachyrium scoparium, Solidago canadensis, Sorghastrum nutans
# plot_n_sub$dominant_species_code2 <- gsub("^.*_","", plot_n_sub$dominant_species_code)
# 
# dominant_species_list <- plot_n_sub %>%
#   group_by(dominant_species_code2) %>%
#   summarize(n = n())
# 
# plot_sub_lrr9_dom2 <- plot_n_sub %>%
#   filter(dominant_species_code2 == "Andropogon gerardii" | dominant_species_code2 == "Andropogon gerardi" | dominant_species_code2 == "Elymus repens" | dominant_species_code2 == "Poa pratensis" | dominant_species_code2 == "Schizachyrium scoparium" | dominant_species_code2 == "Solidago canadensis" | dominant_species_code2 == "Sorghastrum nutans")
# 
# plot_sub_lrr9_dom2 <- plot_sub_lrr9_dom2 %>%
#   group_by(uniqueid) %>%
#   arrange(uniqueid, year) %>%
#   mutate(prior_year_dom = lag(dominant_relative_abund),
#          prior_year_type = lag(spei9_category)) %>%
#   ungroup()
# normal_biomass9sub_dom2 <- plot_sub_lrr9_dom2  %>%
#   group_by(uniqueid) %>%
#   filter(spei9_category == "Normal" | spei9_category == "Moderate wet" | spei9_category == "Moderate dry") %>%
#   summarize(average_normal_dom = mean(dominant_relative_abund, na.rm=T))
# plot_sub_lrr9_dom2 <- left_join(plot_sub_lrr9_dom2, normal_biomass9sub_dom2)
# 
# lrr9sub_dom2 <- plot_sub_lrr9_dom2 %>%
#   # filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
#   filter(spei9_category == "Extreme wet" | spei9_category == "Extreme dry") %>%
#   mutate(LRR = log(dominant_relative_abund / average_normal_dom))  # Calculate the log response ratio
# 
# lrr.lm9sub_dom2 <- lmer(LRR ~ spei9_category*nitrogen + (1|site/experiment/uniqueid) + (1|year), data = lrr9sub_dom2)
# summary(lrr.lm9sub_dom2)
# simres <- simulateResiduals(lrr.lm9sub_dom2)
# plot(simres)
# 
# # LRR plot for dominance of top seven species
# lrr9sub_dom2 %>%
#   ggplot(aes(x = spei9_category, y = LRR, col = nitrogen)) +
#   stat_summary(fun.data = mean_cl_boot, position = position_dodge(0.2)) + 
#   geom_hline(yintercept=0, linetype='dashed', col = 'black')+
#   theme_bw() +
#   scale_color_manual(values = cols) +
#   labs(x = "Event type", col = "Fertilizer")
# 
# # LRR for evenness
# plot_sub_lrr9_ev <- plot_n_sub %>%
#   group_by(uniqueid) %>%
#   arrange(uniqueid, year) %>%
#   mutate(prior_year_ev = lag(Evar),
#          prior_year_type = lag(spei9_category)) %>%
#   ungroup()
# normal_biomass9sub_ev <- plot_sub_lrr9_ev  %>%
#   group_by(uniqueid) %>%
#   filter(spei9_category == "Normal" | spei9_category == "Moderate wet" | spei9_category == "Moderate dry") %>%
#   summarize(average_normal_ev = mean(Evar, na.rm=T))
# plot_sub_lrr9_ev <- left_join(plot_sub_lrr9_ev, normal_biomass9sub_ev)
# 
# lrr9sub_ev <- plot_sub_lrr9_ev %>%
#   # filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
#   filter(spei9_category == "Extreme wet" | spei9_category == "Extreme dry") %>%
#   mutate(LRR = log(Evar / average_normal_ev))  # Calculate the log response ratio
# 
# lrr.lm9sub_ev <- lmer(LRR ~ spei9_category*nitrogen + (1|site/experiment/uniqueid) + (1|year), data = lrr9sub_ev)
# summary(lrr.lm9sub_ev)
# simres <- simulateResiduals(lrr.lm9sub_ev)
# plot(simres)
# 
# ## Final LRR plot for evenness
# lrr9sub_ev %>%
#   ggplot(aes(x = spei9_category, y = LRR, col = nitrogen)) +
#   stat_summary(fun.data = mean_cl_boot, position = position_dodge(0.2)) + 
#   geom_hline(yintercept=0, linetype='dashed', col = 'black')+
#   theme_bw()
# 
# ## Final LRR plot version 1
# anpp_plot + rich_plot + dom_plot & plot_annotation(tag_levels = 'A')
# 
# # Using SPEI12
# plot_sub_lrr12 <- plot_n_sub %>%
#   group_by(uniqueid) %>%
#   arrange(uniqueid, year) %>%
#   mutate(prior_year_biomass = lag(plot_biomass),
#          prior_year_type = lag(spei12_category)) %>%
#   ungroup()
# normal_biomass12sub <- plot_sub_lrr12  %>%
#   group_by(uniqueid) %>%
#   filter(spei12_category == "Normal" | spei12_category == "Moderate wet" | spei12_category == "Moderate dry") %>%
#   summarize(average_normal_biomass = mean(plot_biomass, na.rm=T))
# plot_sub_lrr12 <- left_join(plot_sub_lrr12, normal_biomass12sub)
# 
# lrr12sub <- plot_sub_lrr12 %>%
#   filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
#   filter(spei12_category == "Extreme wet" | spei12_category == "Extreme dry") %>%
#   mutate(LRR = log(plot_biomass / average_normal_biomass))  # Calculate the log response ratio
# 
# # Take absolute value of SPEI12
# lrr12sub <- lrr12sub %>%
#   mutate(abs_spei12 = abs(spei12))
# 
# lrr.lm12sub <- lmer(LRR ~ abs_spei12*spei12_category*nitrogen + (1|site/experiment/uniqueid) + (1|year), data = lrr12sub)
# summary(lrr.lm12sub)
# Anova(lrr.lm12sub, type = "III")
# simres <- simulateResiduals(lrr.lm12sub)
# plot(simres)
# 
# # Plot the LRR model
# plot_model(
#   lrr.lm12sub,
#   type = "pred",
#   terms= c("abs_spei12", "spei12_category", "nitrogen"),
#   show.data = TRUE
# )
# 
# 
# 
# # LRR control, N, NPK
# plot_npk <- plot_trt %>%
#   filter(nutrients_added == "NPK+" | nutrients_added == "N" | nutrients_added == "NPK" | nutrients_added == "none" | nutrients_added == "no_fertilizer")
# 
# plot_npk %>% group_by(nutrients_added) %>% summarize(n = n())
# 
# # Add column with nitrogen, NPK, or N
# plot_npk <- plot_npk %>%
#   mutate(nitrogen = dplyr::recode(nutrients_added, "NPK+" = "NPK", "none" = "no_fertilizer")) %>%
#   mutate(nitrogen = factor(nitrogen, levels = c("NPK","N","Oceania","no_fertilizer")))
# 
# # Using SPEI3
# plot_npk_lrr3 <- plot_npk %>%
#   group_by(uniqueid) %>%
#   arrange(uniqueid, year) %>%
#   mutate(prior_year_biomass = lag(plot_biomass),
#          prior_year_type = lag(spei3_category)) %>%
#   ungroup()
# normal_biomass3npk <- plot_npk_lrr3  %>%
#   group_by(uniqueid) %>%
#   filter(spei3_category == "Normal" | spei3_category == "Moderate wet" | spei3_category == "Moderate dry") %>%
#   summarize(average_normal_biomass = mean(plot_biomass, na.rm=T))
# plot_npk_lrr3 <- left_join(plot_npk_lrr3, normal_biomass3npk)
# 
# lrr3npk <- plot_npk_lrr3 %>%
#   filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
#   filter(spei3_category == "Extreme wet" | spei3_category == "Extreme dry") %>%
#   mutate(LRR = log(plot_biomass / average_normal_biomass))  # Calculate the log response ratio
# 
# # Take absolute value of SPEI3
# lrr3npk <- lrr3npk %>%
#   mutate(abs_spei3 = abs(spei3))
# 
# lrr.lm3npk <- lmer(LRR ~ abs_spei3*spei3_category*nitrogen + (1|site/experiment/uniqueid) + (1|year), data = lrr3npk)
# summary(lrr.lm3npk)
# Anova(lrr.lm3npk, type = "III")
# simres <- simulateResiduals(lrr.lm3npk)
# plot(simres)
# 
# # Plot the LRR model
# plot_model(
#   lrr.lm3npk,
#   type = "pred",
#   terms= c("abs_spei3", "spei3_category", "nitrogen"),
#   show.data = TRUE
# )
# 
# # Using SPEI6
# plot_npk_lrr6 <- plot_npk %>%
#   group_by(uniqueid) %>%
#   arrange(uniqueid, year) %>%
#   mutate(prior_year_biomass = lag(plot_biomass),
#          prior_year_type = lag(spei6_category)) %>%
#   ungroup()
# normal_biomass6npk <- plot_npk_lrr6  %>%
#   group_by(uniqueid) %>%
#   filter(spei6_category == "Normal" | spei6_category == "Moderate wet" | spei6_category == "Moderate dry") %>%
#   summarize(average_normal_biomass = mean(plot_biomass, na.rm=T))
# plot_npk_lrr6 <- left_join(plot_npk_lrr6, normal_biomass6npk)
# 
# lrr6npk <- plot_npk_lrr6 %>%
#   filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
#   filter(spei6_category == "Extreme wet" | spei6_category == "Extreme dry") %>%
#   mutate(LRR = log(plot_biomass / average_normal_biomass))  # Calculate the log response ratio
# 
# # Take absolute value of SPEI6
# lrr6npk <- lrr6npk %>%
#   mutate(abs_spei6 = abs(spei6))
# 
# lrr.lm6npk <- lmer(LRR ~ abs_spei6*spei6_category*nitrogen + (1|site/experiment/uniqueid) + (1|year), data = lrr6npk)
# summary(lrr.lm6npk)
# Anova(lrr.lm6npk, type = "III")
# simres <- simulateResiduals(lrr.lm6npk)
# plot(simres)
# 
# # Plot the LRR model
# plot_model(
#   lrr.lm6npk,
#   type = "pred",
#   terms= c("abs_spei6", "spei6_category", "nitrogen"),
#   show.data = TRUE
# )
# 
# # Using SPEI9
# plot_npk_lrr9 <- plot_npk %>%
#   group_by(uniqueid) %>%
#   arrange(uniqueid, year) %>%
#   mutate(prior_year_biomass = lag(plot_biomass),
#          prior_year_type = lag(spei9_category)) %>%
#   ungroup()
# normal_biomass9npk <- plot_npk_lrr9  %>%
#   group_by(uniqueid) %>%
#   filter(spei9_category == "Normal" | spei9_category == "Moderate wet" | spei9_category == "Moderate dry") %>%
#   summarize(average_normal_biomass = mean(plot_biomass, na.rm=T))
# plot_npk_lrr9 <- left_join(plot_npk_lrr9, normal_biomass9npk)
# 
# lrr9npk <- plot_npk_lrr9 %>%
#   filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
#   filter(spei9_category == "Extreme wet" | spei9_category == "Extreme dry") %>%
#   mutate(LRR = log(plot_biomass / average_normal_biomass))  # Calculate the log response ratio
# 
# # Take absolute value of SPEI9
# lrr9npk <- lrr9npk %>%
#   mutate(abs_spei9 = abs(spei9))
# 
# lrr.lm9npk <- lmer(LRR ~ abs_spei9*spei9_category*nitrogen + (1|site/experiment/uniqueid) + (1|year), data = lrr9npk)
# summary(lrr.lm9npk)
# Anova(lrr.lm9npk, type = "III")
# simres <- simulateResiduals(lrr.lm9npk)
# plot(simres)
# 
# # Plot the LRR model
# plot_model(
#   lrr.lm9npk,
#   type = "pred",
#   terms= c("abs_spei9", "spei9_category", "nitrogen"),
#   show.data = TRUE
# )
# 
# # Using SPEI12
# plot_npk_lrr12 <- plot_npk %>%
#   group_by(uniqueid) %>%
#   arrange(uniqueid, year) %>%
#   mutate(prior_year_biomass = lag(plot_biomass),
#          prior_year_type = lag(spei12_category)) %>%
#   ungroup()
# normal_biomass12npk <- plot_npk_lrr12  %>%
#   group_by(uniqueid) %>%
#   filter(spei12_category == "Normal" | spei12_category == "Moderate wet" | spei12_category == "Moderate dry") %>%
#   summarize(average_normal_biomass = mean(plot_biomass, na.rm=T))
# plot_npk_lrr12 <- left_join(plot_npk_lrr12, normal_biomass12npk)
# 
# lrr12npk <- plot_npk_lrr12 %>%
#   filter(prior_year_type != "Extreme wet" | prior_year_type != "Extreme dry") %>%
#   filter(spei12_category == "Extreme wet" | spei12_category == "Extreme dry") %>%
#   mutate(LRR = log(plot_biomass / average_normal_biomass))  # Calculate the log response ratio
# 
# # Take absolute value of SPEI12
# lrr12npk <- lrr12npk %>%
#   mutate(abs_spei12 = abs(spei12))
# 
# lrr.lm12npk <- lmer(LRR ~ abs_spei12*spei12_category*nitrogen + (1|site/experiment/uniqueid) + (1|year), data = lrr12npk)
# summary(lrr.lm12npk)
# Anova(lrr.lm12npk, type = "III")
# simres <- simulateResiduals(lrr.lm12npk)
# plot(simres)
# 
# # Plot the LRR model
# plot_model(
#   lrr.lm12npk,
#   type = "pred",
#   terms= c("abs_spei12", "spei12_category", "nitrogen"),
#   show.data = TRUE
# )
# 
# # Resistance
# # Weird NAs... resistance values that we don't have plots for (maybe filtered out plots)
# ece <- distinct(ece)
# ece <-ece %>%
#   rename(year = ex_year)
# ece$year <- as.factor(ece$year)
# plot_ece <- right_join(plot_n_sub, ece)
# plot_ece <- plot_ece %>%
#   filter(!is.na(spei9_category))
# 
# resist.lm9 <- lmer(log(resistance) ~ spei9_category*nitrogen + (1|site/experiment/uniqueid) + (1|year), data = plot_ece)
# summary(resist.lm9)
# emm_resist <- emmeans(resist.lm9, pairwise ~ spei9_category * nitrogen, infer = T)
# summary(emm_resist)
# simres <- simulateResiduals(resist.lm9)
# plot(simres)
# 
# plot_model(
#   resist.lm9,
#   type = "pred",
#   terms= c("spei9_category", "nitrogen"),
#   show.data = F
# ) + theme_bw()
# 
# plot_ece %>%
#   ggplot(aes(x = spei9_category, y = log(resistance), col = spei9_category)) +
#   stat_summary(fun.data = mean_cl_boot, position = position_dodge(0.2), aes(shape = nitrogen)) + 
#   theme_bw() +
#   scale_color_manual(values = cols)
# 
# resist_plot <- plot_ece %>%
#   mutate(nitrogen = relevel(as.factor(nitrogen), 'no_fertilizer', 'nitrogen')) %>%
#   ggplot(aes(x = spei9_category, y = log(resistance), col = spei9_category)) +
#   stat_summary(fun.data = mean_cl_boot, position = position_dodge(0.2), aes(shape = nitrogen)) + 
#   theme_bw() +
#   labs(x = "Event type", y = "ln(resistance)") +
#   theme(legend.position="none") +
#   scale_color_manual(values = cols)
# 
# # Resilience
# resil.lm9 <- lmer(log(resilience) ~ spei9_category*nitrogen + (1|site/experiment/uniqueid) + (1|year), data = plot_ece)
# summary(resil.lm9)
# emm_resil <- emmeans(resil.lm9, pairwise ~ spei9_category * nitrogen, infer = T)
# summary(emm_resil)
# simres <- simulateResiduals(resil.lm9)
# plot(simres)
# 
# plot_model(
#   resil.lm9,
#   type = "pred",
#   terms= c("spei9_category", "nitrogen"),
#   show.data = F
# ) + theme_bw()
# 
# plot_ece %>%
#   ggplot(aes(x = spei9_category, y = log(resilience), col = spei9_category)) +
#   stat_summary(fun.data = mean_cl_boot, position = position_dodge(0.2), aes(shape = nitrogen)) + 
#   theme_bw() +
#   scale_color_manual(values = cols)
# 
# resil_plot <- plot_ece %>%
#   mutate(nitrogen = relevel(as.factor(nitrogen), 'no_fertilizer', 'nitrogen')) %>%
#   ggplot(aes(x = spei9_category, y = log(resilience), col = spei9_category)) +
#   stat_summary(fun.data = mean_cl_boot, position = position_dodge(0.2), aes(shape = nitrogen)) + 
#   theme_bw() +
#   labs(x = "Event type", y = "ln(resilience)") +
#   theme(legend.position="none") +
#   scale_color_manual(values = cols)
# 
# # Final LRR plot version 2
# anpp_plot + resist_plot + resil_plot + rich_plot + dom_plot + guide_area() + plot_layout(guides = 'collect') & plot_annotation(tag_levels = 'A') 


