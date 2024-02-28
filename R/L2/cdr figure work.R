# TITLE:         LTER Grassland Rock: CDR - exploring diversity-stability 
# AUTHOR:        Max Zaret
# COLLABORATORS: LTER synthesis group
# DATA INPUT: 
cdr_data #pulled this data fom CDR_explore_L1.R script
# DATA OUTPUT:
# PROJECT:       LTER Grassland Rock synthesis group
# DATE:          2/28/2024

cdr_data

library(ggpubr) #for quick visual stats on plots
# super basic questions - 

# across all expts and trts does diversity increase ANPP?
ggplot(cdr_data, aes (x = richness, y = plot_biomass)) +  # richness
  geom_point() +
  geom_smooth(method="lm") + 
  stat_cor() +
  xlab("Richness") + ylab("Biomass (g/m2)") +
  theme_classic () 

################################
# look at 2012 drought
################################
cdr_data_dr <- cdr_data %>% 
  filter(year > 2010 & year < 2014) # subset 2011, 2012, 2013 data (before, during, after)
head(cdr_data_dr)
table(cdr_data_dr$year)



# convert to wide format 
cdr_data_dr_wide <- cdr_data_dr %>% 
  ungroup() %>%
  select (c (year, uniqueid, plot_biomass)) %>%  # select columns I will use
  pivot_wider(names_from= year, values_from = plot_biomass) %>%  # pivot wider so each year is own column
  mutate (perc_change_dr = ((`2012`- `2011`) / `2011`) * 100) %>% # add column for perc change during drought
  mutate (perc_change_recov = ((`2013`- `2011`) / `2011`) * 100) # add column for perc recovery back to 2011


cdr_rich_2011 <- cdr_data %>% 
  filter(year == 2011) # just get richness for 2011

# merger 2011 richness for each plot with the wide-format biomass data fror 
cdr_data_dr_wide <- merge(cdr_data_dr_wide, cdr_rich_2011 , by = "uniqueid", all=TRUE)


# drought resistance
ggplot(cdr_data_dr_wide, aes (x = richness, y = perc_change_dr)) +  
  geom_point()+
  #scale_color_manual(values = c("blue", "red","green")) + 
  geom_hline(yintercept=0, col = "red") + 
  geom_smooth(method = "lm", se=FALSE) + 
  stat_cor() +
  xlab("rich") + ylab("perc change with drought")
# can add facet by treatment here! such as nutrients_added
# does N addition change sensitivity to dr

# drought recovery
ggplot(cdr_data_dr_wide, aes (x = richness, y = perc_change_recov )) + 
  geom_point()+
  scale_color_manual(values = c("blue", "red","green")) + 
  geom_hline(yintercept=0, col = "red") + 
  geom_smooth(method = "lm") + 
  stat_cor() +
  xlab("rich") + ylab("perc recovery") 

# convert to wide format 
cdr_data_dr_wide <- cdr_data_dr %>% 
  ungroup() %>%
  select (c (year, uniqueid, plot_biomass)) %>%  # select columns I will use
  pivot_wider(names_from= year, values_from = plot_biomass) %>%  # pivot wider so each year is own column
  mutate (resistance = `2011`/abs(`2012` - `2011`)) %>% # resistance based on Isbell 2015
  mutate (resilience = abs((`2012` - `2011`)/(`2013`- `2011`))) #resilience based on Isbell 2015 PNAS



cdr_rich_2011 <- cdr_data %>% 
  filter(year == 2011) # just get richness for 2011

# merger 2011 richness for each plot with the wide-format biomass data fror 
cdr_data_dr_wide <- merge(cdr_data_dr_wide, cdr_rich_2011 , by = "uniqueid", all=TRUE)


# drought resistance
ggplot(cdr_data_dr_wide, aes (x = richness, y = log(resistance) )) +  
  geom_point()+
  #scale_color_manual(values = c("blue", "red","green")) + 
  geom_hline(yintercept=0, col = "red") + 
  geom_smooth(method = "lm") + 
  stat_cor() +
  xlab("rich") + ylab("log scale resistance") + 
  labs(title="2012 Drought")
# can add facet by treatment here! such as nutrients_added
# does N addition change sensitivity to dr

# drought resilience
ggplot(cdr_data_dr_wide, aes (x = richness, y = log(resilience) )) + 
  geom_point()+
  #scale_color_manual(values = c("blue", "red","green")) + 
  geom_hline(yintercept=0, col = "red") + 
  geom_smooth(method = "lm") + 
  stat_cor() +
  xlab("rich") + ylab("log resilience") + 
  labs(title="2012 Drought")


#1988 drought?????
cdr_data_dr <- cdr_data %>% 
  filter(year > 1986 & year < 1990) # subset 2011, 2012, 2013 data (before, during, after)
head(cdr_data_dr)
table(cdr_data_dr$year)

# convert to wide format 
cdr_data_dr_wide <- cdr_data_dr %>% 
  ungroup() %>%
  select (c (year, uniqueid, plot_biomass)) %>%  # select columns I will use
  pivot_wider(names_from= year, values_from = plot_biomass) %>%  # pivot wider so each year is own column
  mutate (resistance = `1987`/abs(`1988` - `1987`)) %>% # resistance based on Isbell 2015
  mutate (resilience = abs((`1988` - `1987`)/(`1989`- `1987`))) #resilience based on Isbell 2015 PNAS



cdr_rich_2011 <- cdr_data %>% 
  filter(year == 1987) # just get richness for 2011

# merger 2011 richness for each plot with the wide-format biomass data fror 
cdr_data_dr_wide <- merge(cdr_data_dr_wide, cdr_rich_2011 , by = "uniqueid", all=TRUE)


# drought resistance
ggplot(cdr_data_dr_wide, aes (x = richness, y = log(resistance) )) +  
  geom_point()+
  #scale_color_manual(values = c("blue", "red","green")) + 
  geom_hline(yintercept=0, col = "red") + 
  geom_smooth(method = "lm") + 
  stat_cor() +
  xlab("rich") + ylab("log scale resistance") + 
  labs(title="1988 Drought")

# can add facet by treatment here! such as nutrients_added
# does N addition change sensitivity to dr

# drought resilience
ggplot(cdr_data_dr_wide, aes (x = richness, y = log(resilience) )) + 
  geom_point()+
  #scale_color_manual(values = c("blue", "red","green")) + 
  geom_hline(yintercept=0, col = "red") + 
  geom_smooth(method = "lm") + 
  stat_cor() +
  xlab("rich") + ylab("log resilience") + 
  labs(title="1988 Drought")

#try averaging normal years and comparing to climate events then compare results

cdr_normal <- cdr_data %>%
  filter(year != 2012 | year != 1988) %>%
  group_by(uniqueid) %>%
  summarize(norm_biomass = mean(plot_biomass),
            norm_richness = mean(richness))

cdr_dr <- cdr_data %>%
  ungroup() %>%
  filter(year == 2012) %>%
  mutate(drought_biomass = plot_biomass) %>%
  select(uniqueid, drought_biomass)

cdr_recovery <- cdr_data %>%
  ungroup() %>%
  filter(year == 2013) %>%
  mutate(recovery_biomass = plot_biomass) %>%
  select(uniqueid, recovery_biomass)

cdr_compare <- merge(cdr_normal, cdr_dr, by="uniqueid", all=TRUE)
cdr_compare <- merge(cdr_compare, cdr_recovery, by="uniqueid", all=TRUE)

cdr_compare <- cdr_compare %>%
  mutate (resistance = norm_biomass/abs(drought_biomass - norm_biomass)) %>% # resistance based on Isbell 2015
  mutate (resilience = abs((drought_biomass - norm_biomass)/(recovery_biomass- norm_biomass)))

# drought resistance
ggplot(cdr_compare, aes (x = norm_richness, y = log(resistance) )) +  
  geom_point()+
  #scale_color_manual(values = c("blue", "red","green")) + 
  geom_hline(yintercept=0, col = "red") + 
  geom_smooth(method = "lm") + 
  stat_cor() +
  xlab("rich") + ylab("log scale resistance") + 
  #scale_y_continuous(limits=c(0,6))
  labs(title="2012 Drought with 'normal year' biomass metric")

# can add facet by treatment here! such as nutrients_added
# does N addition change sensitivity to dr

# drought resilience
ggplot(cdr_compare, aes (x = norm_richness, y = log(resilience))) + 
  geom_point()+
  #scale_color_manual(values = c("blue", "red","green")) + 
  geom_hline(yintercept=0, col = "red") + 
  geom_smooth(method = "lm") + 
  stat_cor() +
  xlab("rich") + ylab("log resilience") + 
  labs(title="2012 Drought with 'normal year' biomass metric")
