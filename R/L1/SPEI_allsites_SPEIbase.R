# TITLE:        LTER Grassland Rock: SPEI from SPEIbase for all sites
# AUTHORS:      Caitlin Broderick, Ashley Darst
# COLLABORATORS:  
# DATA INPUT:   Data imported as csv files from shared Google drive L0 folder
# DATA OUTPUT:  (not yet, but generate a list of years /spei/spei categories to match up with plant data)
# PROJECT:      LTER Grassland Rock
# DATE:         April 2024



# SPEIbase: https://spei.csic.es/database.html  <- same as Isbell paper. 

#Previous results suggest that primary productivity
#responds to approximately annual water balances in temperate grasslands

# and isbell paper used different spei values - i am using SPEI 12 as the default
# this was the first they considered - annual water balances
#note i left the lat long in the file names for easy retrievel of different SPEI metric


# Clear all existing data
rm(list=ls())

# Load packaages
library(tidyverse)
library(lubridate)

# Set working directory 
L0_dir <- Sys.getenv("L0DIR")
L1_dir <- Sys.getenv("L1DIR")
list.files(L0_dir)

# Read in datasets
# SPEI-12
kbs <- read.csv(file.path(L0_dir, "KBS/-85.250000_42.250000.csv"), stringsAsFactors = FALSE)
knz <- read.csv(file.path(L0_dir, "KNZ/-96.750000_39.250000.csv"), stringsAsFactors = FALSE)
cdr <- read.csv(file.path(L0_dir, "CDR/-93.250000_45.250000.csv"), stringsAsFactors = FALSE)

#SPEI-9
kbs9 <- read.csv(file.path(L0_dir, "KBS/-85.250000_42.250000_spei9.csv"), stringsAsFactors = FALSE)
knz9 <- read.csv(file.path(L0_dir, "KNZ/-96.750000_39.250000_spei9.csv"), stringsAsFactors = FALSE)
cdr9 <- read.csv(file.path(L0_dir, "CDR/-93.250000_45.250000_spei9.csv"), stringsAsFactors = FALSE)

#SPEI-6
kbs6 <- read.csv(file.path(L0_dir, "KBS/-85.250000_42.250000_spei6.csv"), stringsAsFactors = FALSE)
knz6 <- read.csv(file.path(L0_dir, "KNZ/-96.750000_39.250000_spei6.csv"), stringsAsFactors = FALSE)
cdr6 <- read.csv(file.path(L0_dir, "CDR/-93.250000_45.250000_spei6.csv"), stringsAsFactors = FALSE)

#SPEI-3
kbs3 <- read.csv(file.path(L0_dir, "KBS/-85.250000_42.250000_spei3.csv"), stringsAsFactors = FALSE)
knz3 <- read.csv(file.path(L0_dir, "KNZ/-96.750000_39.250000_spei3.csv"), stringsAsFactors = FALSE)
cdr3 <- read.csv(file.path(L0_dir, "CDR/-93.250000_45.250000_spei3.csv"), stringsAsFactors = FALSE)

# separate columns
# have to do it here! no maatter what I do, opening this in excel goofs it all up
cdr <- cdr %>%  separate_wider_delim(dates.spei12, delim = ";", names = c("date", "spei12"))
kbs <- kbs %>%  separate_wider_delim(dates.spei12, delim = ";", names = c("date", "spei12"))
knz <- knz %>%  separate_wider_delim(dates.spei12, delim = ";", names = c("date", "spei12"))

cdr9 <- cdr9 %>%  separate_wider_delim(dates.spei09, delim = ";", names = c("date", "spei9"))
kbs9 <- kbs9 %>%  separate_wider_delim(dates.spei09, delim = ";", names = c("date", "spei9"))
knz9 <- knz9 %>%  separate_wider_delim(dates.spei09, delim = ";", names = c("date", "spei9"))

cdr6 <- cdr6 %>%  separate_wider_delim(dates.spei06, delim = ";", names = c("date", "spei6"))
kbs6 <- kbs6 %>%  separate_wider_delim(dates.spei06, delim = ";", names = c("date", "spei6"))
knz6 <- knz6 %>%  separate_wider_delim(dates.spei06, delim = ";", names = c("date", "spei6"))

cdr3 <- cdr3 %>%  separate_wider_delim(dates.spei03, delim = ";", names = c("date", "spei3"))
kbs3 <- kbs3 %>%  separate_wider_delim(dates.spei03, delim = ";", names = c("date", "spei3"))
knz3 <- knz3 %>%  separate_wider_delim(dates.spei03, delim = ";", names = c("date", "spei3"))

# merge spei by date
cdr <- merge(cdr, cdr3, by = "date")
cdr <- merge(cdr, cdr6, by = "date")
cdr <- merge(cdr, cdr9, by = "date")

kbs <- merge(kbs, kbs3, by = "date")
kbs <- merge(kbs, kbs6, by = "date")
kbs <- merge(kbs, kbs9, by = "date")

knz <- merge(knz, knz3, by = "date")
knz <- merge(knz, knz6, by = "date")
knz <- merge(knz, knz9, by = "date")

# lubridate the date columns
kbs$date <- lubridate::ymd(as.character(kbs$date))
knz$date <- lubridate::ymd(as.character(knz$date))
cdr$date <- lubridate::ymd(as.character(cdr$date))

kbs$spei12 <- as.numeric (as.character(kbs$spei12))
knz$spei12 <- as.numeric (as.character(knz$spei12))
cdr$spei12 <- as.numeric (as.character(cdr$spei12))

kbs$spei3 <- as.numeric (as.character(kbs$spei3))
knz$spei3 <- as.numeric (as.character(knz$spei3))
cdr$spei3 <- as.numeric (as.character(cdr$spei3))

kbs$spei6 <- as.numeric (as.character(kbs$spei6))
knz$spei6 <- as.numeric (as.character(knz$spei6))
cdr$spei6 <- as.numeric (as.character(cdr$spei6))

kbs$spei9 <- as.numeric (as.character(kbs$spei9))
knz$spei9 <- as.numeric (as.character(knz$spei9))
cdr$spei9 <- as.numeric (as.character(cdr$spei9))

# add year column for graphing purposes
kbs$year <- year(kbs$date)
knz$year <- year(knz$date)
cdr$year <- year(cdr$date)

max(kbs$year)
max(knz$year)
max(cdr$year)
# isbell first filtered for peak biomass harvest at each site. 
  # also, only want 100 years? can add that later 
kbs.aug <- kbs %>% filter(month(date) == 8 & !is.nan(spei12))
knz.aug <- knz %>% filter(month(date) == 8 & !is.nan(spei12))
cdr.aug <- cdr %>% filter(month(date) == 8 & !is.nan(spei12))


# mean? just curious. should be about zero
mean(kbs.aug$spei12) 
mean(knz.aug$spei12) 
mean(cdr.aug$spei12)  # looking good. 





# calculate interquartile range - "normal range"
# this is the diffeerence between upper and lower quartile ranges. 
# Calculate 1 in 10 year event 
# the isbell paper does this a little differently FYI
# UPDATE: SPEI is already inherently doing this. We need to pick values to stay the same across all sites.
# Let's use Isbell's:
#   >= 1.28 = extreme wet
#   < 1.28 - >= 0.67 = moderate dry
#   > -0.67 - < 0.67 = normal
#   > -1.28 - <= -0.67 = moderate dry
#   <= -1.28 = extreme dry

quantile(kbs.aug$spei12, c(0.1, 0.25, 0.75, 0.9))
quantile(knz.aug$spei12, c(0.1, 0.25, 0.75, 0.9))
quantile(cdr.aug$spei12, c(0.1, 0.25, 0.75, 0.9))

# add column for classifying spei3 years based on isbell's categories

kbs.aug <-kbs.aug   %>% mutate(spei3_category = case_when(
  kbs.aug$spei3 <= -1.28 ~ "Extreme dry", 
  kbs.aug$spei3 > -1.28 & kbs.aug$spei3 <= -0.67 ~ "Moderate dry", 
  kbs.aug$spei3 > -0.67 & kbs.aug$spei3 < 0.67 ~ "Normal", 
  kbs.aug$spei3 >= 0.67 & kbs.aug$spei3 < 1.28 ~ "Moderate wet", 
  kbs.aug$spei3 >= 1.28 ~ "Extreme wet"
))


knz.aug <-knz.aug   %>% mutate(spei3_category = case_when(
  knz.aug$spei3 <= -1.28 ~ "Extreme dry", 
  knz.aug$spei3 > -1.28 & knz.aug$spei3 <= -0.67 ~ "Moderate dry", 
  knz.aug$spei3 > -0.67 & knz.aug$spei3 < 0.67 ~ "Normal", 
  knz.aug$spei3 >= 0.67 & knz.aug$spei3 < 1.28 ~ "Moderate wet", 
  knz.aug$spei3 >= 1.28 ~ "Extreme wet"
))


cdr.aug <-cdr.aug   %>% mutate(spei3_category = case_when(
  cdr.aug$spei3 <= -1.28 ~ "Extreme dry", 
  cdr.aug$spei3 > -1.28 & cdr.aug$spei3 <= -0.67 ~ "Moderate dry", 
  cdr.aug$spei3 > -0.67 & cdr.aug$spei3 < 0.67 ~ "Normal", 
  cdr.aug$spei3 >= 0.67 & cdr.aug$spei3 < 1.28 ~ "Moderate wet", 
  cdr.aug$spei3 >= 1.28 ~ "Extreme wet"
))

# add column for classifying spei6 years based on isbell's categories

kbs.aug <-kbs.aug   %>% mutate(spei6_category = case_when(
  kbs.aug$spei6 <= -1.28 ~ "Extreme dry", 
  kbs.aug$spei6 > -1.28 & kbs.aug$spei6 <= -0.67 ~ "Moderate dry", 
  kbs.aug$spei6 > -0.67 & kbs.aug$spei6 < 0.67 ~ "Normal", 
  kbs.aug$spei6 >= 0.67 & kbs.aug$spei6 < 1.28 ~ "Moderate wet", 
  kbs.aug$spei6 >= 1.28 ~ "Extreme wet"
))


knz.aug <-knz.aug   %>% mutate(spei6_category = case_when(
  knz.aug$spei6 <= -1.28 ~ "Extreme dry", 
  knz.aug$spei6 > -1.28 & knz.aug$spei6 <= -0.67 ~ "Moderate dry", 
  knz.aug$spei6 > -0.67 & knz.aug$spei6 < 0.67 ~ "Normal", 
  knz.aug$spei6 >= 0.67 & knz.aug$spei6 < 1.28 ~ "Moderate wet", 
  knz.aug$spei6 >= 1.28 ~ "Extreme wet"
))


cdr.aug <-cdr.aug   %>% mutate(spei6_category = case_when(
  cdr.aug$spei6 <= -1.28 ~ "Extreme dry", 
  cdr.aug$spei6 > -1.28 & cdr.aug$spei6 <= -0.67 ~ "Moderate dry", 
  cdr.aug$spei6 > -0.67 & cdr.aug$spei6 < 0.67 ~ "Normal", 
  cdr.aug$spei6 >= 0.67 & cdr.aug$spei6 < 1.28 ~ "Moderate wet", 
  cdr.aug$spei6 >= 1.28 ~ "Extreme wet"
))

# add column for classifying spei9 years based on isbell's categories

kbs.aug <-kbs.aug   %>% mutate(spei9_category = case_when(
  kbs.aug$spei9 <= -1.28 ~ "Extreme dry", 
  kbs.aug$spei9 > -1.28 & kbs.aug$spei9 <= -0.67 ~ "Moderate dry", 
  kbs.aug$spei9 > -0.67 & kbs.aug$spei9 < 0.67 ~ "Normal", 
  kbs.aug$spei9 >= 0.67 & kbs.aug$spei9 < 1.28 ~ "Moderate wet", 
  kbs.aug$spei9 >= 1.28 ~ "Extreme wet"
))


knz.aug <-knz.aug   %>% mutate(spei9_category = case_when(
  knz.aug$spei9 <= -1.28 ~ "Extreme dry", 
  knz.aug$spei9 > -1.28 & knz.aug$spei9 <= -0.67 ~ "Moderate dry", 
  knz.aug$spei9 > -0.67 & knz.aug$spei9 < 0.67 ~ "Normal", 
  knz.aug$spei9 >= 0.67 & knz.aug$spei9 < 1.28 ~ "Moderate wet", 
  knz.aug$spei9 >= 1.28 ~ "Extreme wet"
))


cdr.aug <-cdr.aug   %>% mutate(spei9_category = case_when(
  cdr.aug$spei9 <= -1.28 ~ "Extreme dry", 
  cdr.aug$spei9 > -1.28 & cdr.aug$spei9 <= -0.67 ~ "Moderate dry", 
  cdr.aug$spei9 > -0.67 & cdr.aug$spei9 < 0.67 ~ "Normal", 
  cdr.aug$spei9 >= 0.67 & cdr.aug$spei9 < 1.28 ~ "Moderate wet", 
  cdr.aug$spei9 >= 1.28 ~ "Extreme wet"
))

# add column for classifying spei12 years based on isbell's categories

kbs.aug <-kbs.aug   %>% mutate(spei12_category = case_when(
  kbs.aug$spei12 <= -1.28 ~ "Extreme dry", 
  kbs.aug$spei12 > -1.28 & kbs.aug$spei12 <= -0.67 ~ "Moderate dry", 
  kbs.aug$spei12 > -0.67 & kbs.aug$spei12 < 0.67 ~ "Normal", 
  kbs.aug$spei12 >= 0.67 & kbs.aug$spei12 < 1.28 ~ "Moderate wet", 
  kbs.aug$spei12 >= 1.28 ~ "Extreme wet"
))


knz.aug <-knz.aug   %>% mutate(spei12_category = case_when(
  knz.aug$spei12 <= -1.28 ~ "Extreme dry", 
  knz.aug$spei12 > -1.28 & knz.aug$spei12 <= -0.67 ~ "Moderate dry", 
  knz.aug$spei12 > -0.67 & knz.aug$spei12 < 0.67 ~ "Normal", 
  knz.aug$spei12 >= 0.67 & knz.aug$spei12 < 1.28 ~ "Moderate wet", 
  knz.aug$spei12 >= 1.28 ~ "Extreme wet"
))


cdr.aug <-cdr.aug   %>% mutate(spei12_category = case_when(
  cdr.aug$spei12 <= -1.28 ~ "Extreme dry", 
  cdr.aug$spei12 > -1.28 & cdr.aug$spei12 <= -0.67 ~ "Moderate dry", 
  cdr.aug$spei12 > -0.67 & cdr.aug$spei12 < 0.67 ~ "Normal", 
  cdr.aug$spei12 >= 0.67 & cdr.aug$spei12 < 1.28 ~ "Moderate wet", 
  cdr.aug$spei12 >= 1.28 ~ "Extreme wet"
))


# order factor
kbs.aug$spei3_category <- factor(kbs.aug$spei3_category, levels = c("Extreme dry", "Moderate dry", "Normal", "Moderate wet", "Extreme wet"))
knz.aug$spei3_category <- factor(knz.aug$spei3_category, levels = c("Extreme dry", "Moderate dry", "Normal", "Moderate wet", "Extreme wet"))
cdr.aug$spei3_category <- factor(cdr.aug$spei3_category, levels = c("Extreme dry", "Moderate dry", "Normal", "Moderate wet", "Extreme wet"))

kbs.aug$spei6_category <- factor(kbs.aug$spei6_category, levels = c("Extreme dry", "Moderate dry", "Normal", "Moderate wet", "Extreme wet"))
knz.aug$spei6_category <- factor(knz.aug$spei6_category, levels = c("Extreme dry", "Moderate dry", "Normal", "Moderate wet", "Extreme wet"))
cdr.aug$spei6_category <- factor(cdr.aug$spei6_category, levels = c("Extreme dry", "Moderate dry", "Normal", "Moderate wet", "Extreme wet"))

kbs.aug$spei9_category <- factor(kbs.aug$spei9_category, levels = c("Extreme dry", "Moderate dry", "Normal", "Moderate wet", "Extreme wet"))
knz.aug$spei9_category <- factor(knz.aug$spei9_category, levels = c("Extreme dry", "Moderate dry", "Normal", "Moderate wet", "Extreme wet"))
cdr.aug$spei9_category <- factor(cdr.aug$spei9_category, levels = c("Extreme dry", "Moderate dry", "Normal", "Moderate wet", "Extreme wet"))

kbs.aug$spei12_category <- factor(kbs.aug$spei12_category, levels = c("Extreme dry", "Moderate dry", "Normal", "Moderate wet", "Extreme wet"))
knz.aug$spei12_category <- factor(knz.aug$spei12_category, levels = c("Extreme dry", "Moderate dry", "Normal", "Moderate wet", "Extreme wet"))
cdr.aug$spei12_category <- factor(cdr.aug$spei12_category, levels = c("Extreme dry", "Moderate dry", "Normal", "Moderate wet", "Extreme wet"))

# Plot SPEI6
ggplot(kbs.aug , aes (x = year, y = spei6)) +
  geom_line(alpha = 0.5) + 
  geom_point( aes(color = spei6_category), size = 1.25) + theme_bw() +
  scale_x_continuous(breaks = seq(1900, 2025, by = 10)) +
  annotate( "text", label = "KBS",
            x = 1910, y = 2, size = 5, colour = "black") + 
  scale_color_manual(values = c("#F5191CFF", "#E78200FF", "#EAC728FF", "#81BB95FF", "#3B99B1FF"))

ggplot(knz.aug , aes (x = year, y = spei6)) +
  geom_line(alpha = 0.5) + 
  geom_point( aes(color = spei6_category), size = 1.25) + theme_bw() +
  scale_x_continuous(breaks = seq(1900, 2025, by = 10)) +
  annotate( "text", label = "KNZ",
            x = 1910, y = 2, size = 5, colour = "black") + 
  scale_color_manual(values = c("#F5191CFF", "#E78200FF", "#EAC728FF", "#81BB95FF", "#3B99B1FF"))

ggplot(cdr.aug , aes (x = year, y = spei6)) +
  geom_line(alpha = 0.5) + 
  geom_point( aes(color = spei6_category), size = 1.25) + theme_bw() +
  scale_x_continuous(breaks = seq(1900, 2025, by = 10)) +
  annotate( "text", label = "CDR",
            x = 1910, y = 2, size = 5, colour = "black") + 
  scale_color_manual(values = c("#F5191CFF", "#E78200FF", "#EAC728FF", "#81BB95FF", "#3B99B1FF"))

# Plot SPEI12
kbs_spei <- ggplot(kbs.aug , aes (x = year, y = spei12)) +
  geom_line(alpha = 0.5) + 
  geom_point( aes(color = spei12_category), size = 1.25) + theme_bw() +
  scale_x_continuous(breaks = seq(1900, 2025, by = 10)) +
  annotate( "text", label = "KBS",
    x = 1910, y = 2, size = 5, colour = "black") + 
  scale_color_manual(values = c("#F5191CFF", "#E78200FF", "#EAC728FF", "#81BB95FF", "#3B99B1FF"))
kbs_spei
ggsave("KBS_SPEI.png",plot =kbs_spei, dpi = 300, width =8, height = 4, units = "in")


knz_spei <- ggplot(knz.aug , aes (x = year, y = spei12)) +
  geom_line(alpha = 0.5) + 
  geom_point( aes(color = spei12_category), size = 1.25) + theme_bw() +
  scale_x_continuous(breaks = seq(1900, 2025, by = 10)) +
  annotate( "text", label = "KNZ",
            x = 1910, y = 2, size = 5, colour = "black") + 
  scale_color_manual(values = c("#F5191CFF", "#E78200FF", "#EAC728FF", "#81BB95FF", "#3B99B1FF"))
knz_spei
ggsave("KNZ_SPEI.png",plot =knz_spei, dpi = 300, width =8, height = 4, units = "in")


cdr_spei <- ggplot(cdr.aug , aes (x = year, y = spei12)) +
  geom_line(alpha = 0.5) + 
  geom_point( aes(color = spei12_category), size = 1.25) + theme_bw() +
  scale_x_continuous(breaks = seq(1900, 2025, by = 10)) +
  annotate( "text", label = "CDR",
            x = 1910, y = 2, size = 5, colour = "black") + 
  scale_color_manual(values = c("#F5191CFF", "#E78200FF", "#EAC728FF", "#81BB95FF", "#3B99B1FF"))
cdr_spei
ggsave("CDR_SPEI.png",plot =cdr_spei, dpi = 300, width =8, height = 4, units = "in")




#combine all sites

cdr.aug$site = "CDR"
kbs.aug$site = "KBS"
knz.aug$site = "KNZ"

allsite.aug <- rbind (cdr.aug, kbs.aug, knz.aug)
head(allsite.aug)
allsite.aug

allsite.aug <- allsite.aug %>% select (-date) # remove date column, do not think we need it. just will join based on year and site.

head(allsite.aug )



write.csv(allsite.aug, file.path(L1_dir, "./SPEI_12_allsites.csv"), row.names=F)
