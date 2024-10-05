library(tidyverse)

L2_dir <- Sys.getenv("L2DIR")

list.files(L2_dir)

species_abundance_SPEI <- read.csv(file.path(L2_dir, "species_abundance_SPEI.csv"), stringsAsFactors = F)
plot_metrics_SPEI <- read.csv(file.path(L2_dir, "plot_metrics_SPEI.csv"), stringsAsFactors = F)

species_abundance_SPEI_Metric <- species_abundance_SPEI %>% 
  group_by(year, site, dataset, plot, higher_order_organization, uniqueid,
           cover_method, area_sampled_bio, area_sampled_cover,
           spei12, spei3, spei6, spei9, spei12_category) %>% 
  summarize(Berger_Parker = max(relative_abundance),
            Richness = n())

plot_metrics_SPEI_diversity <- plot_metrics_SPEI %>% 
  left_join(., species_abundance_SPEI_Metric, by = c("year", "site",  "higher_order_organization", "plot", "uniqueid",
                                                     "spei12", "spei3", "spei6", "spei9", "spei12_category"))


write.csv(species_abundance_SPEI_Metric, file.path(L2_dir, "./plot_metrics_SPEI_diversity.csv"), row.names=F)
