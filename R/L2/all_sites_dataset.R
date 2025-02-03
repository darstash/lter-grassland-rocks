# TITLE:        LTER Grassland Rock: Merge CDR, KBS, and KNZ datasets
# AUTHORS:      Ashley Darst
# COLLABORATORS:  
# DATA INPUT:   Data imported as cleaned site-specific csv files from shared Google drive L1 folder
# DATA OUTPUT:    
# PROJECT:      LTER Grassland Rock
# DATE:         July 2024

# this code is to harmonize plant biomass and plant composition datasets from CDR, KBS, and KNZ

# Clear all existing data
rm(list=ls())

# Load packages
library(tidyverse)
library(gtools)
library(lubridate)
library(TNRS)

# Set working directory 
L0_dir <- Sys.getenv("L0DIR")
L1_dir <- Sys.getenv("L1DIR")
L2_dir <- Sys.getenv("L2DIR")
list.files(L1_dir)

# Read in CSV files ----
# CDR
CDR_metadata <- read.csv(file.path(L1_dir, "CDR_metadata.csv"))
CDR_species_abundance <- read.csv(file.path(L1_dir, "CDR_specieslevel_abundance.csv"))
CDR_plot_metrics <- read.csv(file.path(L1_dir, "CDR_plotlevel_metrics.csv"))

# KBS
KBS_metadata <- read.csv(file.path(L1_dir, "KBS_metadata.csv"))
KBS_species_abundance <- read.csv(file.path(L1_dir, "KBS_species_level_abundance.csv"))
KBS_plot_metrics <- read.csv(file.path(L1_dir, "KBS_plot_level_metrics.csv"))

# KNZ 
KNZ_metadata <- read.csv(file.path(L1_dir, "KNZ_metadata.csv"))
KNZ_species_abundance <- read.csv(file.path(L1_dir, "KNZ_specieslevel_abundance.csv"))
KNZ_plot_metrics <- read.csv(file.path(L1_dir, "KNZ_plotlevel_metrics.csv"))

# Merge datasets ----
# Plot metrics

# Make plot metric datasets compatible for merge
CDR_plot_metrics$measurement_scale_biomass <- gsub("\\m.*", "", CDR_plot_metrics$measurement_scale_biomass)
CDR_plot_metrics$measurement_scale_biomass <- as.double(CDR_plot_metrics$measurement_scale_biomass)
CDR_plot_metrics$measurement_scale_cover <- gsub("\\m.*", "", CDR_plot_metrics$measurement_scale_cover)
CDR_plot_metrics$measurement_scale_cover <- as.double(CDR_plot_metrics$measurement_scale_cover)
CDR_plot_metrics$plot <- as.character(CDR_plot_metrics$plot)
KNZ_plot_metrics$plot <- as.character(KNZ_plot_metrics$plot)

plot_metrics <- full_join(CDR_plot_metrics, KBS_plot_metrics)
plot_metrics <- full_join(plot_metrics, KNZ_plot_metrics)

plot_metrics <- plot_metrics %>%
  select(-c(source, shannon, evenness, richness, X, dataset))

write.csv(plot_metrics, file.path(L2_dir, "./plot_metrics.csv"), row.names=F)

# Metadata

# Make metadata datasets compatible for merge
CDR_metadata$plot <- as.character(CDR_metadata$plot)
KNZ_metadata$plot <- as.character(KNZ_metadata$plot)

metadata <- full_join(CDR_metadata, KBS_metadata)
metadata <- full_join(metadata, KNZ_metadata) # Two different nutrient added columns for KNZ

# Merge temperature columns
metadata <- metadata %>%
  mutate(meantemp = coalesce(meantemp, temperature)) %>%
  select(-temperature)
metadata <- metadata %>%
  mutate(annualprecip = coalesce(annualprecip, precipitation)) %>%
  select(-c(precipitation, growtemp, experiment))
metadata <- metadata %>%
  mutate(growprecip = coalesce(growprecip, growing_precipitation)) %>%
  select(-c(growing_precipitation, X))
  
write.csv(metadata, file.path(L2_dir, "./metadata.csv"), row.names=F)

# Species abundance

# Make species abundance datasets compatible for merge
CDR_species_abundance$plot <- as.character(CDR_species_abundance$plot)
KNZ_species_abundance$plot <- as.character(KNZ_species_abundance$plot)

species_abundance <- full_join(CDR_species_abundance, KBS_species_abundance)
species_abundance <- full_join(species_abundance, KNZ_species_abundance)

species_abundance <- species_abundance %>%
  select(-X)

######### UNFINISHED SPECIES NAMES FIXING #####################
# make species information dataset
# idea: at this stage, don't touch names in the dataset anymore, but have a 
# dataframe that assigns each weirdo and non-weirdo entry in the species column
# a cleaned up version of names 

## Try TNRS

# TO DO HERE: loop TNRS confirmed species names (in results) back into the 
# species_abundance data set. Create lists of things that are non-plants and
# different identification levels. Filter out unwanted things from dataset 
# (probably non-plant things and spp things)

# species_abundance <- species_abundance %>%
#   mutate(
#     species_raw = species,
#     species = species_raw %>%
#       # fix some common typo patterns in the data
#       str_trim() %>%
#       str_squish() %>%
#       str_to_sentence() %>%
#       gsub(pattern = " \\(\\*\\)", replacement = "") %>%
#       gsub(pattern = " \\(l\\.\\)", replacement = "_") %>%
#       gsub(pattern = " l\\.", replacement = "_") %>%
#       gsub(pattern = "Unk ", replacement = "Unknown ") %>%
#       gsub(pattern = "Unk_", replacement = "Unknown ") %>%
#       gsub(pattern = "_$", replacement = "",),
#     
#     # fix names that are not resolved by TNRS
#     species = case_when(species %in% "Acalyp virgi"                   ~ "Acalypha virginica",
#                         species %in% "Achillea millefolium(lanulosa)" ~ "Achillea millefolium lanulosa",
#                         species %in% "Agerat altis"                   ~ "Ageratina altissima",
#                         species %in% "Agrost hyema"                   ~ "Agrostis hyemalis",
#                         species %in% "Ambros artem"                   ~ "Ambrosia artemisiifolia",
#                         species %in% "Ambros psilo"                   ~ "Ambrosia psilostachya",
#                         species %in% "Ambrosia art"                   ~ "Ambrosia artemisiifolia",
#                         species %in% "Ambrosia psilo"                 ~ "Ambrosia psilostachya",
#                         species %in% "Ammann cocci"                   ~ "Ammannia coccinea",              # KNZ, is that correct?
#                         species %in% "Amorph canes"                   ~ "Amorpha canescens",
#                         species %in% "Amorph fruti"                   ~ "Amorpha fruticosa",              # KNZ, is that correct?
#                         species %in% "Anaphalis margaritacea_ Benth. & Hook. f." ~ "Anaphalis margaritacea",
#                         species %in% "Androp gerar"                   ~ "Andropogon gerardii",
#                         species %in% "Antenn negle"                   ~ "Antennaria neglecta",
#                         species %in% "Antennaria plantaginifolia_ Richards." ~ "Antennaria plantaginifolia",  
#                         species %in% "Apocyn canna"                   ~ "Apocynum cannabinum",
#                         species %in% "Arabidopsis thaliana_ Heynh."   ~ "Arabidopsis thaliana",
#                         species %in% "Arabis glabra_ Bernh."          ~ "Arabis glabra",
#                         species %in% "Arrhenatherum elatius_ Beauv. ex j. & C. Presl" ~ "Arrhenatherum elatius",
#                         species %in% "Artemi ludov"                   ~ "Artemesia ludoviciana",
#                         species %in% "Artemisia (caudata) campestris" ~ "Artemisia campestris",
#                         species %in% "Asclep sulli"                   ~ "Artemisia ludoviciana",
#                         species %in% "Asclep syria"                   ~ "Asclepias syriaca",
#                         species %in% "Asclep tuber"                   ~ "Asclepias tuberosa",
#                         species %in% "Asclep verti"                   ~ "Asclepias verticillata",
#                         species %in% "Asclep virdf"                   ~ "Asclepias viridiflora",
#                         species %in% "Asclep virds"                   ~ "Asclepias viridis",
#                         species %in% "Asplenium platyneuron_ Oakes"   ~ "Asplenium platyneuron",
#                         species %in% "Aster basal leaves"             ~ "Aster sp.",
#                         species %in% "Astrag crass"                   ~ "Astragalus crassicarpus",
#                         species %in% "Bothri bladh"                   ~ "Bothri bladh",                    # KNZ ?!
#                         species %in% "Boutel curti"                   ~ "Bouteloua curtipendula",
#                         species %in% "Boutel dacty"                   ~ "Bouteloua dactyloides",
#                         species %in% "Boutel graci"                   ~ "Bouteloua gracilis",
#                         species %in% "Boutel hirsu"                   ~ "Bouteloua hirsuta",
#                         species %in% "Bricke eupat"                   ~ "Brickellia eupatorioides",
#                         species %in% "Brickellia (kuhnia) eupatoriodes" ~ "Brickellia eupatorioides",
#                         species %in% "Capsella bursa-pastoris_ Medicus" ~ " Capsella bursa-pastoris",
#                         species %in% "Ceanot herba"                   ~ "Ceanothus herbaceus",
#                         species %in% "Celtis occid"                   ~ "Celtis occidentalis",
#                         species %in% "Centaurea stoebe_ ssp. micranthos (gugler) hayek" ~ "Centaurea stoebe micranthos",
#                         species %in% "Cercis canad"                   ~ "Cercis canadensis",               # KNZ, is that correct?
#                         species %in% "Chamae fasci"                   ~ "Chamaecrista fasciculata",        # KNZ, is that correct?
#                         species %in% "Cirsiu altis"                   ~ "Cirsium altissimum",
#                         species %in% "Cirsiu undul"                   ~ "Cirsium undulatum",
#                         species %in% "Cirsium altissimum_ Spreng."    ~ "Cirsium altissimum",
#                         species %in% "Cirsium arvense_ Scop."         ~ "Cirsium arvense",
#                         species %in% "Conyza canad"                   ~ "Conyza canadensis",
#                         species %in% "Conyza canadensis_ Cronq."      ~ "Conyza canadensis",
#                         species %in% "Cornus alternifolia_f."         ~ "Cornus alternifolia",
#                         species %in% "Cornus drumm"                   ~ "Cornus drummondii",
#                         species %in% "Crepis capillaris_ Wallr."      ~ "Crepis capillaris",
#                         species %in% "Croton monan"                   ~ "Croton monanthogynus",
#                         species %in% "Cuscut glome"                   ~ "Cuscuta glomerata",
#                         species %in% "Cuscuta"                        ~ "Cuscuta sp.",
#                         species %in% "Cynanc laeve"                   ~ "Cynanchum laeve",                 # KNZ, is that correct?
#                         species %in% "Cyperu lupul"                   ~ "Cyperus lupulinus",
#                         species %in% "Cyperus_spp"                    ~ "Cyperus spp.",
#                         species %in% "Cyperus_spp."                   ~ "Cyperus spp.",
#                         species %in% "Dalea multi"                    ~ "Dalea multiflora",
#                         species %in% "Dalea purpu"                    ~ "Dalea purpurea",
#                         species %in% "Desman illin"                   ~ "Desmanthus illinoensis",
#                         species %in% "Desmod illin"                   ~ "Desmodium illinoense",
#                         species %in% "Desmod panic"                   ~ "Desmodium paniculatum",
#                         species %in% "Desmodium canescens_ Dc."       ~ "Desmodium canescens",
#                         species %in% "Desmodium marilandicum_ Dc."    ~ "Desmodium marilandicum",
#                         species %in% "Desmodium paniculatum_ Dc."     ~ "Desmodium paniculatum",
#                         species %in% "Dichan oligo"                   ~ "Dichanthelium oligosanthes",
#                         species %in% "Dichan ovale"                   ~ "Dichanthelium ovale",
#                         species %in% "Digita cogna"                   ~ "Digitaria cognata",
#                         species %in% "Digitaria sanguinalis_ Scop."   ~ "Digitaria sanguinalis",
#                         species %in% "Echina angus"                   ~ "Echinacea angustifolia",
#                         species %in% "Echinacea purpurea_ Moench"     ~ "Echinacea purpurea",
#                         species %in% "Echinochloa crus-galli_ Beauv." ~ "Echinochloa crus-galli",
#                         species %in% "Elymus canad"                   ~ "Elymus canadensis",
#                         species %in% "Elymus repens_ Gould"           ~ "Elymus repens",
#                         species %in% "Elymus virgi"                   ~ "Elymus virginicus",
#                         species %in% "Eragro spect"                   ~ "Eragrostis spectabilis",
#                         species %in% "Eriger phila"                   ~ "Erigeron philadelphicus",
#                         species %in% "Eriger strig"                   ~ "Erigeron strigosus",
#                         species %in% "Erigeron"                       ~ "Erigeron sp.",
#                         species %in% "Erigeron annuus_ Pers."         ~ "Erigeron annuus",
#                         species %in% "Eupato altis"                   ~ "Eupatorium altissimum",
#                         species %in% "Eupatorium"                     ~ "Eupatorium sp.",
#                         species %in% "Eupatorium eupertorioides"      ~ "Eupatorium eupertorioides",     # KNZ?!
#                         species %in% "Euphor cyath"                   ~ "Euphorbia cyathophora",         # KNZ, is that correct?
#                         species %in% "Euphor denta"                   ~ "Euphorbia dentata",
#                         species %in% "Euphor margi"                   ~ "Euphorbia marginata",
#                         species %in% "Euphor nutan"                   ~ "Euphorbia nutans",
#                         species %in% "Euphor serpe"                   ~ "Euphorbia serpens",
#                         species %in% "Euphor spath"                   ~ "Euphorbia spathulata",
#                         species %in% "Euphorbia (supina) maculata"    ~ "Euphorbia maculata",
#                         species %in% "Euthamia graminifolia_ Nutt."   ~ "Euthamia graminifolia",
#                         species %in% "Geum canad"                     ~ "Geum canadense",
#                         species %in% "Gledit triac"                   ~ "Gleditsia triacanthos",
#                         species %in% "Helian annuu"                   ~ "Helianthus annuus",
#                         species %in% "Helian pauci"                   ~ "Helianthus pauciflorus",
#                         species %in% "Heliopsis helianthoides_ Sweet" ~ "Heliopsis helianthoides",
#                         species %in% "Kummer stipu"                   ~ "Kummer stipu",                  # KNZ?!
#                         species %in% "Lactuc serri"                   ~ "Lactuca serriola",
#                         species %in% "Lactuca serriola scariola"      ~ "Lactuca serriola (scariola)",
#                         species %in% "Lepidi densi"                   ~ "Lepidium densiflorum",
#                         species %in% "Lepidium campestre_ R.br."      ~ "Lepidium campestre",
#                         species %in% "Lesped capit"                   ~ "Lespedeza capitata",
#                         species %in% "Lesped viola"                   ~ "Lespedeza violacea",
#                         species %in% "Leucos multi"                   ~ "Leucos multi",                  # KNZ?!
#                         species %in% "Liatri punct"                   ~ "Liatris punctata",
#                         species %in% "Lithos incis"                   ~ "Lithospermum incisum",
#                         species %in% "Lychnis latifolia ssp. Alba"    ~ "Lychnis alba",
#                         species %in% "Lycopu ameri"                   ~ "Lycopus americanus",            # KNZ, is that correct?
#                         species %in% "Melilo offic"                   ~ "Melilotus officinalis",
#                         species %in% "Melilotus officinalis_ Lam."    ~ "Melilotus officinalis",
#                         species %in% "Mimosa (schrankii) quadrivavlis"~ "Mimosa quadrivalvis",           # KNZ, is that correct?
#                         species %in% "Mimosa quadr"                   ~ "Mimosa quadrivalvis",
#                         species %in% "Mirabi albid"                   ~ "Mirabilis albida",
#                         species %in% "Mirabi linea"                   ~ "Mirabilis linearis",
#                         species %in% "Monard fistu"                   ~ "Monarda fistulosa",
#                         species %in% "Muhlen racem"                   ~ "Muhlenbergia racemosa",
#                         species %in% "Muhlenbergia species"           ~ "Muhlenbergia sp.",
#                         species %in% "Oenoth curti"                   ~ "Oenothera curtiflora",
#                         species %in% "Packer platt"                   ~ "Packera plattensis",
#                         species %in% "Panicu capil"                   ~ "Panicum capillare",
#                         species %in% "Panicu virga"                   ~ "Panicum virgatum",
#                         species %in% "Pariet pensy"                   ~ "Parietaria pensylvanica",
#                         species %in% "Parthe quinq"                   ~ "Parthenocissus quinquefolia",
#                         species %in% "Parthenocissus quinquefolia_ Planch." ~ "Parthenocissus quinquefolia",
#                         species %in% "Pediom escul"                   ~ "Pediomelum esculentum",
#                         species %in% "Physal heter"                   ~ "Physalis heterophylla",
#                         species %in% "Physal virgi"                   ~ "Physalis virginiana",
#                         species %in% "Physalis"                       ~ "Physalis sp.",
#                         species %in% "Planta virgi"                   ~ "Plantago virginica",
#                         species %in% "Plantago (purshii) patagonica"  ~ "Plantago purshii patagonica",
#                         species %in% "Poa prate"                      ~ "Poa pratensis",
#                         species %in% "Polyga verti"                   ~ "Polygala verticillata",
#                         species %in% "Populu delto"                   ~ "Populus deltoides",
#                         species %in% "Prunus ameri"                   ~ "Prunus americana",
#                         species %in% "Quercus borealis-ellipsoidalis" ~ "Quercus sp.",
#                         species %in% "Ratibi colum"                   ~ "Ratibida columnifera",
#                         species %in% "Rubus occid"                    ~ "Rubus occidentalis",
#                         species %in% "Rudbeckia (hirta) serotina"     ~ "Rudbeckia hirta",
#                         species %in% "Rudbeckia hirta_ var. pulcherrima farw." ~ "Rudbeckia hirta var. pulcherrima",
#                         species %in% "Ruellia"                        ~ "Ruellia sp.",
#                         species %in% "Rumex altis"                    ~ "Rumex altissimus",            # KNZ, is that correct?
#                         species %in% "Rumex species"                  ~ "Rumex sp.",
#                         species %in% "Schiza scopa"                   ~ "Schizachyrium scoparium",
#                         species %in% "Scutel parvu"                   ~ "Scutellaria parvula",         # KNZ, is that correct?
#                         species %in% "Senna maril"                    ~ "Senna marilandica",           # KNZ, is that correct?
#                         species %in% "Setaria viridis_ Beauv."        ~ "Setaria viridi",
#                         species %in% "Silene antir"                   ~ "Silene antirrhina",
#                         species %in% "Silphi integ"                   ~ "Silphium integrifolium",      # KNZ, is that correct?
#                         species %in% "Silphi lacin"                   ~ "Silphium laciniatum",
#                         species %in% "Sisymbrium officinale_ Scop."   ~ "Sisymbrium officinale",
#                         species %in% "Sisyri campe"                   ~ "Sisyrinchium campestre",
#                         species %in% "Solanu ptyca"                   ~ "Solanum ptycanthum",
#                         species %in% "Solanu rostr"                   ~ "Solanum rostratum",
#                         species %in% "Solida misso"                   ~ "Solidago missouriensis",
#                         species %in% "Sonchus asper_ Hill"            ~ "Sonchus asper",
#                         species %in% "Sorghastrum nutans_ Nash ex small" ~ "Sorghastrum nutans",
#                         species %in% "Spiran verna"                   ~ "Spiranthes vernalis",
#                         species %in% "Sporob compo"                   ~ "Sporobolus compositus",
#                         species %in% "Sporob crypt"                   ~ "Sporobolus cryptandrus",
#                         species %in% "Sporobolus (compositus) asper"  ~ "Sporobolus compositus",
#                         species %in% "Stellaria media_ Vill."         ~ "Stellaria media",
#                         species %in% "Sympho orbic"                   ~ "Symphocarpus orbiculatus",
#                         species %in% "Symphy spp."                    ~ "Symphyotrichum ssp.",
#                         species %in% "Symphy drumm"                   ~ "Symphyotrichum drummondii",   # KNZ, is that correct?
#                         species %in% "Symphy erico"                   ~ "Symphyotrichum ericoides",
#                         species %in% "Symphy oblon"                   ~ "Symphyotrichum_oblongifolium",
#                         species %in% "Symphyotrichym (aster) ericoides" ~ "Symphyotrichum oblongifolium",
#                         species %in% "Symphyotrichym (aster) oblongifolium" ~ "Symphyotrichum oblongifolium",
#                         species %in% "Taraxa offic"                   ~ "Taraxacum officinale",
#                         species %in% "Toxico nutta"                   ~ "Toxicoscordion nuttallii",
#                         species %in% "Toxico radic"                   ~ "Toxicodendron radicans",
#                         species %in% "Toxicodendron radicans_ Ktze."  ~ "Toxicodendron radicans",
#                         species %in% "Tragopogon dubius (major)"      ~ "Tragopogon dubius major",
#                         species %in% "Tridens flavus_ A.s.hitchc."    ~ "	Tridens flavus",
#                         species %in% "Trifolium"                      ~ "Trifolium sp.",
#                         species %in% "Tripsa dacty"                   ~ "Tripsacum dactyloides",
#                         species %in% "Ulmus ameri"                    ~ "Ulmus americana",
#                         species %in% "Velvet leaf"                    ~ "Abutilon theophrasti",
#                         species %in% "Verbes alter"                   ~ "Verbes alter",                # KNZ?!
#                         species %in% "Vernon baldw"                   ~ "Vernonia baldwinii",
#                         species %in% "Viola nephr"                    ~ "Viola nephrophylla",
#                         species %in% "Viola petidifida"               ~ "Viola pedatifida",            # KNZ, is that correct?
#                         species %in% "Vitis sp. (kbs.us)"             ~ "Vitis sp.",
#                         species %in% "Zantho ameri"                   ~ "Zanthoxylum americanum",
#                         species %in% "Zigadenus"                      ~ "Zigadenus sp.",
#                         
#                         # unify things that are identified at a lower taxonomic level than species or genus
#                         species %in% c("1st year woody", 
#                                        "Miscellaneous woody",
#                                        "Unknown_woody",
#                                        "Woody")                      ~ "Unknown woody sp.",
#                         species %in% c("Miscellaneous woody plants") ~ "Unknown woody spp.",
#                         species %in% c("Annual forb", 
#                                        "Another unknown dicot", 
#                                        "Forb", 
#                                        "Forb sp.",  
#                                        "Miscellaneous forb",
#                                        "Miscellaneous forb2",  
#                                        "Psoralea?",   
#                                        "Unknown _comandra-like",
#                                        "Unknown annual forb",
#                                        "Unknown dicot",
#                                        "Unknown forb",
#                                        "Unknown forb 1",
#                                        "Unknown forb 3",
#                                        "Unknown_forb",
#                                        "Unknown_forb seedling",
#                                        "Composite basal leaves")    ~ "Unknown forb sp.",
#                         species %in% c("Miscellaneous herbs",
#                                        "Forb seedlings",
#                                        "Dicots",
#                                        "Forb seedlings",
#                                        "Forbes",
#                                        "Forbs",
#                                        "Miscellaneous forbs",
#                                        "Miscellaneous herbs")       ~ "Unknown forb spp.",
#                         species %in% c("Miscellaneous grasses 2",
#                                        "Miscellaneous grass",
#                                        "Unknown grass",
#                                        "Unknown grass 1",
#                                        "Unknown sedge")             ~ "Unknown grassXXX sp.", # X needed to not get wrong match with TNRS
#                         species %in% c("Miscellaneous grasses",
#                                        "Miscellaneous rushes",
#                                        "Carex spp.",
#                                        "C3 grasses",
#                                        "C4 grasses",
#                                        "Miscellaneous grasses",
#                                        "Miscellaneous rushes",
#                                        "Miscellaneous sedges",
#                                        "Sedges")                    ~ "Unknown grassXXX spp.",
#                         species %in% c("Miscellaneous sp.",
#                                        "Miscellaneous sp.",
#                                        "Miscellaneous sp. 2",
#                                        "Unknown",
#                                        "Unknown 1",
#                                        "Unknown 2",
#                                        "Unknown 3",
#                                        "Unknown sp.",
#                                        "Unknown monocot")           ~ "Unknown sp.",
#                         species %in% c("Monocots",
#                                        "Unsorted")                  ~ "Unknown spp.",
#                         species %in% c("Unknown fabaceae",
#                                        "Unkown fabaceae",
#                                        "Unknown legume")            ~ "Fabaceae sp.",
#                         species %in% c("Miscellaneous legumes",
#                                        "Legumes")                   ~ "Fabaceae spp.",
#                         species %in% "Unknown asteraceae"           ~ "Asteraceae sp.",
#                         species %in% "Unknown brassicaceae"         ~ "Brassicaceae sp.",
#                         species %in% "Unknown cupressaceae sp."     ~ "Cupressaceae sp.",
#                         species %in% "Unknown elm"                  ~ "Ulmus sp.",
#                         species %in% "Unknown lamiaceae"            ~ "Lamiaceae sp.",
#                         species %in% "Unknown orchidaceae"          ~ "Orchidaceae sp.",
#                         species %in% "Unknown rosacae"              ~ "Rosaceae sp.",
#                         species %in% "Unknown solanaceae"           ~ "Solanaceae sp.",
#                         
#                         # unify non-living plant things names
#                         species %in% c("Bare_ground", 
#                                        "Ground")                    ~ "Bare ground",
#                         species %in% c("Litter", 
#                                        "Miscellaneous woody litter",
#                                        "Other litter",
#                                        "Woody litter")              ~ "LitterX",            # X needed to not get wrong match with TNRS
#                         species %in% "Other" ~"Something else",
#                         .default = as.character(species)),
#     species = gsub(species, pattern = "_", replacement = " ")
#   )
#                                  
# 
# 
# species_list <- species_abundance$species %>% unique()
# results <- TNRS(taxonomic_names = species_list)
# 
# 
# results %>% 
#   filter(!Name_matched_rank %in% c("species", "variety", "subspecies")) %>%
#   select(Name_submitted, Name_matched, Name_matched_rank, Genus_submitted, Genus_matched, Accepted_family, Unmatched_terms, WarningsEng) %>%
#   arrange(Name_submitted) 













## CDR list of things
genus_sp_in_biomass <- c("Acer sp.",
                         "Allium sp.",
                         "Alnus sp.",
                         "Apocynum sp.",
                         "Arabis sp.",
                         "Aristida sp.",
                         "Asclepias sp.",
                         "Aster sp.",
                         "Bromus sp.",
                         "Calamovilfa sp.",
                         "Carex sp.",
                         "Chenopodium sp.",
                         "Cirsium sp.",
                         "Cyperus sp.",
                         "Digitaria sp.",
                         "Digitaria sp.",
                         "Equisetum sp.",
                         "Eragrostis sp.",
                         "Erigeron sp.",
                         "Galium sp.",
                         "Helianthus sp.",
                         "Hieracium sp.",
                         "Juncus sp.",
                         "Lactuca sp.",
                         "Liatris sp.",
                         "Lithospermum sp.",
                         "Melilotus sp.",
                         "Oenothera sp.",
                         "Oxalis sp.",
                         "Panicum sp.",
                         "Parthenocissus sp.",
                         "Penstemon sp.",
                         "Pinus sp.",
                         "Plantago sp.",
                         "Poa sp.",
                         "Polygala sp.",
                         "Polygonatum sp.",
                         "Potentilla sp.",
                         "Prunus sp.",
                         "Quercus borealis-ellipsoidalis",
                         "Quercus sp.",
                         "Ranunculus sp.",
                         "Rhus sp.",
                         "Rudbeckia sp.",
                         "Rubus sp.",
                         "Rumex sp.",
                         "Salix sp.",
                         "Sedges",
                         "Senecio sp.",
                         "Setaria sp.",
                         "Silene sp.",
                         "Solidago sp.",
                         "Sporobolus sp.",
                         "Tradescantia sp.",
                         "Tragopogon sp.",
                         "Trifolium sp.",
                         "Viola sp.")

non_plant_things_in_biomass <- c("Corn litter", 
                                 "Fungi",
                                 "Fungi sp.",
                                 "Ground",
                                 "Miscellaneous litter",
                                 "Mosses",
                                 "Mosses & lichens",
                                 "Mosses & lichens 2",
                                 "moses & lichens",
                                 "Mosses and lichens",
                                 "Lichen",
                                 "Lichens",
                                 "Other",
                                 "Other animal diggings",
                                 "Other litter",
                                 "Pine litter",
                                 "Pine cones",
                                 "pine needles",
                                 "Pine needles",
                                 "Pine twigs",
                                 "woody debris",
                                 "Woody debris",
                                 "Leaves")

maybe_plant_things_in_biomass <- c("Grass seedlings",
                                   "1st year woody",
                                   "Bryophyte",
                                   "C3 grasses",
                                   "C4 grasses",
                                   "Forb",
                                   "Forb seedlings",
                                   "Forb sp.",
                                   "Forbes",
                                   "Legumes",
                                   "Miscellaneous forb",
                                   "Miscellaneous forbs",
                                   "Miscellaneous forb 1",
                                   "Miscellaneous forb 2",
                                   "Miscellaneous grass",
                                   "Miscellaneous grasses",
                                   "Miscellaneous grasses 2",
                                   "Miscellaneous herb",
                                   "Miscellaneous herbs",
                                   "Miscellaneous herbs 2",
                                   "Miscellaneous legumes",
                                   "Miscellaneous liter",
                                   "Miscellaneous litter",
                                   "Miscellaneous rushes",
                                   "Miscellaneous sedges",
                                   "Miscellaneous sp.",
                                   "Miscellaneous sp. 2",
                                   "Miscellaneous woody tree",
                                   "Miscellaneous  woody",
                                   "Miscellaneous Woody",
                                   "Miscellaneous woody 1",
                                   "Miscellaneous woody 2",
                                   "Miscellaneous woody plants",
                                   "Miscellaneous woody plants 1",
                                   "Miscellaneous woody plants 2",
                                   "Miscellaneous woody litter",
                                   "Unknown",
                                   "Unknown cupressaceae sp.",
                                   "Unknown fabaceae",
                                   "Unknown lamiaceae",
                                   "Unknown sp.",
                                   "Sedges",
                                   "Woody")



species_list <- cbind.data.frame(
  species = species_abundance$species %>%
    factor() %>%
    levels(),
  species_clean = species_abundance$species %>%
    factor() %>%
    levels() %>%
    str_trim() %>%
    str_squish() %>%
    str_to_sentence() %>%
    gsub(pattern = " \\(\\*\\)", replacement = "") %>%
    gsub(pattern = " \\(l\\.\\)", replacement = "_") %>%
    gsub(pattern = " l\\.", replacement = "_") %>%
    gsub(pattern = "Unk ", replacement = "Unknown ") %>%
    gsub(pattern = "Unk_", replacement = "Unknown ") %>%
    gsub(pattern = " ", replacement = "_") 
) %>%
  mutate(
    category = case_when(startsWith(as.character(species), "Unknown")       ~ "unidentified_plant_things",
                         startsWith(as.character(species), "Miscellaneous") ~ "unidentified_plant_things",
                         species_clean %in% genus_sp_in_biomass           | species %in% genus_sp_in_biomass ~"identified_to_genus",
                         species_clean %in% non_plant_things_in_biomass   | species %in% non_plant_things_in_biomass ~ "non_plant_things",
                         species_clean %in% maybe_plant_things_in_biomass | species %in% maybe_plant_things_in_biomass  ~ "unidentified_plant_things")
  ) %>%
  arrange(species_clean)

species_list %>% select(species_clean, category) %>% unique()

species_abundance %>% filter(species %in% (species_list %>% filter(str_detect(species, " \\(\\*\\)")))$species) %>% select(site, higher_order_organization) %>% unique()

write.csv(species_abundance, file.path(L2_dir, "./species_abundance.csv"), row.names=F)

