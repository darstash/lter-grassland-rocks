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

species_abundance2 <- species_abundance %>%
  mutate(
    species_raw = species,
    species = species_raw %>%
      # fix some common typo patterns in the data
      str_trim() %>%
      str_squish() %>%
      str_to_sentence() %>%
      gsub(pattern = " \\(\\*\\)", replacement = "") %>%
      gsub(pattern = " \\(l\\.\\)", replacement = "_") %>%
      gsub(pattern = " l\\.", replacement = "_") %>%
      gsub(pattern = "Unk ", replacement = "Unknown ") %>%
      gsub(pattern = "Unk_", replacement = "Unknown ") %>%
      gsub(pattern = "_$", replacement = "",),

    # fix names that are not resolved by TNRS
    species = case_when(species %in% "Acalyp virgi"                   ~ "Acalypha virginica",
                        species %in% "Achill mille"                   ~ "Achillea millefolium",
                        species %in% "Achillea millefolium(lanulosa)" ~ "Achillea millefolium lanulosa",
                        species %in% "Agerat altis"                   ~ "Ageratina altissima",
                        species %in% "Agrost hyema"                   ~ "Agrostis hyemalis",
                        species %in% "Ambros artem"                   ~ "Ambrosia artemisiifolia",
                        species %in% "Ambros psilo"                   ~ "Ambrosia psilostachya",
                        species %in% "Ambrosia art"                   ~ "Ambrosia artemisiifolia",
                        species %in% "Ambrosia psilo"                 ~ "Ambrosia psilostachya",
                        species %in% "Ammann cocci"                   ~ "Ammannia coccinea",              # KNZ, is that correct?
                        species %in% "Amorph canes"                   ~ "Amorpha canescens",
                        species %in% "Amorph fruti"                   ~ "Amorpha fruticosa",              # KNZ, is that correct?
                        species %in% "Amphicarpa bracteata"           ~ "Amphicarpaea bracteata",
                        species %in% "Anaphalis margaritacea_ Benth. & Hook. f." ~ "Anaphalis margaritacea",
                        species %in% "Androp gerar"                   ~ "Andropogon gerardii",
                        species %in% "Antenn negle"                   ~ "Antennaria neglecta",
                        species %in% "Antennaria plantaginifolia_ Richards." ~ "Antennaria plantaginifolia",
                        species %in% "Apocyn canna"                   ~ "Apocynum cannabinum",
                        species %in% "Arabidopsis thaliana_ Heynh."   ~ "Arabidopsis thaliana",
                        species %in% "Arabis glabra_ Bernh."          ~ "Arabis glabra",
                        species %in% "Arrhenatherum elatius_ Beauv. ex j. & C. Presl" ~ "Arrhenatherum elatius",
                        species %in% "Artemi ludov"                   ~ "Artemesia ludoviciana",
                        species %in% "Artemisia (caudata) campestris" ~ "Artemisia campestris",
                        species %in% "Asclep sulli"                   ~ "Asclepias sullivantii",
                        species %in% "Asclep syria"                   ~ "Asclepias syriaca",
                        species %in% "Asclep tuber"                   ~ "Asclepias tuberosa",
                        species %in% "Asclep verti"                   ~ "Asclepias verticillata",
                        species %in% "Asclep virdf"                   ~ "Asclepias viridiflora",
                        species %in% "Asclep virds"                   ~ "Asclepias viridis",
                        species %in% "Asplenium platyneuron_ Oakes"   ~ "Asplenium platyneuron",
                        species %in% "Aster basal leaves"             ~ "Aster sp.",
                        species %in% "Astrag crass"                   ~ "Astragalus crassicarpus",
                        species %in% "Bothri bladh"                   ~ "Bothriochloa bladhii",
                        species %in% "Boutel curti"                   ~ "Bouteloua curtipendula",
                        species %in% "Boutel dacty"                   ~ "Bouteloua dactyloides",
                        species %in% "Boutel graci"                   ~ "Bouteloua gracilis",
                        species %in% "Boutel hirsu"                   ~ "Bouteloua hirsuta",
                        species %in% "Bricke eupat"                   ~ "Brickellia eupatorioides",
                        species %in% "Brickellia (kuhnia) eupatoriodes" ~ "Brickellia eupatorioides",
                        species %in% "Calylophus serrulatus"          ~ "Oenothera serrulata",
                        species %in% "Capsella bursa-pastoris_ Medicus" ~ "Capsella bursa-pastoris",
                        species %in% "Ceanot herba"                   ~ "Ceanothus herbaceus",
                        species %in% "Celtis occid"                   ~ "Celtis occidentalis",
                        species %in% "Centaurea stoebe_ ssp. micranthos (gugler) hayek" ~ "Centaurea stoebe micranthos",
                        species %in% "Cercis canad"                   ~ "Cercis canadensis",
                        species %in% "Chamae fasci"                   ~ "Chamaecrista fasciculata",
                        species %in% "Cirsiu altis"                   ~ "Cirsium altissimum",
                        species %in% "Cirsiu undul"                   ~ "Cirsium undulatum",
                        species %in% "Cirsium altissimum_ Spreng."    ~ "Cirsium altissimum",
                        species %in% "Cirsium arvense_ Scop."         ~ "Cirsium arvense",
                        species %in% "Conyza canad"                   ~ "Conyza canadensis",
                        species %in% "Conyza canadensis_ Cronq."      ~ "Conyza canadensis",
                        species %in% "Cornus alternifolia_f."         ~ "Cornus alternifolia",
                        species %in% "Cornus drumm"                   ~ "Cornus drummondii",
                        species %in% "Crepis capillaris_ Wallr."      ~ "Crepis capillaris",
                        species %in% "Croton monan"                   ~ "Croton monanthogynus",
                        species %in% "Cuscut glome"                   ~ "Cuscuta glomerata",
                        species %in% "Cuscuta"                        ~ "Cuscuta sp.",
                        species %in% "Cynanc laeve"                   ~ "Cynanchum laeve",
                        species %in% "Cyperu lupul"                   ~ "Cyperus lupulinus",
                        species %in% "Cyperus_spp"                    ~ "Cyperus spp.",
                        species %in% "Cyperus_spp."                   ~ "Cyperus spp.",
                        species %in% "Dalea multi"                    ~ "Dalea multiflora",
                        species %in% "Dalea purpu"                    ~ "Dalea purpurea",
                        species %in% "Desman illin"                   ~ "Desmanthus illinoensis",
                        species %in% "Desmod illin"                   ~ "Desmodium illinoense",
                        species %in% "Desmod panic"                   ~ "Desmodium paniculatum",
                        species %in% "Desmodium canescens_ Dc."       ~ "Desmodium canescens",
                        species %in% "Desmodium marilandicum_ Dc."    ~ "Desmodium marilandicum",
                        species %in% "Desmodium paniculatum_ Dc."     ~ "Desmodium paniculatum",
                        species %in% "Dichan oligo"                   ~ "Dichanthelium oligosanthes",
                        species %in% "Dichan ovale"                   ~ "Dichanthelium ovale",
                        species %in% "Digita cogna"                   ~ "Digitaria cognata",
                        species %in% "Digitaria sanguinalis_ Scop."   ~ "Digitaria sanguinalis",
                        species %in% "Echina angus"                   ~ "Echinacea angustifolia",
                        species %in% "Echinacea purpurea_ Moench"     ~ "Echinacea purpurea",
                        species %in% "Echinochloa crus-galli_ Beauv." ~ "Echinochloa crus-galli",
                        species %in% "Elymus canad"                   ~ "Elymus canadensis",
                        species %in% "Elymus repens_ Gould"           ~ "Elymus repens",
                        species %in% "Elymus virgi"                   ~ "Elymus virginicus",
                        species %in% "Eragro spect"                   ~ "Eragrostis spectabilis",
                        species %in% "Eriger phila"                   ~ "Erigeron philadelphicus",
                        species %in% "Eriger strig"                   ~ "Erigeron strigosus",
                        species %in% "Erigeron"                       ~ "Erigeron sp.",
                        species %in% "Erigeron annuus_ Pers."         ~ "Erigeron annuus",
                        species %in% "Eupato altis"                   ~ "Eupatorium altissimum",
                        species %in% "Eupatorium"                     ~ "Eupatorium sp.",
                        species %in% "Eupatorium eupertorioides"      ~ "Asteraceae sp.",     # KNZ?! Don't think so. there is Eupatorium altissimum and Brickellia eupatorioidies
                        species %in% "Euphor cyath"                   ~ "Euphorbia cyathophora",
                        species %in% "Euphor denta"                   ~ "Euphorbia dentata",
                        species %in% "Euphor margi"                   ~ "Euphorbia marginata",
                        species %in% "Euphor nutan"                   ~ "Euphorbia nutans",
                        species %in% "Euphor serpe"                   ~ "Euphorbia serpens",
                        species %in% "Euphor spath"                   ~ "Euphorbia spathulata",
                        species %in% "Euphorbia (supina) maculata"    ~ "Euphorbia maculata",
                        species %in% "Euthamia graminifolia_ Nutt."   ~ "Euthamia graminifolia",
                        species %in% "Geum canad"                     ~ "Geum canadense",
                        species %in% "Gledit triac"                   ~ "Gleditsia triacanthos",
                        species %in% "Helian annuu"                   ~ "Helianthus annuus",
                        species %in% "Helian pauci"                   ~ "Helianthus pauciflorus",
                        species %in% "Heliopsis helianthoides_ Sweet" ~ "Heliopsis helianthoides",
                        species %in% "Koeler macra"                   ~ "Koeleria macrantha",
                        species %in% "Kummer stipu"                   ~ "Kummerowia stipulacea",
                        species %in% "Lactuc serri"                   ~ "Lactuca serriola",
                        species %in% "Lactuca serriola scariola"      ~ "Lactuca serriola (scariola)",
                        species %in% "Lepidi densi"                   ~ "Lepidium densiflorum",
                        species %in% "Lepidium campestre_ R.br."      ~ "Lepidium campestre",
                        species %in% "Lesped capit"                   ~ "Lespedeza capitata",
                        species %in% "Lesped viola"                   ~ "Lespedeza violacea",
                        species %in% "Leucos multi"                   ~ "Leucospora multifida",
                        species %in% "Liatri punct"                   ~ "Liatris punctata",
                        species %in% "Linum sulca"                    ~ "Linum sulcatum",
                        species %in% "Lithos incis"                   ~ "Lithospermum incisum",
                        species %in% "Lychnis latifolia ssp. Alba"    ~ "Lychnis alba",
                        species %in% "Lycopu ameri"                   ~ "Lycopus americanus",
                        species %in% "Melilo offic"                   ~ "Melilotus officinalis",
                        species %in% "Melilotus officinalis_ Lam."    ~ "Melilotus officinalis",
                        species %in% "Mimosa (schrankii) quadrivavlis"~ "Mimosa quadrivalvis",
                        species %in% "Mimosa quadr"                   ~ "Mimosa quadrivalvis",
                        species %in% "Mirabi albid"                   ~ "Mirabilis albida",
                        species %in% "Mirabi linea"                   ~ "Mirabilis linearis",
                        species %in% "Monard fistu"                   ~ "Monarda fistulosa",
                        species %in% "Muhlen racem"                   ~ "Muhlenbergia racemosa",
                        species %in% "Muhlenbergia species"           ~ "Muhlenbergia sp.",
                        species %in% "Oenoth curti"                   ~ "Oenothera curtiflora",
                        species %in% "Packer platt"                   ~ "Packera plattensis",
                        species %in% "Panicu capil"                   ~ "Panicum capillare",
                        species %in% "Panicu virga"                   ~ "Panicum virgatum",
                        species %in% "Pariet pensy"                   ~ "Parietaria pensylvanica",
                        species %in% "Parthe quinq"                   ~ "Parthenocissus quinquefolia",
                        species %in% "Parthenocissus quinquefolia_ Planch." ~ "Parthenocissus quinquefolia",
                        species %in% "Pediom escul"                   ~ "Pediomelum esculentum",
                        species %in% "Petalostemum purpureum"         ~ "Dalea purpurea",
                        species %in% "Petalostemum sp."               ~ "Dalea sp.",
                        species %in% "Physal heter"                   ~ "Physalis heterophylla",
                        species %in% "Physal longi"                   ~ "Physalis longifolia",
                        species %in% "Physal virgi"                   ~ "Physalis virginiana",
                        species %in% "Physalis"                       ~ "Physalis sp.",
                        species %in% "Planta virgi"                   ~ "Plantago virginica",
                        species %in% "Plantago (purshii) patagonica"  ~ "Plantago purshii patagonica",
                        species %in% "Poa prate"                      ~ "Poa pratensis",
                        species %in% "Polyga verti"                   ~ "Polygala verticillata",
                        species %in% "Populu delto"                   ~ "Populus deltoides",
                        species %in% "Prunus ameri"                   ~ "Prunus americana",
                        species %in% "Quercus borealis-ellipsoidalis" ~ "Quercus sp.",
                        species %in% "Ratibi colum"                   ~ "Ratibida columnifera",
                        species %in% "Rhus armoa"                     ~ "Rhus aromatica",
                        species %in% "Rubus occid"                    ~ "Rubus occidentalis",
                        species %in% "Rubus pensi"                    ~ "Rubus pensilvanicus",
                        species %in% "Rudbeckia (hirta) serotina"     ~ "Rudbeckia hirta",
                        species %in% "Rudbeckia hirta_ var. pulcherrima farw." ~ "Rudbeckia hirta var. pulcherrima",
                        species %in% "Ruellia"                        ~ "Ruellia sp.",
                        species %in% "Rumex altis"                    ~ "Rumex altissimus",
                        species %in% "Rumex species"                  ~ "Rumex sp.",
                        species %in% "Schiza scopa"                   ~ "Schizachyrium scoparium",
                        species %in% "Scutel parvu"                   ~ "Scutellaria parvula",
                        species %in% "Senna maril"                    ~ "Senna marilandica",
                        species %in% "Setaria viridis_ Beauv."        ~ "Setaria viridis",
                        species %in% "Silene antir"                   ~ "Silene antirrhina",
                        species %in% "Silphi integ"                   ~ "Silphium integrifolium",
                        species %in% "Silphi lacin"                   ~ "Silphium laciniatum",
                        species %in% "Sisymbrium officinale_ Scop."   ~ "Sisymbrium officinale",
                        species %in% "Sisyri campe"                   ~ "Sisyrinchium campestre",
                        species %in% "Solanu ptyca"                   ~ "Solanum ptycanthum",
                        species %in% "Solanu carol"                   ~ "Solanum carolinense",
                        species %in% "Solanu rostr"                   ~ "Solanum rostratum",
                        species %in% "Solida altis"                   ~ "Solidago altissima",
                        species %in% "Solida misso"                   ~ "Solidago missouriensis",
                        species %in% "Sonchus asper_ Hill"            ~ "Sonchus asper",
                        species %in% "Sorgha nutan"                   ~ "Sorghastrum nutans",
                        species %in% "Sorghastrum nutans_ Nash ex small" ~ "Sorghastrum nutans",
                        species %in% "Spiran verna"                   ~ "Spiranthes vernalis",
                        species %in% "Sporob compo"                   ~ "Sporobolus compositus",
                        species %in% "Sporob crypt"                   ~ "Sporobolus cryptandrus",
                        species %in% "Sporobolus (compositus) asper"  ~ "Sporobolus compositus",
                        species %in% "Stellaria media_ Vill."         ~ "Stellaria media",
                        species %in% "Sympho orbic"                   ~ "Symphoricarpos orbiculatus",
                        species %in% "Symphy spp."                    ~ "Symphyotrichum ssp.",
                        species %in% "Symphy drumm"                   ~ "Symphyotrichum drummondii",
                        species %in% "Symphy erico"                   ~ "Symphyotrichum ericoides",
                        species %in% "Symphy oblon"                   ~ "Symphyotrichum_oblongifolium",
                        species %in% "Symphyotrichym (aster) ericoides" ~ "Symphyotrichum ericoides",
                        species %in% "Symphyotrichym (aster) oblongifolium" ~ "Symphyotrichum oblongifolium",
                        species %in% "Taraxa offic"                   ~ "Taraxacum officinale",
                        species %in% "Toxico nutta"                   ~ "Toxicoscordion nuttallii",
                        species %in% "Toxico radic"                   ~ "Toxicodendron radicans",
                        species %in% "Toxicodendron radicans_ Ktze."  ~ "Toxicodendron radicans",
                        species %in% "Tragopogon dubius (major)"      ~ "Tragopogon dubius major",
                        species %in% "Tridens flavus_ A.s.hitchc."    ~ "	Tridens flavus",
                        species %in% "Trifolium"                      ~ "Trifolium sp.",
                        species %in% "Tripsa dacty"                   ~ "Tripsacum dactyloides",
                        species %in% "Ulmus ameri"                    ~ "Ulmus americana",
                        species %in% "Velvet leaf"                    ~ "Abutilon theophrasti",
                        species %in% "Verbes alter"                   ~ "Verbesina alternifolia",
                        species %in% "Vernon baldw"                   ~ "Vernonia baldwinii",
                        species %in% "Viola nephr"                    ~ "Viola nephrophylla",
                        species %in% "Viola petidifida"               ~ "Viola pedatifida",
                        species %in% "Viola_rafinesquii"              ~ "Viola rafinesquei",
                        species %in% "Viola rafinesquii"              ~ "Viola rafinesquei", 
                        species %in% "Vitis sp. (kbs.us)"             ~ "Vitis sp.",
                        species %in% "Zantho ameri"                   ~ "Zanthoxylum americanum",
                        species %in% "Zigadenus"                      ~ "Zigadenus sp.",

                        # unify things that are identified at a lower taxonomic level than species or genus
                        species %in% c("1st year woody",
                                       "Miscellaneous woody",
                                       "Unknown_woody",
                                       "Woody")                      ~ "Unknown woodyXXX sp.",
                        species %in% c("Miscellaneous woody plants") ~ "Unknown woodyXXX spp.", # X needed to not get wrong match with TNRS
                        species %in% c("Annual forb",
                                       "Another unknown dicot",
                                       "Forb",
                                       "Forb sp.",
                                       "Miscellaneous forb",
                                       "Miscellaneous forb2",
                                       "Psoralea?",
                                       "Unknown _comandra-like",
                                       "Unknown annual forb",
                                       "Unknown dicot",
                                       "Unknown forb",
                                       "Unknown forb 1",
                                       "Unknown forb 3",
                                       "Unknown_forb",
                                       "Unknown_forb seedling",
                                       "Composite basal leaves")    ~ "Unknown forb sp.",
                        species %in% c("Miscellaneous herbs",
                                       "Forb seedlings",
                                       "Dicots",
                                       "Forb seedlings",
                                       "Forbes",
                                       "Forbs",
                                       "Miscellaneous forbs",
                                       "Miscellaneous herbs")       ~ "Unknown forb spp.",
                        species %in% c("Miscellaneous grasses 2",
                                       "Miscellaneous grass",
                                       "Unknown grass",
                                       "Unknown grass 1",
                                       "Unknown sedge")             ~ "Unknown grassXXX sp.", # X needed to not get wrong match with TNRS
                        species %in% c("Miscellaneous grasses",
                                       "Miscellaneous rushes",
                                       "Carex spp.",
                                       "C3 grasses",
                                       "C4 grasses",
                                       "Miscellaneous grasses",
                                       "Miscellaneous rushes",
                                       "Miscellaneous sedges",
                                       "Sedges")                    ~ "Unknown grassXXX spp.",
                        species %in% c("Miscellaneous sp.",
                                       "Miscellaneous sp.",
                                       "Miscellaneous sp. 2",
                                       "Unknown",
                                       "Unknown 1",
                                       "Unknown 2",
                                       "Unknown 3",
                                       "Unknown sp.",
                                       "Unknown monocot")           ~ "Unknown sp.",
                        species %in% c("Monocots",
                                       "Unsorted")                  ~ "Unknown spp.",
                        species %in% c("Unknown fabaceae",
                                       "Unkown fabaceae",
                                       "Unknown legume")            ~ "Fabaceae sp.",
                        species %in% c("Miscellaneous legumes",
                                       "Legumes")                   ~ "Fabaceae spp.",
                        species %in% "Unknown asteraceae"           ~ "Asteraceae sp.",
                        species %in% "Unknown brassicaceae"         ~ "Brassicaceae sp.",
                        species %in% "Unknown cupressaceae sp."     ~ "Cupressaceae sp.",
                        species %in% "Unknown elm"                  ~ "Ulmus sp.",
                        species %in% "Unknown lamiaceae"            ~ "Lamiaceae sp.",
                        species %in% "Unknown orchidaceae"          ~ "Orchidaceae sp.",
                        species %in% "Unknown rosacae"              ~ "Rosaceae sp.",
                        species %in% "Unknown solanaceae"           ~ "Solanaceae sp.",
                        species %in% "Unknown pteridophyta"         ~ "Pteridophyta sp.",

                        # unify non-living plant things names
                        species %in% c("Bare_ground",
                                       "Ground")                    ~ "Bare ground",
                        species %in% c("Litter",
                                       "Miscellaneous woody litter",
                                       "Other litter",
                                       "Woody litter")              ~ "LitterXXX",          # X needed to not get wrong match with TNRS
                        species %in% "Other" ~"Something else",
                        species %in% "Fungi sp."                    ~ "FungiXXX",           # X needed to not get wrong match with TNRS
                        species %in% "Other animal diggings"        ~ "Animal diggings",
                        .default = as.character(species)),
    species = gsub(species, pattern = "_", replacement = " ")
  )



species_list <- species_abundance2$species %>% unique()
results <- TNRS(taxonomic_names = species_list)

species_name_match_list <- species_abundance2 %>%
  merge(results %>% select(Name_submitted, Accepted_name, Accepted_species, Accepted_name_rank),
        by.x = "species", by.y = "Name_submitted") %>%
  mutate(species_final = ifelse(Accepted_name %in% "", species,
                                ifelse(Accepted_name_rank %in% c("variety", "subspecies", "species"), Accepted_species, 
                                       ifelse(Accepted_name_rank %in% c("genus", "family"), species, Accepted_name
                                              )
                                       )
                                ),
         Accepted_name_rank = ifelse(Accepted_name_rank %in% c("subspecies", "variety"), "species", Accepted_name_rank)) %>%
  select(species_raw, species, species_final, Accepted_name_rank) %>%
  unique() %>%
  rename("species_precleaned"   = "species",
         "species_harmonized"   = "species_final", 
         "identification_level" = "Accepted_name_rank") %>%
  mutate(species_raw = factor(species_raw),
         species_precleaned = factor(species_precleaned),
         
         species_harmonized = gsub(species_harmonized, pattern = "XXX", replacement = ""),
         species_harmonized = factor(species_harmonized),
         
         identification_level = case_when(species_precleaned %in% c("Animal diggings", "Bare ground") ~ "", .default = identification_level),
         identification_level = factor(identification_level),
         
         broad_category = identification_level,
         broad_category = case_when(species_precleaned %in% c("Bare ground",
                                                              "LitterXXX",
                                                              "FungiXXX",
                                                              "Dung",
                                                              "Lichen",
                                                              "Something else",
                                                              "Animal diggings") ~ "non_living_plant_things",
                                    
                                    species_precleaned %in% c("Acer spp.",
                                                              "Conyza spp.",
                                                              "Crataegus spp.",
                                                              "Cyperus spp.",
                                                              "Euphorbia spp.",
                                                              "Fabaceae spp.",
                                                              "Lonicera spp.",
                                                              "Malus spp.",
                                                              "Morus spp.",
                                                              "Prunus spp.",
                                                              "Unknown forb spp.",
                                                              "Unknown grassXXX spp.",
                                                              "Unknown spp.",
                                                              "Unknown woodyXXX spp.") ~ "unidentified_mixed_batches",
                                    
                                    species_precleaned %in% c("Unknown woodyXXX sp.",
                                                              "Unknown forb sp.",
                                                              "Unknown grassXXX sp.",
                                                              "Unknown sp.",
                                                              "Bryophyte",
                                                              "Pteridophyta sp.") ~ "unidentified_individual_species",
                                    .default = broad_category),
         
         broad_category = factor(broad_category, levels = c("species", "genus", "family", "unidentified_individual_species", "unidentified_mixed_batches", "non_living_plant_things")))



species_abundance <- species_abundance %>%
  merge(species_name_match_list %>%
          select(!species_precleaned),
        by.x = "species",
        by.y = "species_raw") 

write.csv(species_abundance, file.path(L2_dir, "./species_abundance_L1.csv"), row.names=F)
write.csv(species_name_match_list, file.path(L2_dir, "./species_name_list_L1.csv"), row.names=F) 
