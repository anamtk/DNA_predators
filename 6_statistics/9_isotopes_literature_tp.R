###########################
# Isotopes to literature trophic position
# July 22, 2020
# Ana Miller-ter Kuile
###########################

#corrected isotope data for the predators that have DNA too
#families in the guts of those spiders
#trophic positions of those families

###########################
# Load packages
library(here)
library(tidyverse)
library(ggplot2)
library(glmmTMB)
library(emmeans)
library(DHARMa)
library(effects)
###########################

###########################
# Load data
###########################
pal <- read.csv(here("data", "outputs",
                     "5_rarefied_taxonomic_sort",
                     "fam_prey_DNA.csv"))

all_fams <- read.csv(here("data", "outputs", "7_all_webs", 
                          "diet_families_tp.csv"))

meta <- read.csv(here("data", "Sample_metadata.csv"))

isotopes <- read.csv(here('data', "isotopes", "2015 Palmyra Cane Spider Isotopes.csv"))

baseline <- read.csv(here('data', "isotopes", "2009-2012 Palmyra baseline isotopes.csv"))

###########################
# Manipulate data
###########################

meta <- meta %>%
  filter(Isotope_ID != "NA") 

samples <- as.vector(meta$Extraction.ID)

pal$sample <- str_sub(pal$sample, end=-2)

pal <- pal %>%
  filter(reads > 0) %>%
  filter(sample %in% samples) %>%
  filter(Class != "Mammalia") %>%
  dplyr::select(sample, Class, Order, Family) %>%
  left_join(all_fams, by = c("Family" = "taxon_Family")) %>%
  dplyr::select(-X)

###########################
# Correct Isotopes to food chain length
###########################

#TP = 1+(d15Nconsumer - d15Nbase)/delta

#d15Nbase = d15Nplants*alpha + (d15Nmarinewrack(1-alpha))/delta

#alpha = (d13Cconsumer - d13C marine wrack)/(d13Cplant -d13Cmarine wrack)

#delta = 3.4 per mil, or 0.0034

#for each consumer, new rows:
#d15N plants, d15Nmarine, d13C marine, d13c plant

baseline$source <- ifelse(baseline$Species == "Algae", "marine", 
                          ifelse(baseline$Species == "Guano", "birds", "terrestrial"))

base <- baseline %>%
  group_by(source, Island.name) %>%
  summarise(d15N = mean(d15N, na.rm=T),
            d13C = mean(d13C, na.rm=T)) %>%
  filter(source != "birds") %>%
  pivot_wider(
    names_from = source,
    values_from = c(d15N, d13C)
  ) 

marineC <- base$d13C_marine[1]

marineN <- base$d15N_marine[1]

terrestrial <- base %>%
  dplyr::select(-d13C_marine, -d15N_marine) %>%
  filter(Island.name != "")

isotopes$Island <- sapply(str_split(isotopes$ID, " "), function(x){return(x[[1]])})

isotopes <- isotopes %>%
  dplyr::select(ID, d15N, d13C, Island) %>%
  rename("d15N_consumer" = "d15N",
         'd13C_consumer' = 'd13C') %>%
  left_join(terrestrial, by = c("Island" = "Island.name")) %>%
  mutate(d13C_marine = marineC,
         d15N_marine = marineN)

#TP = 1+(d15Nconsumer - d15Nbase)/delta

#d15Nbase = d15Nplants*alpha + (d15Nmarinewrack(1-alpha))/delta

#alpha = (d13Cconsumer - d13C marine wrack)/(d13Cplant -d13Cmarine wrack)

#delta = 3.4 per mil, or 0.0034

iso <- isotopes %>%
  mutate(alpha = (d13C_consumer - d13C_marine)/(d13C_terrestrial - d13C_marine),
         d15N_base = ((d15N_terrestrial*alpha) + (d15N_marine*(1-alpha))/3.4),
         iso_TP = 1 + (d15N_consumer = d15N_base)/3.4) %>%
  left_join(meta, by = c("ID" = "Isotope_ID")) %>%
  filter(!is.na(Extraction.ID))

###########################
# TL of diet
###########################

pal$tl <- ifelse(pal$broad_tp == "basal", 2,
                 ifelse(pal$broad_tp == "omnivorous", 3, 4))

tp_df <- pal %>%
  group_by(sample) %>%
  summarise(mean_tl = mean(tl, na.rm=T),
            sd= sd(tl, na.rm=T),
            total = n(),
            se = sd/sqrt(total)) %>%
  left_join(iso, by = c("sample" = "Extraction.ID")) %>%
  filter(!is.na(d15N_consumer)) 

tp_df$se[is.na(tp_df$se)] <- 0

ggplot(tp_df, aes(x = iso_TP, y = mean_tl)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_tl - se, ymax = mean_tl + se)) +
  ylim(2,5) +
  xlim(2,5) +
  geom_segment(aes(x = 2, y = 2, xend = 5, yend = 5, colour = "segment")) +
  theme_bw() +
  labs(x = "Isotopic Trophic Position", y = "Literature Trophic Position") +
  theme(legend.position = "none")

View(tp_df)  

