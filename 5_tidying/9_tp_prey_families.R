###########################
# Links families Trophic Position Assignment
# July 22, 2020
# Ana Miller-ter Kuile
###########################

#going to be looking at the trophic positions
#of diet items across different webs
#this script imports family data, exports
#for manual assignment via the literature
#and then exports the final document with
#various different trophic levels

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


pub <- read.csv(here("data", "outputs", 
                     "6_pub_webs", 
                     "pub_webs_predation.csv"))

###########################
# Unique families
###########################

pal_fams <- pal %>%
  distinct(Family) %>%
  rename("taxon_Family" = "Family")

all_fams <- pub %>%
  distinct(taxon_Family) %>%
  bind_rows(pal_fams)

###########################
# OUtput to Google trophic groups
###########################

write.csv(all_fams, here("data", "outputs", "7_all_webs", 
                         "diet_families.csv"))

###########################
# Import the updated version
###########################

all_fams <- read.csv(here("data", "outputs", "7_all_webs", 
                          "diet_families_tp.csv"))

###########################
# Combine Pal and published web data
###########################
pub <- pub %>%
  dplyr::select(-X, -type, -taxonomy.name, -taxon_Kingdom,
                -taxon_Class, -taxon_Order, )

pal %>%
  distinct(Family) %>%
  tally(name = "family_richness")

pal <- pal %>%
  group_by(pred_ID, Family) %>%
  summarise(reads = sum(reads)) %>%
  filter(reads > 0) %>%
  mutate(web = "Palmyra",
         species_richness = 409,
         family_richness = 69,
         coll_method = "HTS molecular",
         pub_year = 2020,
         original_name = Family,
         resource = Family) %>%
  rename("consumer" = "pred_ID",
         "taxon_Family" = "Family") %>%
  dplyr::select(-reads)

all_fam <- bind_rows(pal, pub) %>%
  left_join(all_fams, by = "taxon_Family") %>%
  dplyr::select(-X) %>%
  unite(tp, trophic_group, predatory_omnivory, sep = "", remove = FALSE)

###########################
# Create trophic position variables
###########################

all_fam$tp <- factor(all_fam$tp, levels = c("herbivore", "detritovore", "omnivoreno",
                                            "omnivoreyes", "parasite", "parasitoid", "predator"))

all_fam$broad_tp <- ifelse(all_fam$tp == "herbivore" | all_fam$tp == "detritovore" | all_fam$tp == "omnivoreno", "basal",
                           ifelse(all_fam$tp == "omnivoreyes", "omnivorous", 
                           ifelse(is.na(all_fam$tp), NA, "predatory")))

###########################
# Export
###########################

write.csv(all_fam, here("data", "outputs", "7_all_webs", "all_interactions_and_tp.csv"))
