###########################
# Links families Trophic Position Assignment
# October 13, 2020
# Ana Miller-ter Kuile
###########################

#going to be looking at the trophic positions
#of diet items 
#this script imports family data, exports
#for manual assignment via the literature
#and then exports the final document with
#various different trophic levels

###########################
# Load packages
library(here)
library(tidyverse)
###########################

###########################
# Load data
###########################

pal <- read.csv(here("data", "outputs",
                     "5_rarefied_taxonomic_sort",
                     "fam_prey_DNA.csv"))

sizes <- read.csv(here("data", "outputs", "6_prey_sizes",
                       "DNA_interaction_pred_prey_sz.csv"))

###########################
# Unique families
###########################

pal_fams <- pal %>%
  distinct(Family) %>%
  rename("taxon_Family" = "Family")

###########################
# OUtput to Google trophic groups
###########################

write.csv(pal_fams, here("data", "outputs", "7_prey_tp", 
                         "diet_families.csv"))

###########################
# Import the updated version
###########################

all_fams <- read.csv(here("data", "outputs", "7_prey_tp", 
                          "diet_families_tp.csv"))

###########################
# Combine TP and body size data
###########################

sizes <- sizes %>%
  left_join(all_fams, by = c("Family" = "taxon_Family"))

###########################
# Export
###########################

write.csv(sizes, here("data", "outputs", "8_final_dataset", "pred_prey_sizes_tp_DNAinteractions.csv"))

