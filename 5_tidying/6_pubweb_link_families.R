###########################
# Links families
# July 22, 2020
# Ana Miller-ter Kuile
###########################

#going to be looking at how many families are in our
#links across webs. 

###########################
# Load packages
library(here)
library(tidyverse)
library(ggplot2)
###########################

###########################
# Load data
###########################

pal <- read.csv(here("data", "outputs",
                     "5_rarefied_taxonomic_sort",
                     "fam_prey_DNA.csv"))


pub <- read.csv(here("data", "outputs", "6_pub_webs", "pub_webs_predation.csv"))

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
