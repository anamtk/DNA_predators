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


