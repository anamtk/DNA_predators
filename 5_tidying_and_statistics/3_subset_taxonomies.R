###########################
#Subset taxonomies
#Ana Miller-ter Kuile
#June 8, 2020
###########################

#Need to subset prey taxonomies per predator from the rarefied dataset. 
#This will include giving the predators names that match the ASV ID, 
#so that I can filter out the ASVs where these match from any diet 
#analysis

###########################
#Load Packages ####
###########################
library(here)
library(tidyverse)
library(fuzzyjoin) #fuzzy inner join for string detection join

###########################
#Load Data ####
###########################
#Prey DNA taxonomies
taxa <- read.csv(here("data", "outputs", "1_taxonomic_assignment", "ASV_taxonomies.csv"))

preds <- read.csv(here("data", "Predator_IDs.csv"))
#load community data
comm <- read.csv(here("data", "outputs", "4_rarefied", "community_rare.csv"))

comm <- comm %>%
  rename("ASV" = "X.1") %>%
  dplyr::select(-X) 

###########################
#Manipulate Data to Long and Assign Predator IDs ####
###########################

comm_long <- comm %>%
  gather(sample, reads, HEV01a:HEV99d)

#Heteropoda venatoria
#Neoscona theisi
#Scyotodes Striatipes
#Centipede - Geophilomorpha, Tygarrup javanicus
#Phisis holdhausi, Tettigoniidae
#Smeringopus palidus
#Euborellia annulipes
#Pantala flavescens
#Keijia mneon, Platnickina mneon

#Geophilomorpha
#Keijia mneon
#Neoscona theisi
#Neoscona
#Pantala flavescens
#Smeringopus pallidus
#Sparassidae
#Tygarrup javanicus
#Tettigoniidae
comm_long <- comm_long %>%
  fuzzy_inner_join(preds, by = c("sample" = "sample_str"), match_fun = str_detect)
