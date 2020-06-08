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
#Manipulate Data to Long and Attach Predator IDs ####
###########################
#make long by sample and read abundance
comm_long <- comm %>%
  gather(sample, reads, HEV01a:HEV99d)

#this  binds to the predator ID DF, detecting the string from the sample_str
#column in the predator DF in each sample name. 
comm_long <- comm_long %>%
  fuzzy_inner_join(preds, by = c("sample" = "sample_str"), match_fun = str_detect)

###########################
#Attach prey IDs for each ASV ####
###########################
taxa_comm <- comm_long %>%
  left_join(taxa, by = "ASV")


###########################
#subset to sort ####
###########################
#Right now  this dataset includes a few kinds of data
#1. ASVs that are not predators
#2. ASVs where pred_ID and unique_ID match
#3. ASVs where pred_ID and some level of taxonomic assignment match
#4. ASVs where pred_ID and some level of BOLD assignment match

#ASVs that ARE predators (some of these are duplicate right now)

a <- taxa_comm %>%
  dplyr::filter(pred_ID == unique_ID)

b <- taxa_comm %>%
  filter(pred_ID == ID_bold)

c <- taxa_comm %>%
  filter(pred_Genus == unique_ID)

d <- taxa_comm %>%
  filter(pred_Family == unique_ID)

a %>%
  anti_join(b)

#these both led to no new taxonomies
#taxa_comm %>%
#  filter(pred_Genus == ID_bold & ID_bold != "")

#taxa_comm %>%
#  filter(pred_Family == ID_bold & ID_bold != "")



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

