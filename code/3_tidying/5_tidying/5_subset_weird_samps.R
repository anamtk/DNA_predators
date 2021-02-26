###########################
#Subset ISO
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
taxa <- read.csv(here("data", "outputs", "1_taxonomic_assignment", "ASV_taxonomies_summed_wIndiv.csv"))
taxa <- taxa %>%
  dplyr::select(-X)

preds <- read.csv(here("data", "Predator_IDs.csv"))
#load community data
iso <- read.csv(here("data", "outputs", "4_rarefied", "iso_rare.csv"))

iso <- iso %>%
  rename("ASV" = "X")

###########################
#Manipulate Data to Long and Attach Predator IDs ####
###########################
#make long by sample and read abundance
iso_long <- iso %>%
  gather(sample, reads, ISO10c:ISO9c)

#this  binds to the predator ID DF, detecting the string from the sample_str
#column in the predator DF in each sample name. 
iso_long <- iso_long %>%
  fuzzy_inner_join(preds, by = c("sample" = "sample_str"), match_fun = str_detect)

###########################
#Attach prey IDs for each ASV ####
###########################
taxa_iso <- iso_long %>%
  left_join(taxa, by = "ASV") %>%
  filter(Domain == "Eukaryota") #remove NA taxonomies from this community

###########################
#sort taxnomies ####
###########################

#subset predator DNA from this:
ISO_pred <- taxa_iso %>%
  ungroup() %>%
  filter(Order == "Scorpiones") %>% #only order for scorpions
  filter(Family %in% c("Buthidae", "")) #either matched to the family or no family (being conservative)


ISO_prey <- taxa_iso %>%
  anti_join(ISO_pred, by = c("sample", "ASV"))

#based on non-zero hits here, my intuition is that most of the DNA in these 
#samples is that of other predators from the sequencing run these were run on
#we can keep this as an option, but it just seems questionable to me?