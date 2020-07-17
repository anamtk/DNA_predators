###########################
#Subset species taxonomies
#Ana Miller-ter Kuile
#July 17, 2020
###########################

#This is the species-level data only from the taxonomic 
#assignments (fewer ASVS)

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
taxa <- read.csv(here("data", "outputs", "1_taxonomic_assignment", "species_taxonomies.csv"))
taxa <- taxa %>%
  dplyr::select(-X)

preds <- read.csv(here("data", "Predator_IDs.csv"))
#load community data
comm <- read.csv(here("data", "outputs", "4_rarefied", "community_rare.csv"))

comm <- comm %>%
  rename("ASV" = "X")

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
  left_join(taxa, by = "ASV") %>%
  filter(!is.na(Species)) %>% #remove ASVs with no species-level taxonomy
  mutate(run = substr(sample, nchar(sample)-1+1, nchar(sample))) #indicates the 
#run it was run on

###########################
#find those where pred_ID = taxonomy ####
###########################

pred0 <- taxa_comm %>%
  filter(pred_ID == Species | pred_ID == Order)

#being conservative about things that *may* be predator DNA
pred1 <- taxa_comm %>%
  filter(sample_str == "LRS" & Family == "Oonopidae")
  
pred2 <- taxa_comm %>%
  filter(sample_str == "SME" & Genus == "Carapoia") 

pred3 <- taxa_comm %>%
  filter(sample_str == "HEV" & Family == "Sparassidae")

pred4 <- taxa_comm %>%
  filter(sample_str == "LRS" & Family == "Theridiidae")

pred5 <- taxa_comm %>%
  filter(sample_str == "NEO" & Family == "Theridiidae")

#bind them all together
pred <- pred0 %>%
  bind_rows(pred1) %>%
  bind_rows(pred2) %>%
  bind_rows(pred3) %>%
  bind_rows(pred4) %>%
  bind_rows(pred5)
  
#anti-join for all the prey
not_pred <- taxa_comm %>%
  anti_join(pred, by = c("sample", "ASV"))

#look at the species = 51
not_pred %>%
  filter(reads > 0) %>%
  distinct(Species)

###########################
#Export ####
###########################

write.csv(pred, here("data", "outputs", "5_rarefied_taxonomic_sort", "species_predator_DNA.csv"))

write.csv(not_pred, here("data", "outputs", "5_rarefied_taxonomic_sort", "species_prey_DNA.csv"))
