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

e <- taxa_comm %>%
  filter(pred_Order == unique_ID & pred_Order == "Geophilomorpha")

f <- taxa_comm %>%
  filter(pred_Order == ID_bold & pred_Order == "Odonata")

g <- taxa_comm %>%
  filter(pred_Order == ID_bold & pred_Order == "Dermaptera")

predator_ASVs <- a %>%
  bind_rows(b) %>%
  bind_rows(c) %>%
  bind_rows(d) %>%
  bind_rows(e) %>%
  bind_rows(f) %>%
  bind_rows(g) %>%
  group_by(ASV, sample, reads, pred_ID, unique_ID) %>%
  distinct()

#230 total samples
all_taxa <- taxa_comm %>% 
  ungroup() %>%
  distinct(sample)

#176 of 230 samples got a predator taxonomy... i don't believe it...
assigned <- predator_ASVs %>%
  ungroup() %>%
  distinct(sample)

#these are those samples
unassigned <- taxa_comm %>%
  anti_join(assigned, by = "sample")

#MY approach will be: remove the predator ASVs from other samples
#Then export the per ASV read abundance, fit a distribution (like Jerde)
#and then use this to inform a high-end cutoff for the samples
#that didn't get an assignmen

###########################
#subset by species to sort ####
###########################

#For each predator, may need to revisit some of the high-read things that
#DON'T match to family or lower, as these are most likely predator DNA
#My thought for correcting for this is to create a master predator DNA
#DF, and then fit a distribution to these data, then use  this (ala Jerde)
#to predict if ASV reads above a certain amount should be subset as predator
#if they match to order or Class of the predator species in question

###########################
##HEV ####
###########################

#subset HEV predators
HEV <- taxa_comm %>%
  filter(sample_str == "HEV")

#ID matched to species?
#HEV %>% 
#  filter(pred_ID == unique_ID)

#ID matched to genus?
#HEV %>% 
#  filter(pred_Genus == unique_ID)

#ID matched to family?
HEV %>% 
  filter(pred_Family == unique_ID)

#ID matched on bold?
HEV %>% 
  filter(pred_ID == ID_bold)

###########################
##NEO ####
###########################
#NEO
NEO <- taxa_comm %>%
  filter(sample_str == "NEO")

#ID matched to species?
NEO %>% 
  filter(pred_ID == unique_ID)

#ID matched to genus?
NEO %>% 
  filter(pred_Genus == unique_ID)

#ID matched to family?
#NEO %>% 
 # filter(pred_Family == unique_ID)

#ID matched to bold?
NEO %>% 
  filter(pred_ID == ID_bold)

###########################
##SCY ####
###########################
#SCY
SCY <- taxa_comm %>%
  filter(sample_str == "SCY")

#ID matched to species?
#SCY %>% 
#  filter(pred_ID == unique_ID)

#ID matched to genus?
#SCY %>% 
#  filter(pred_Genus == unique_ID)

#ID matched to family?
#SCY %>% 
#filter(pred_Family == unique_ID)

#ID matched to bold?
SCY %>% 
  filter(pred_ID == ID_bold)

###########################
##CEN ####
###########################
#CEN 
#subset centipede
CEN <- taxa_comm %>%
  filter(sample_str == "CEN")

#ID matched to order?
CEN %>% 
  filter(pred_ID == unique_ID)

#ID matched to bold?
#CEN %>% 
#  filter(pred_ID == ID_bold)

###########################
##PHH ####
###########################

#PHH subset 
PHH <- taxa_comm %>%
  filter(sample_str == "PHH")

#ID matched to species?
#PHH %>% 
#  filter(pred_ID == unique_ID)

#ID matched to genus?
##PHH %>% 
 #filter(pred_Genus == unique_ID)

#ID matched to family?
#PHH %>% 
# filter(pred_Family == unique_ID)

#ID matched to order? - not sure i trust this, want to be able to distinguish pred from non
PHH %>%
  filter(pred_Order == unique_ID)

#ID matched to bold order?
PHH %>% 
  filter(pred_Order == ID_bold)

###########################
##SME ####
###########################

#SME subset 
SME <- taxa_comm %>%
  filter(sample_str == "SME")

#ID matched to species?
SME %>% 
  filter(pred_ID == unique_ID)

#ID matched to genus?
#SME %>% 
# filter(pred_Genus == unique_ID)

#ID matched to family?
#SME %>% 
# filter(pred_Family == unique_ID)

#ID matched to bold?
SME %>% 
  filter(pred_ID == ID_bold)

###########################
##EUB ####
###########################
#EUB subset 
EUB <- taxa_comm %>%
  filter(sample_str == "EUB")

#ID matched to species?
#EUB %>% 
#  filter(pred_ID == unique_ID)

#ID matched to genus?
##EUB %>% 
# filter(pred_Genus == unique_ID)

#ID matched to family?
#EUB %>% 
# filter(pred_Family == unique_ID)

#Order?
#EUB %>%
#  filter(pred_Order == unique_ID)

#ID matched to bold?
#EUB %>% 
#  filter(pred_ID == ID_bold)

#Order matched bold ID?
EUB %>% 
  filter(pred_Order == ID_bold)

###########################
##LRS ####
###########################

#LRS subset 
LRS <- taxa_comm %>%
  filter(sample_str == "LRS")

#ID matched to species?
LRS %>% 
  filter(pred_ID == unique_ID)

#ID matched to genus?
##LRS %>% 
 #filter(pred_Genus == unique_ID)

#ID matched to family?
#LRS %>% 
# filter(pred_Family == unique_ID)

#ID matched to bold?
#LRS %>% 
#  filter(pred_ID == ID_bold)

###########################
##PAN ####
###########################
#PAN subset 
PAN <- taxa_comm %>%
  filter(sample_str == "PAN")

#ID matched to species?
PAN %>% 
  filter(pred_ID == unique_ID)

#ID matched to genus?
#PAN %>% 
# filter(pred_Genus == unique_ID)

#ID matched to family?
##PAN %>% 
#filter(pred_Family == unique_ID)

#Order?
#PAN %>%
#  filter(pred_Order == unique_ID)

#ID matched to bold?
PAN %>% 
  filter(pred_ID == ID_bold)

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

