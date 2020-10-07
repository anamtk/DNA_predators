###########################
# Master list of body sizes
# Ana Miller-ter Kuile
# October 7, 2020
###########################

# this script combines data from the Palmyra food web project,
# my undergrad thesis, and literature sources to make a 
# full list of all the body size (mass and length) for 
# all species in the predator and prey lists for my
# project. it also outputs a set of prey families for 
# which we have no size for literature sleuthing

###########################
# Load packages
package.list <- c("here", "tidyverse")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}
#############################

#############################
# Size datasets load --------
#############################

#undergrad size data for two species
ana_ug <- read.csv(here("data", "size_data", "Pal_UG_mass_length.csv"))

#palmyra node and stage data from FW project
pal_nodes <- read.csv(here("data", "size_data", "pal_nodes_mass_length.csv"))
pal_stages <- read.csv(here("data", "size_data", "pal_stages_mass_length.csv"))

#literature
pantala_lit <- read.csv(here("data", "size_data", "su_pantala_mass_length.csv"))
sohlstroem <- read.csv(here("data", "size_data", "Sohlstroem_mass_length.csv"))

#############################
# Compile size data --------
#############################
#compiled dataframe needs:
#Class, Order, Family, Genus, Species, Length_mm, Mass_mg, Source

ana_ug <- ana_ug %>%
  dplyr::select(Order, Family, Genus, Species, Length_mm, Weight_mg) %>%
  rename("Mass_mg" = "Weight_mg") %>%
  mutate(Source = "Ana_UG",
         Class = ifelse(Order == "Orthoptera", "Insecta", "Arachnida"))
 
pal_nodes <- pal_nodes %>%
  dplyr::select(Class, Order.1, Family, Genus, specific_epithet, 
                Body_Length_Mean_mm, Body_Mass_Mean_mg) %>%
  rename("Order" = "Order.1",
         "Species" = "specific_epithet",
         "Length_mm" = "Body_Length_Mean_mm",
         "Mass_mg" = "Body_Mass_Mean_mg") %>%
  filter(!is.na(Mass_mg)) %>%
  mutate(Source = "Pal_nodes")

#for this i prioritized individuals from predators since
#others are averaged for the nodes list above
pal_stages %>%
  dplyr::select(Class, Family, Order, Species_Name, 
                Length_mm, Mass_mg) %>%
  rename("Species" = "Species_Name") %>%
  filter(!is.na(Mass_mg)) %>%
  filter(Species %in% c("Heteropoda_venatoria", 
                             "Neoscona_theisi",
                             "Scytodes_longipes",
                             "Scytodes_sp_2",
                             "Oonopidae_sp1",
                             "Oonopidae_sp2",
                             "Oonopidae_sp3",
                             "Oonopidae_sp4",
                             "Oonopidae_sp5",
                             "Oonopidae_sp6",
                             "Opopaea_deserticola",
                             "Phisis_holdhausi")) %>%
  mutate(Genus = "",
         Source = "Pal_stages")

#master prey family list
pal <- read.csv(here("data", "outputs",
                     "5_rarefied_taxonomic_sort",
                     "fam_prey_DNA_conservative.csv"))

#remove the sequencing run from sample name
pal$sample <- str_sub(pal$sample, end=-2)

unique(pal$Family)

