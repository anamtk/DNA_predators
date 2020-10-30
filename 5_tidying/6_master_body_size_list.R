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
sohlstrom <- read.csv(here("data", "size_data", "Sohlstroem_mass_length.csv"))

#############################
# Compile size data --------
#############################
#compiled dataframe needs:
#Class, Order, Family, Genus, Species, Length_mm, Mass_mg, Source

#subset and organize undergrad data
ana_ug <- ana_ug %>%
  dplyr::select(Order, Family, Genus, Species, Length_mm, Weight_mg) %>%
  rename("Mass_mg" = "Weight_mg") %>%
  mutate(Source = "Ana_UG",
         Class = ifelse(Order == "Orthoptera", "Insecta", "Arachnida"))
 
#subset and organize node-level data from palmyra
pal_nodes <- pal_nodes %>%
  dplyr::select(Class, Order.1, Family, Genus, specific_epithet, 
                Body_Length_Mean_mm, Body_Mass_Mean_mg) %>%
  rename("Order" = "Order.1",
         "Species" = "specific_epithet",
         "Length_mm" = "Body_Length_Mean_mm",
         "Mass_mg" = "Body_Mass_Mean_mg") %>%
  filter(!is.na(Mass_mg)) %>%
  mutate(Source = "Pal_nodes") 
  

#subset and organize individual level data from palmyra
#for this i prioritized individuals from predators since
#others are averaged for the nodes list above
pal_stages <- pal_stages %>%
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
         Source = "Pal_stages") %>%
  mutate(Mass_mg = as.numeric(Mass_mg))

#organize and subset the literature data
#Class, Order, Family, Genus, Species, Length_mm, Mass_mg, Source
pantala_lit <- pantala_lit %>%
  rename("Mass_mg" = "mass_mg",
         "Length_mm" = "length_mm") %>%
  mutate(Class = "Insecta")

sohlstrom <- sohlstrom %>%
  filter(Geographic_Region == "tropical") %>%
  dplyr::select(Order, Family, Live_mass_mg, length_mm) %>%
  rename("Mass_mg" = "Live_mass_mg",
         "Length_mm" = "length_mm") %>%
  mutate(Class = "",
         Source = "Sohlstrom et al. 2018") %>%
  mutate(Order = str_to_title(Order), #make first character uppercase!!
         Family = str_to_title(Family))

#final compilation of ALL         
all_bs <- ana_ug %>%
  bind_rows(pal_nodes) %>%
  bind_rows(pal_stages) %>%
  bind_rows(pantala_lit) %>%
  bind_rows(sohlstrom)

#############################
# Subset only predators --------
#############################
bs_pred <- all_bs %>%
  filter(Genus %in% c("Phisis", "Heteropoda","Opopaea", "Scytodes",
                        "Neoscona", "Pantala") |
           Family %in% c("Oonopidae",  "Scytodidae", "Pholcidae", "Libellulidae", "Anisolabididae") |
  Species %in% c("Neoscona_theisi", "Phisis_holdhausi", "Heteropoda_venatoria") |
  Order == "Geophilomorpha") 

#############################
# Subset only prey in DNA --------
#############################

#master prey family list
pal <- read.csv(here("data", "outputs",
                     "5_rarefied_taxonomic_sort",
                     "fam_prey_DNA_conservative.csv"))

#remove the sequencing run from sample name
pal$sample <- str_sub(pal$sample, end=-2)

prey_ids <- pal %>%
  distinct(Class, Order, Family)

bs_prey <- all_bs %>%
  filter(Family %in% prey_ids$Family |
           Order %in% c("Entomobryomorpha", "Poduromorpha", "Isopoda", "Psocoptera",
                        "Geophilomorpha", "Thysanoptera", "Psocoptera")) %>%
  add_row(Order = "", Family= "", Genus = "", Species = "", Length_mm = NA,
          Mass_mg = 0.000633, Source = "Yaninek 1993", Class = "Acari") %>%
   add_row(Order = "", Family= "", Genus = "", Species = "", Length_mm = NA,
          Mass_mg = 0.007035, Source = "Yaninek 1993", Class = "Acari")

#From Yaninek 1993, wet weight for an Acari:
#0.000633
#0.007035
#############################
# Export --------
#############################

write.csv(bs_pred, here("data", "outputs", "6_prey_sizes", "pred_mass_length.csv"))
write.csv(bs_prey, here("data", "outputs", "6_prey_sizes", "prey_mass_length.csv"))

