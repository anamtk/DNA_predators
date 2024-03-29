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
ana_ug <- read.csv(here("data",
                        "raw_data",
                        "4_body_size_data",
                        "Pal_UG_mass_length.csv"))

#palmyra raw body size data
pal <- read.csv(here("data",
                     "raw_data",
                     "4_body_size_data", 
                     "palmyra_food_web_sizes.csv"))

pal <- pal %>%
  dplyr::select(Species_Name,
                Length_mm,
                Total_Mass_mg) %>%
  rename(Mass_mg = Total_Mass_mg) %>%
  filter(!is.na(Mass_mg))

pal_IDs <- read.csv(here("data",
                         "raw_data",
                         "4_body_size_data",
                         "palmyra_nodes.csv"))

pal_IDs <- pal_IDs %>%
  dplyr::select(Species_Stage_Name,
                Phylum,
                Class,
                Order.1,
                Family,
                Genus,
                specific_epithet) %>%
  rename(Species_Name = Species_Stage_Name,
         Order = Order.1)

#how many species per family on Palmyra (justification of family-level assignment):
pal_IDs %>%
  group_by(Family) %>%
  tally(name= "total") %>%
  filter(Family != "") %>%
  summarise(mean= mean(total),
            n = n(),
            sd = sd(total),
            se = sd/sqrt(n),
            max = max(total))

#here is where some data go away
pal <- pal %>%
  left_join(pal_IDs, by = "Species_Name")

#literature
pantala_lit <- read.csv(here("data",
                             "raw_data",
                             "4_body_size_data",
                             "su_pantala_mass_length.csv"))
sohlstrom <- read.csv(here("data",
                           "raw_data",
                           "4_body_size_data",
                           "Sohlstroem_mass_length.csv"))

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
 
#give source for Palmyra data
pal <- pal %>%
  mutate(Source = "Palmyra") %>%
  rename("Species" = "Species_Name") %>%
  filter(Species != "Heteropoda_venatoria")

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
  bind_rows(pal) %>%
  bind_rows(pantala_lit) %>%
  bind_rows(sohlstrom)

all_bs %>%
  group_by(Source) %>%
  tally() 

#############################
# Subset only predators --------
#############################
bs_pred <- all_bs %>%
  filter(Genus %in% c("Phisis", "Heteropoda","Opopaea", "Scytodes",
                        "Neoscona", "Pantala") |
           Family %in% c("Mecistocephalidae", "Sparassidae", "Oonopidae",  "Scytodidae", "Pholcidae", "Libellulidae", "Anisolabididae") |
  Species %in% c("Neoscona_theisi", "Phisis_holdhausi", "Heteropoda_venatoria"))

bs_pred %>%
  distinct(Family)

bs_pred %>%
  distinct(Family, Source) %>%
  group_by(Source) %>%
  tally()
#############################
# Subset only prey in DNA --------
#############################

#master prey family list
pal_fams <- read.csv(here("data", 
                          "outputs", 
                          "2e_rarefied_taxonomic_sort", 
                          "fam_prey_DNA_conservative.csv"))

#remove the sequencing run from sample name
pal_fams$sample <- str_sub(pal_fams$sample, end=-2)

prey_ids <- pal_fams %>%
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
bs_prey %>%
  distinct(Family, Source) %>% 
  group_by(Source) %>%
  tally()
#############################
# Export --------
#############################

write.csv(bs_pred, here("data", 
                        "outputs", 
                        "2g_master_body_size_lists", 
                        "pred_mass_length.csv"))

write.csv(bs_prey, here("data", 
                        "outputs", 
                        "2g_master_body_size_lists",  
                        "prey_mass_length.csv"))

write.csv(all_bs, here("data",
                       "outputs",
                       "2g_master_body_size_lists",
                       "all_mass_length.csv"))




