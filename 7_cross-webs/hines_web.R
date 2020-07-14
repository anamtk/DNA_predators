###########################
# Hines et al. 2019 web cleaning
# July 13, 2020
# Ana Miller-ter Kuile
###########################

#This code takes data from the Hines et al. 2019 food web
#from rmangal and cleans it up into several usable datasets
#for my analyses. 

#Specifically, outputs will be:
#1. Predation links of predators with prey concatenated at 
#family level
#2. Predation links per predator species 

#I've exported the hines data from rmangal already in another
#script, and have used those node data to assign family-level
#taxonomies to each species in that dataset. 

###########################
# Load packages
library(here)
library(tidyverse)
library(ggplot2)
###########################


###########################
# Load data
###########################

interactions <- read.csv(here("Published_webs", "final_cut",
                              "hines_2019", "interactions_hines.csv"))

interactions <- interactions %>%
  dplyr::select(node_from, node_to, type, method)

nodes <- read.csv(here("Published_webs", "final_cut",
                       "hines_2019", "nodes_hines.csv"))

nodes <- nodes %>%
  dplyr::select(node_id, original_name, taxon_Kingdom,
                taxon_Class, taxon_Order, taxon_Family)

###########################
# All predation links
###########################

predation <- interactions %>%
  filter(type == "predation") %>%
  left_join(nodes, by =c("node_to" = "node_id"))

###########################
# Per predator links
###########################
per_pred <- predation %>%
  group_by(node_from, taxon_Family) %>%
  tally(name = "links") %>%
  group_by(node_from) %>%
  summarise(links = n())

per_pred %>%
  summarise(min = min(links),
            max = max(links),
            mean = mean(links), 
            total = n(), 
            sd = sd(links),
            se = sd/sqrt(total))

###########################
# Export data
###########################
write.csv(predation, here("data", "outputs", 
                          "6_pub_webs", "hines_predation.csv"))

write.csv(per_pred, here("data", "outputs", 
                          "6_pub_webs", "hines_per_pred.csv"))






