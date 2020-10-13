###########################
# Laigle et al. 2017 web cleaning
# July 13, 2020
# Ana Miller-ter Kuile
###########################

#This code takes data from the Laigle et al. 2017 food web
#from Dryad and cleans it up into several usable datasets
#for my analyses. 

#Specifically, outputs will be:
#1. Predation links of predators with prey concatenated at 
#family level
#2. Predation links per predator species 
#3. keep predator and prey body sizes somewhere?
#4. interaction frequency somewhere from the cross-web stuff?

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
                              "Laigle_etal_2017",
                              "interactions_Laigle.csv"))

nodes <- read.csv(here("Published_webs", "final_cut",
                              "Laigle_etal_2017",
                              "nodes_Laigle.csv"))

prey_ids <- nodes %>%
  dplyr::select(X, Genus, Family, Order, Class, Phylum,
                Kingdom, sp)

###########################
# All predation links
###########################
predation <- interactions %>%
  left_join(nodes, by = c("resource" = "X")) %>%
  filter(Kingdom == "Animalia") %>%
  dplyr::select(-Web, -Poison, -herbivore, -carnivore,
                -fungivore, -detritivore, -above, -below,
                -RS1, -RS2) %>%
  filter(interaction > 0)

###########################
# Per predator links by prey family
###########################
per_pred <- predation %>%
  group_by(consumer, Family) %>%
  summarise(links = sum(interaction)) %>%
  group_by(consumer) %>%
  summarise(links = n())

per_pred %>%
  summarise(min = min(links),
            max = max(links),
            mean = mean(links), 
            total = n(), 
            sd = sd(links),
            se = sd/sqrt(total))

###########################
# Export predation links
###########################

write.csv(predation, here("data", "outputs", 
                          "6_pub_webs", "laigle_predation.csv"))

write.csv(per_pred, here("data", "outputs", 
                         "6_pub_webs", "laigle_per_pred.csv"))

###########################
# Predator/prey body sizes
###########################

View(interactions)

pred_size <- nodes %>%
  dplyr::select(X, mass) %>%
  rename("predator" = "X",
         "pred_mass" = "mass")

prey_size <- nodes %>%
  dplyr::select(X, mass) %>%
  rename("prey" = "X",
         "prey_mass" = "mass")

sizes <- interactions %>%
  left_join(pred_size, by = c("consumer" = "predator")) %>%
  left_join(prey_size, by = c("resource" = "prey"))

###########################
# Cross-web interaction frequency
###########################

frequency <- sizes %>%
  filter(interaction > 0) %>%
  group_by(consumer, resource, pred_mass, prey_mass) %>%
  summarise(frequency = sum(interaction))

write.csv(frequency, here("data", "outputs", 
                      "6_pub_webs", "laigle_size_and_frequency.csv"))
