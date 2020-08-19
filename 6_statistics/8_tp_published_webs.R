#############################
#Trophic positions of published vs. HTS links
#August 18, 2020
#Ana Miller-ter Kuile
#############################

#This script examines whether the trophic positions of diet 
#items assigned via different published methods differ
#from those from the HTS data
#This is important because of the importance of intraguild
#predation in food web ecology

#############################
#Load packages
library(here)
library(tidyverse)
library(ggplot2)
library(glmmTMB)
library(DHARMa)
library(MuMIn)
library(effects)
library(emmeans)
library(ggeffects)
library(vegan)
#############################

#note, it may be important to break this analysis up
#by first: more specific trophic groups and then
#by courser scale - e.g. are they eating plants, vs. 
#other animals, or both.

#############################
#Load data
#############################

links <- read.csv(here("data", "outputs", 
                       "7_all_webs", 
                       "all_interactions_and_tp.csv"))

#############################
#Summarize by specific TP
#############################

links_tp <- links %>%
  unite(web_consumer, consumer, web, remove = F) %>%
  group_by(web_consumer, web, family_richness, coll_method,
           pub_year, tp) %>%
  tally(name = "link_number") 

#############################
#Summarize by broad TP
#############################

links_btp <- links %>%
  unite(web_consumer, consumer, web, remove = F) %>%
  group_by(web_consumer, web, family_richness, coll_method,
           pub_year, broad_tp) %>%
  tally(name = "link_number") 

###########################
# other vis
###########################
links_btp <- links_btp %>%
  filter(broad_tp != "")

links_sp <- links %>%
  unite(web_consumer, consumer, web, remove = F) %>%
  group_by(web_consumer, web, family_richness, coll_method,
           pub_year) %>%
  tally(name = "total_link_number") %>%
  ungroup() %>%
  dplyr::select(web_consumer, total_link_number)

web_species_tp <- links %>%
  group_by(web, broad_tp) %>%
  distinct(resource) %>%
  summarise(total = n()) %>%
  filter(broad_tp != "") %>%
  pivot_wider(names_from = "broad_tp",
              values_from = "total")

web_species_tp[is.na(web_species_tp)] <- 0

links_btp <- links_btp %>%
  left_join(links_sp, by = "web_consumer") %>%
  left_join(web_species_tp, by = "web")

ggplot(links_btp, aes(x = coll_method, y = link_number/total_link_number, fill = broad_tp)) +
  geom_boxplot() + theme_bw()
