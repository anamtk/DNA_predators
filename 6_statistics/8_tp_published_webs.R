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
  group_by(consumer, web, family_richness, coll_method,
           pub_year, tp) %>%
  tally(name = "link_number")

meta_tp <- links_tp %>%
  dplyr::select(consumer, web, family_richness, coll_method,
                pub_year)

tp_matrix <- links_tp %>%
  ungroup() %>%
  dplyr::select(consumer, tp, link_number) %>%
  pivot_wider(names_from = consumer,
              values_from = link_number) %>%
  filter(!is.na(tp))

tp_matrix[is.na(tp_matrix)] <- 0
