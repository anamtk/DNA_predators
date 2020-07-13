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