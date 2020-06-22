###########################
#Species Stage Structure
#Ana Miller-ter Kuile
#June 22, 2020
###########################

#This code will look at how diet changes/doesn't change with
#species body size. 

#For this script, it is probably a good idea to start explorations
#with the species for which we have large sample sizes and/or
#large ranges of body size

#In this logic - Euborellia, Phisis, and Heteropoda are probably first

##########################
#Load packages####
###########################
library(here)
library(tidyverse)
library(ggplot2)
library(vegan)

#########################
#Load data####
###########################
#data from the taxonomic sort, which is already concatenated
#by unique ID within a predator individual
dna_fam <- read.csv(here("data", "outputs", "5_rarefied_taxonomic_sort", "fam_prey_DNA.csv"))

#sample metadata for sample sizes
meta <- read.csv(here("data", "Sample_metadata.csv"))

#########################
#Heteropoda venatoria####
###########################
HEV <- dna_fam %>%
  filter(sample_str == "HEV")

