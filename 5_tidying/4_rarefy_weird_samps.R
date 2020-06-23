###########################
#Rarefying weird samples
#Ana Miller-ter Kuile
#June 4, 2020
###########################

#This code will rarefy the weird smaples that didn't
#sequence deeply. Perhaps they still have useful diet 
#info in them? Worth exploring

#for this, will probably only look at the ISO samples
#since I can't really compare the others to low 
#sequenced samples


###########################
#Load packages
###########################
library(here)
library(tidyverse)
library(vegan)

###########################
#Load data
###########################

data <- read.csv(here("data", "outputs", "3_depth_corrected", "all_samples.csv"))

iso <- data %>%
  dplyr::select(ASV, ISO10c, ISO11c, ISO12c, ISO13c, ISO14c, ISO15c, ISO16c, 
                ISO17c,  ISO18c,  ISO19c , ISO1c,   ISO20c,  ISO21c , ISO22c,
                ISO23c , ISO24c, ISO25c,  ISO26c , ISO27c , ISO2c ,  ISO30c ,
                ISO31c,  ISO32c , ISO33c,  ISO34c ,ISO35c,  ISO36c,  ISO37c,
                ISO4c,   ISO5c,   ISO6c,   ISO7c ,  ISO8c,   ISO9c) %>%
  column_to_rownames(var = "ASV")

###########################
#Assess sequencing depth variation
###########################

min(colSums(iso)) #31
max(colSums(iso)) #8843
hist(colSums(iso))

###########################
#Rarefy
###########################
#we will rarefy based on lowest sample read abundance for each dataset
iso_rare_low <- min(colSums(iso))

#rarefaction is a random process, so set seed for consistent results
set.seed(1)

#Rarefy with rrarefy from vegan package
iso_rare <- as.data.frame(t(rrarefy(t(iso), sample = iso_rare_low)))

###########################
#Export for Subsetting
###########################
write.csv(iso_rare, here("data", "outputs", "4_rarefied", "iso_rare.csv"))
