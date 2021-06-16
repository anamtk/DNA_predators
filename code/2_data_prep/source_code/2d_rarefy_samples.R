###########################
#Rarefying datasets
#Ana Miller-ter Kuile
#June 4, 2020
###########################

#This code will rarefy our datasets which we have already subset 
#for both cross-run comparisons and for community analyses


###########################
#Load packages
###########################
package.list <- c("here", "tidyverse", "vegan")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

###########################
#Load data
###########################

#Cross-run comparsion dataset
cross <- read.csv(here("data", 
                       "outputs", 
                       "2b_depth_corrected",
                       "cross_run_samples.csv"))
cross <- cross %>%
  column_to_rownames(var = "ASV") %>%
  dplyr::select(-X)

#Community dataset
community <- read.csv(here("data", 
                           "outputs", 
                           "2b_depth_corrected",
                           "depth_subset_samples.csv"))

community <- community %>%
  column_to_rownames(var = "ASV") %>%
  dplyr::select(-X)

###########################
#Assess sequencing depth variation
###########################
min(colSums(cross)) #22544
max(colSums(cross)) #158055
hist(colSums(cross)) #wow, pretty normal!


min(colSums(community)) #15954
max(colSums(community)) #236393
hist(colSums(community))

###########################
#Rarefy
###########################
#we will rarefy based on lowest sample read abundance for each dataset
cross_rare_low <- min(colSums(cross))
community_rare_low <- min(colSums(community))

#rarefaction is a random process, so set seed for consistent results
set.seed(1)

#Rarefy with rrarefy from vegan package
cross_rare <- as.data.frame(t(rrarefy(t(cross), sample = cross_rare_low)))
comm_rare <- as.data.frame(t(rrarefy(t(community), sample = community_rare_low)))

###########################
#Export for Subsetting
###########################
write.csv(cross_rare, here("data", 
                           "outputs", 
                           "2d_rarefied", 
                           "cross_run_rare.csv"))

write.csv(comm_rare, here("data", 
                          "outputs", 
                          "2d_rarefied", 
                          "community_rare.csv"))

