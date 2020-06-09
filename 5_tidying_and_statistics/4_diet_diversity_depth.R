###########################
#Diet Diversity and Depth
#Ana Miller-ter Kuile
#June 9, 2020
###########################

#This code will look at how diet diversity by predator is represented
#in our dataset
#Do some predator species have broader diets?
#Do some predator individuals of different species have more in an individual diet?
#Did we capture the majority of the diet of our predators?
#
##########################
#Load packages####
###########################
library(here)
library(tidyverse)
library(ggplot2)
library(iNEXT)
library(glmmTMB)
library(effects)
library(emmeans)
library(DHARMa)
library(MuMIn)

##########################
#Load data####
###########################
#data from the taxonomic sort, which is already concatenated
#by unique ID within a predator individual
dna <- read.csv(here("data", "outputs", "5_rarefied_taxonomic_sort", "prey_DNA.csv"))

meta <- read.csv(here("data", "Sample_metadata.csv"))

##########################
#Tidy data for analyses####
###########################
dna2 <- dna
dna2$sample <- str_sub(dna2$sample, start = 1, end =-2) #remove "a" at end

species_counts <- meta %>%
  distinct(ID, Extraction.ID, No..Individuals) %>%
  group_by(ID) %>%
  summarise(individuals = sum(No..Individuals))
  
#Do some predator species have broader diets?
species <- dna2 %>%
  filter(reads > 0) %>%
  group_by(pred_ID) %>%
  summarise(species = n_distinct(unique_ID)) %>%
  left_join(species_counts, by = c("pred_ID" = "ID")) %>%
  rename("pred_individuals" = "individuals")

#Do some predator individuals of different species have more in an individual diet?
indivs <- dna %>%
  filter(reads > 0) %>%
  group_by(sample, pred_ID) %>%
  tally(name = "species") %>%
  mutate(run = substr(sample, nchar(sample)-1+1, nchar(sample)))

indivs$sample <- str_sub(indivs$sample, start = 1, end =-2) #remove "a" at end

indivs <- indivs %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  distinct(sample, pred_ID, species, run, Method, Island, Habitat, Microhabitat,
           Year, Date.Collected, No..Individuals)

#Did we capture the majority of the diet of our predators?
#NEED: SOmethign that iNext can capture... have to think about what this looks like
#need columns to be summed by predator species, which I think i'll first do
#by total read counts, and then by incidence? yes
pred_abund <- dna %>%
  group_by(pred_ID, unique_ID) %>%
  summarise(reads = sum(reads, na.rm =T)) %>%
  pivot_wider(names_from = pred_ID,
              values_from = reads) %>%
  column_to_rownames(var = "unique_ID")

pred_abund[is.na(pred_abund)] <- 0
pred_abund <- pred_abund[rowSums(pred_abund) != 0,] #removed 16

pred_pres <- dna %>%
  group_by(pred_ID, unique_ID) %>%
  summarise(presence = sum(reads > 0)) %>%
  pivot_wider(names_from = pred_ID,
              values_from = presence) %>%
  column_to_rownames(var = "unique_ID")

pred_pres[is.na(pred_pres)] <- 0
pred_pres <- pred_pres[rowSums(pred_pres) != 0,] #removed 16
##########################
#quick vis####
###########################
ggplot(species, aes(x = pred_ID, y = species/pred_individuals)) +
  geom_point()

ggplot(indivs, aes(x = pred_ID, y = species/No..Individuals)) +
  geom_boxplot() +theme_bw()

##########################
#iNEXT####
###########################
#need to rethink this part a bit... specifically, need to do:
#number of observations, by number of new discovered diet species
abund_depth <- iNEXT(pred_abund, q=0, datatype="abundance") #this determines sequencing depth for each sample

abund_depth$DataInfo$SC
abund <- abund_depth$DataInfo

#graph the interpolated and extrapolated sampling depth per sample
ggiNEXT(abund_depth, type=1, facet.var="none") + 
  # can set se = F to remove shaded regions to see lines better 
  theme_bw() +
  labs(x = "Sampling Depth", y = "Species Abundance", title = "Number of Individuals") +
  theme(axis.text = element_text(size = 20), 
        axis.title = element_text(size = 25))

pres_depth <- iNEXT(pred_pres, q=0, datatype="abundance")

pres_depth$DataInfo$SC
pres <- pres_depth$DataInfo

#graph the interpolated and extrapolated sampling depth per sample
ggiNEXT(pres_depth, type=1, facet.var="none") + 
  # can set se = F to remove shaded regions to see lines better 
  theme_bw() +
  labs(x = "Sampling Depth", y = "Species Presence", title = "Number of Individuals") +
  theme(axis.text = element_text(size = 20), 
        axis.title = element_text(size = 25))

