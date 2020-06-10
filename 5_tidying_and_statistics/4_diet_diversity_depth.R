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
library(vegan)

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
abund_depth$AsyEst
abund <- abund_depth$DataInfo

#graph the interpolated and extrapolated sampling depth per sample
ggiNEXT(abund_depth, type=1, facet.var="none") + 
  # can set se = F to remove shaded regions to see lines better 
  theme_bw() +
  labs(x = "Sampling Depth", y = "Species Richness") +
  theme(axis.text = element_text(size = 20), 
        axis.title = element_text(size = 25), legend.position = "none") 

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

##########################
#vegan####
###########################

sample_counts <- dna2 %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  distinct(ID, sample, No..Individuals) %>%
  dplyr::select(ID, sample, No..Individuals)
  

geo <- dna2 %>%
  filter(pred_ID == "Geophilomorpha") %>%
  dplyr::select(sample, unique_ID, reads) %>%
  group_by(sample, unique_ID) %>%
  summarise(reads = sum(reads)) %>%
  pivot_wider(names_from = "unique_ID",
              values_from = "reads") %>%
  column_to_rownames(var = "sample")
rowSums(geo)

geo1 <- specaccum(geo)
geo2 <- specaccum(geo, "random")

plot(geo1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
boxplot(geo2, col="yellow", add=TRUE, pch="+")

#"Geophilomorpha"       "Euborellia annulipes" "Heteropoda venatoria"
##[4] "Keijia mneon"         "Neoscona theisi"      "Pantala flavescens"  
#[7] "Phisis holdhausi"     "Scytodes longipes"    "Smeringopus pallidus"

eub <- dna2 %>%
  filter(pred_ID == "Euborellia annulipes") %>%
  dplyr::select(sample, unique_ID, reads) %>%
  group_by(sample, unique_ID) %>%
  summarise(reads = sum(reads)) %>%
  pivot_wider(names_from = "unique_ID",
              values_from = "reads") %>%
  column_to_rownames(var = "sample")

eub_samps <- eub %>%
  rownames_to_column(var = "sample") %>%
  left_join(sample_counts, by = "sample") %>%
  dplyr::select(sample, No..Individuals)

eub_wts <- eub_samps$No..Individuals

eub1 <- specaccum(eub)
eub2 <- specaccum(eub, "random", w = eub_wts)

plot(eub1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
boxplot(eub2, col="yellow", add=TRUE, pch="+")

#"Geophilomorpha"       "Euborellia annulipes" "Heteropoda venatoria"
##[4] "Keijia mneon"         "Neoscona theisi"      "Pantala flavescens"  
#[7] "Phisis holdhausi"     "Scytodes longipes"    "Smeringopus pallidus"

hev <- dna2 %>%
  filter(pred_ID == "Heteropoda venatoria") %>%
  dplyr::select(sample, unique_ID, reads) %>%
  group_by(sample, unique_ID) %>%
  summarise(reads = sum(reads)) %>%
  pivot_wider(names_from = "unique_ID",
              values_from = "reads") %>%
  column_to_rownames(var = "sample")

hev_samps <- hev %>%
  rownames_to_column(var = "sample") %>%
  left_join(sample_counts, by = "sample") %>%
  dplyr::select(sample, No..Individuals)
  
hev_wts <- hev_samps$No..Individuals

hev1 <- specaccum(hev)
hev2 <- specaccum(hev, "random", w = hev_wts)

plot(hev1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
boxplot(hev2, col="yellow", add=TRUE, pch="+")

#"Geophilomorpha"       "Euborellia annulipes" "Heteropoda venatoria"
##[4] "Keijia mneon"         "Neoscona theisi"      "Pantala flavescens"  
#[7] "Phisis holdhausi"     "Scytodes longipes"    "Smeringopus pallidus"

lrs <- dna2 %>%
  filter(pred_ID == "Keijia mneon") %>%
  dplyr::select(sample, unique_ID, reads) %>%
  group_by(sample, unique_ID) %>%
  summarise(reads = sum(reads)) %>%
  pivot_wider(names_from = "unique_ID",
              values_from = "reads") %>%
  column_to_rownames(var = "sample")

lrs_samps <- lrs %>%
  rownames_to_column(var = "sample") %>%
  left_join(sample_counts, by = "sample") %>%
  dplyr::select(sample, No..Individuals)

lrs_wts <- lrs_samps$No..Individuals

lrs1 <- specaccum(lrs)
lrs2 <- specaccum(lrs, "random", w = lrs_wts)

plot(lrs1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
boxplot(lrs2, col="yellow", add=TRUE, pch="+")

#"Geophilomorpha"       "Euborellia annulipes" "Heteropoda venatoria"
##[4] "Keijia mneon"         "Neoscona theisi"      "Pantala flavescens"  
#[7] "Phisis holdhausi"     "Scytodes longipes"    "Smeringopus pallidus"

neo <- dna2 %>%
  filter(pred_ID == "Neoscona theisi") %>%
  dplyr::select(sample, unique_ID, reads) %>%
  group_by(sample, unique_ID) %>%
  summarise(reads = sum(reads)) %>%
  pivot_wider(names_from = "unique_ID",
              values_from = "reads") %>%
  column_to_rownames(var = "sample")

neo_samps <- neo %>%
  rownames_to_column(var = "sample") %>%
  left_join(sample_counts, by = "sample") %>%
  dplyr::select(sample, No..Individuals)

neo_wts <- neo_samps$No..Individuals

neo1 <- specaccum(neo)
neo2 <- specaccum(neo, "random", w = neo_wts)

plot(neo1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
boxplot(neo2, col="yellow", add=TRUE, pch="+")

#"Geophilomorpha"       "Euborellia annulipes" "Heteropoda venatoria"
##[4] "Keijia mneon"         "Neoscona theisi"      "Pantala flavescens"  
#[7] "Phisis holdhausi"     "Scytodes longipes"    "Smeringopus pallidus"

pan <- dna2 %>%
  filter(pred_ID == "Pantala flavescens") %>%
  dplyr::select(sample, unique_ID, reads) %>%
  group_by(sample, unique_ID) %>%
  summarise(reads = sum(reads)) %>%
  pivot_wider(names_from = "unique_ID",
              values_from = "reads") %>%
  column_to_rownames(var = "sample")


pan1 <- specaccum(pan)
pan2 <- specaccum(pan, "random")

plot(pan1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
boxplot(pan2, col="yellow", add=TRUE, pch="+")

#########################
#Phisis holdausi
########################

phh <- dna2 %>%
  filter(pred_ID == "Phisis holdhausi") %>%
  dplyr::select(sample, unique_ID, reads) %>%
  group_by(sample, unique_ID) %>%
  summarise(reads = sum(reads)) %>%
  pivot_wider(names_from = "unique_ID",
              values_from = "reads") %>%
  column_to_rownames(var = "sample")

phh_samps <- phh %>%
  rownames_to_column(var = "sample") %>%
  left_join(sample_counts, by = "sample") %>%
  dplyr::select(sample, No..Individuals)

phh_wts <- phh_samps$No..Individuals

phh1 <- specaccum(phh)
phh2 <- specaccum(phh, "random", w=phh_wts)

plot(phh1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
boxplot(phh2, col="yellow", add=TRUE, pch="+")


#"Geophilomorpha"       "Euborellia annulipes" "Heteropoda venatoria"
##[4] "Keijia mneon"         "Neoscona theisi"      "Pantala flavescens"  
#[7] "Phisis holdhausi"     "Scytodes longipes"    "Smeringopus pallidus"

scy <- dna2 %>%
  filter(pred_ID == "Scytodes longipes") %>%
  dplyr::select(sample, unique_ID, reads) %>%
  group_by(sample, unique_ID) %>%
  summarise(reads = sum(reads)) %>%
  pivot_wider(names_from = "unique_ID",
              values_from = "reads") %>%
  column_to_rownames(var = "sample")

scy_samps <- scy %>%
  rownames_to_column(var = "sample") %>%
  left_join(sample_counts, by = "sample") %>%
  dplyr::select(sample, No..Individuals)

scy_wts <- scy_samps$No..Individuals

scy1 <- specaccum(scy)
scy2 <- specaccum(scy, "random", w=scy_wts)

plot(scy1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
boxplot(scy2, col="yellow", add=TRUE, pch="+")


#"Geophilomorpha"       "Euborellia annulipes" "Heteropoda venatoria"
##[4] "Keijia mneon"         "Neoscona theisi"      "Pantala flavescens"  
#[7] "Phisis holdhausi"     "Scytodes longipes"    "Smeringopus pallidus"

sme <- dna2 %>%
  filter(pred_ID == "Smeringopus pallidus") %>%
  dplyr::select(sample, unique_ID, reads) %>%
  group_by(sample, unique_ID) %>%
  summarise(reads = sum(reads)) %>%
  pivot_wider(names_from = "unique_ID",
              values_from = "reads") %>%
  column_to_rownames(var = "sample")

sme_samps <- sme %>%
  rownames_to_column(var = "sample") %>%
  left_join(sample_counts, by = "sample") %>%
  dplyr::select(sample, No..Individuals)

sme_wts <- sme_samps$No..Individuals

sme1 <- specaccum(sme)
sme2 <- specaccum(sme, "random", w=sme_wts)

plot(sme1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
boxplot(sme2, col="yellow", add=TRUE, pch="+")


####################
#Streamline####
###################

#extract predator 
sme <- dna2 %>%
  filter(pred_ID == "Smeringopus pallidus") %>%
  dplyr::select(sample, unique_ID, reads) %>%
  group_by(sample, unique_ID) %>%
  summarise(reads = sum(reads)) %>%
  pivot_wider(names_from = "unique_ID",
              values_from = "reads") %>%
  column_to_rownames(var = "sample")

sme_samps <- sme %>%
  rownames_to_column(var = "sample") %>%
  left_join(sample_counts, by = "sample") %>%
  dplyr::select(sample, No..Individuals)

sme_wts <- sme_samps$No..Individuals

sme1 <- specaccum(sme)
sme2 <- specaccum(sme, "random", w=sme_wts)

plot(sme1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
boxplot(sme2, col="yellow", add=TRUE, pch="+")