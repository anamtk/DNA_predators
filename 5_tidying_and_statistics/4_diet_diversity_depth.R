###########################
#Diet Diversity and Depth
#Ana Miller-ter Kuile
#June 9, 2020
###########################

#This code will look at how diet diversity by predator is represented
#in our dataset

#How much diet diversity did we capture per predator individual?

#Did we capture the majority of the diet of our predators?
##########################
#Load packages####
###########################
library(here)
library(tidyverse)
library(ggplot2)
library(iNEXT)

#########################
#Load data####
###########################
#data from the taxonomic sort, which is already concatenated
#by unique ID within a predator individual
dna <- read.csv(here("data", "outputs", "5_rarefied_taxonomic_sort", "prey_DNA.csv"))

#sample metadata for sample sizes
meta <- read.csv(here("data", "Sample_metadata.csv"))

##########################
#Tidy data for analyses####
###########################
#Per predator diet for visualization
indivs <- dna %>%
  filter(reads > 0) %>% #remove zero reads 
  group_by(sample, pred_ID) %>% #group by sample and keep predator species
  tally(name = "species") %>% #count the number of diet items by species
  mutate(run = substr(sample, nchar(sample)-1+1, nchar(sample))) #indicates the 
#run it was run on

indivs <- indivs %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>% #join to metadata
  distinct(sample, pred_ID, species, run, Method, Island, Habitat, Microhabitat,
           Year, Date.Collected, No..Individuals) #keep only distinct values

#Presence of DNA of each type in each predator species
dna$sample <- str_sub(dna$sample, start = 1, end =-2) #remove "a" at end

#create a DF that summarises the population-level incidience data for 
#each diet item in each predator species
pred_pres <- dna %>%
  group_by(pred_ID, unique_ID) %>% #group each predator's diet by species
  summarise(presence = sum(reads > 0)) %>% #summarize the total number of non-zero reads for each species in pred population
  pivot_wider(names_from = pred_ID, #make wide for iNext making pred_ID the columns and diet species rows
              values_from = presence) %>%
  column_to_rownames(var = "unique_ID") #set the column to rownames for analysis

pred_pres[is.na(pred_pres)] <- 0 #set NA to zero
pred_pres <- pred_pres[rowSums(pred_pres) != 0,] #remove zero sum diet species, removed 16

#create a df that counts the number of individuals of each species, needed
#as first column of the DF for iNEXT with incidence data
species_counts <- dna %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>% #join metadata
  distinct(ID, sample, No..Individuals) %>% #select only the distinct samples
  group_by(ID) %>%
  summarise(count = sum(No..Individuals)) #count the number of individauls of species

#transpose the DF
sp_count <- as.data.frame(t(species_counts))

colnames(sp_count) <- sp_count[1, ] #make the column names the predator species names
sp_count <- sp_count[-1,] #remove the predator species column
row.names(sp_count) <- "total" #set the rowname to "total"

#add the species count column (df) to top of the presence df for iNEXT structure
pres <- sp_count %>%
  bind_rows(mutate_all(pred_pres, as.character))  #had to convert to character

rows <- c("total", rownames(pred_pres)) #rownames disappeared, so add them back in
rownames(pres) <- rows #add row names back in

#convert characters back to numeric for iNEXT
pres <- pres %>%
  mutate_if(is.character, as.numeric)
##########################
#Per individual diet by species####
###########################
#get sum stats
indivs %>%
  ungroup() %>%
  summarise(mean = mean(species), total = n(), sd = sd(species), se= sd/sqrt(total))

#plot per individual diet richness by species
ggplot(indivs, aes(x = pred_ID, y = species/No..Individuals)) +
  geom_hline(yintercept = 3.84, linetype = "dashed") +
  geom_hline(yintercept = 3.84 -0.139, color = "grey") +
  geom_hline(yintercept = 3.84  +0.139, color = "grey") +
  geom_boxplot() +theme_bw() 

##########################
#iNEXT####
###########################
#sampling depth based on incidence data at population level
pres_depth <- iNEXT(pres, q=0, datatype="incidence_freq")

#examine sampling completeness
pres_depth$DataInfo$SC
#look at diversity indices and estimates (may want to extract later):
pres_depth$AsyEst

#graph the interpolated and extrapolated species richness per species
ggiNEXT(pres_depth, type=1, facet.var="none", se=FALSE) + 
  # can set se = F to remove shaded regions to see lines better 
  theme_bw() +
  labs(x = "Number of Individuals", y = "Diet Richness") +
  theme(axis.text = element_text(size = 20), 
        axis.title = element_text(size = 25))

##########################
#SUPPLEMENT: vegan, by species and clunky####
###########################

sample_counts <- dna %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  distinct(ID, sample, No..Individuals) %>%
  dplyr::select(ID, sample, No..Individuals)

##########################
#vegan: geophilomorpha ####
###########################

geo <- dna %>%
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

##########################
#vegan: Euborellia annulipes ####
###########################
eub <- dna %>%
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

##########################
#vegan: "Heteropoda venatoria" ####
###########################
hev <- dna %>%
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

##########################
#vegan: Keijia mneon ####
###########################

lrs <- dna %>%
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

##########################
#vegan: Neoscona theisi ####
###########################
neo <- dna %>%
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

##########################
#vegan: Pantala flavescens ####
###########################

pan <- dna %>%
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

phh <- dna %>%
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

##########################
#vegan: Scytodes longipes ####
###########################

scy <- dna %>%
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

##########################
#vegan: Smeringopus pallidus ####
###########################
sme <- dna %>%
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
