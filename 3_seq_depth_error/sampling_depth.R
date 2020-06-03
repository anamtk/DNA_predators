###########################
#Sampling Depth
#Ana Miller-ter Kuile
#June 3, 2020
###########################

#this is code for looking at sampling depth across samples to ensure that I've
#sufficiently sampled each sample in this dataset. 
#this code looks at sampling depth in community matrices created by dada2 

###########################
#Load Packages ####
###########################
library(tidyverse)
library(here)
library(ggplot2)
#library(ggfortify)
#library(GUniFrac)
library(vegan)
library(iNEXT)
#library(cowplot)
#library(ecodist)
#library(lme4)
#library(broom)
#library(MASS)
#library(ggeffects)
library(glmmTMB)
library(DHARMa)
library(MuMIn)
#library(effects)

###########################
#Load data ####
###########################

#Dada2 combined runs ASVs, read counts by sample
comm <- read.csv(here("data", "outputs", "2_error_corrected", "all_runs.csv"))

#set row names to ASV labels
rownames(comm) <- comm$ASV

###########################
#Subset data for sequencing depth ####
###########################

#select samples minus controls and the ASV column
depth <- comm %>%
  dplyr::select(-X, -ASV, -CL12a, -CL12b, -CL12c, -CL12d, -CL42a, -CL42b, -CL42c,
                -CL42d, -NEGa, -NEGb, -NEGc, -SMEb, -QC1a, -QC1b, -QC1c,
                -QC1d)

#remove any ASVs that have zeros across all samples (probably from removing control)
depth <- depth[rowSums(depth) != 0,] #8 ASVs disappeared

###########################
#Sequencing depth analysis ####
###########################
#Maybe I should attach taxonomy first and then concatenate by taxonomies
#prior to sequencing depth analysis? Maybe will try below. 
#this determines sequencing depth across all samples
seq_depth <- iNEXT(depth, q=0, datatype="abundance") #this determines sequencing depth for each sample

seq_depth$DataInfo$SC
sample_depth <- seq_depth$DataInfo

quantile(sample_depth$n)
boxplot.stats(sample_depth$n)$out

#graph the interpolated and extrapolated sampling depth per sample
ggiNEXT(seq_depth, type=1, facet.var="none", grey = T, se = F) + 
  # can set se = F to remove shaded regions to see lines better 
  theme_bw() +
  labs(x = "Sequencing Depth", y = "ASV Richness", title = "DADA2 Sequencing Depth") +
  theme(legend.position = "none", axis.text = element_text(size = 20), 
        axis.title = element_text(size = 25))

###########################
#Quantiles to find change point in data below which sample removal is needed####
###########################
quants <- as.data.frame(quantile(sample_depth$n, c(.05, .06, .07, .08, 
                                                   .09, .10, .11, .12, .13, .14, 
                                                   .15, .16, .17, .18, .19, .2)))

colnames(quants)
quants <- quants %>%
  rename("reads" = "quantile(sample_depth$n, c(0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2))") %>%
  rownames_to_column(var = "quantile") %>%
  mutate(number = 1:n()) %>%
  mutate(number = number +4)

ggplot(quants, aes(x = number, y = reads)) +
  geom_point() +theme_bw() +
  geom_hline(yintercept = 11211, linetype = "dashed") +
  labs(y = "Sequence reads", x = "Quantile")

###########################
#Sequencing depth analysis by species####
###########################

#this was a previous version of exploration, I feel okay about the 
#quantile cutoff above, so not going to explore this anymore

#I wonder if this would be different if I concatenated the shared taxonomies together
#First. The few samples which have really low read abundance, according to the above
#analysis, are sequenced well (which I don't believe), so I am thinking it might
#be judicious to concatenate shared taxonomic IDs/categories prior to running this
#sequencing depth analysis to see if the cutoff changes. 

###########################
#Attach and Concatenate taxonomies ####
###########################
#Maybe I should attach taxonomy first and then concatenate by taxonomies
#prior to sequencing depth analysis? 

taxonomies <- read.csv(here("data", "outputs", "1_taxonomic_assignment", 
                            "ASV_taxonomies.csv"))

#summarise target DNA
target_comm <- comm %>%
  dplyr::select(-X) %>%
  gather(sample, reads, NEGa:QC1d) %>%
  left_join(taxonomies, by = "ASV") %>%
  group_by(sample, unique_ID) %>%
  summarise(reads = sum(reads)) %>%
  filter(!sample %in% c("CL12a", "CL12b", "CL12c", "CL12d", "CL42a", "CL42b", "CL42c",
                        "CL42d", "NEGa", "NEGb", "NEGc", "QC1a", "QC1b", "QC1c", "QC1d",
                        "SMEb")) %>%
  pivot_wider(names_from = sample, 
              values_from = reads) 

target_comm$unique_ID[is.na(target_comm$unique_ID)] <- "non-target"
  
target_comm <- target_comm %>%
  column_to_rownames(var = "unique_ID")

###########################
#Seq depth on concatenated DNA sequences####
###########################

seq_depth_cat <- iNEXT(target_comm, q=0, datatype="abundance") #this determines sequencing depth for each sample

seq_depth_cat$DataInfo$SC
sample_depth_cat <- seq_depth_cat$DataInfo

#graph the interpolated and extrapolated sampling depth per sample
ggiNEXT(seq_depth_cat, type=1, facet.var="none", grey = T, se = F) + 
  # can set se = F to remove shaded regions to see lines better 
  theme_bw() +
  labs(x = "Sequencing Depth", y = "ASV Richness", title = "DADA2 Sequencing Depth") +
  theme(legend.position = "none", axis.text = element_text(size = 20), 
        axis.title = element_text(size = 25))

###########################
#Non target considerations below####
###########################
tax_categories <- read.csv(here("data", "outputs", "1_taxonomic_assignment", 
                                "all_ASV_tax.csv"))


non_target <- comm %>%
  rename("ASV" = "X") %>%
  gather(sample, reads, CEN10b:SMEb) %>%
  left_join(tax_categories, by = "ASV") %>%
  filter(taxonomy != "target") %>%
  mutate(unique_ID = "non-target")


