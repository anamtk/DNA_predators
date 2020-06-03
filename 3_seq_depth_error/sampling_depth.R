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
#library(vegan)
library(iNEXT)
#library(cowplot)
#library(ecodist)
#library(lme4)
#library(broom)
#library(MASS)
#library(ggeffects)
#library(glmmTMB)
#library(DHARMa)
#library(MuMIn)
#library(effects)

###########################
#Load data ####
###########################

#Dada2 combined runs ASVs, read counts by sample
comm <- read.csv(here("data", "denoised_data", "dada_may", "combined", "ASVs_counts_all.tsv"), sep = "\t")

#rename columns for simplicity
colnames(comm) <- sapply(str_split(colnames(comm), "_"), function(x){return(x[[1]])})
colnames(comm) <- str_remove(colnames(comm), "\\.")

#set row names to ASV labels
rownames(comm) <- comm$X

###########################
#Subset data for sequencing depth ####
###########################

#select samples minus controls and the ASV column
depth <- comm %>%
  dplyr::select(-X, -CL12a, -CL12b, -CL12c, -CL12d, -CL42a, -CL42b, -CL42c,
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

#graph the interpolated and extrapolated sampling depth per sample
ggiNEXT(seq_depth, type=1, facet.var="none", grey = T, se = F) + 
  # can set se = F to remove shaded regions to see lines better 
  theme_bw() +
  labs(x = "Sequencing Depth", y = "ASV Richness", title = "DADA2 Sequencing Depth") +
  theme(legend.position = "none", axis.text = element_text(size = 20), 
        axis.title = element_text(size = 25))

###########################
#Sequencing depth analysis by species####
###########################
#Maybe I should attach taxonomy first and then concatenate by taxonomies
#prior to sequencing depth analysis? 

#This will require also including taxonomies that are non-diet in the initial
#taxonomy DF... ok. or creating a new "non-diet" taxonomy DF, yes, that.


