#######################
#Negative Sample Error correcting
#Ana Miller-ter Kuile
#June 3, 2020
#######################

#I'm fighting with some low sequencing depth values, and I'm wondering if removing
#the error from sequencing (based on the negative sample ASV abundances) prior
#to doing any analyses of sequencing depth is not the best way to go. 
#To this end, I'm going to copy what Austen does for his (since I'm not convinced
#of the Jerde method). Specifically, more abundant ASVs occur in negatives at greater
#abundances, and so I'm going to just erase this number of sequences across all the 
#samples that were a part of that sequencing run. 

#######################
#Load packages ####
#######################
library(here)
library(tidyverse)
library(ggplot2)

#######################
#Load data ####
#######################

#Dada2 combined runs ASVs, read counts by sample
comm <- read.csv(here("data", "denoised_data", "dada_may", "combined", "ASVs_counts_all.tsv"), sep = "\t")

#rename columns for simplicity
colnames(comm) <- sapply(str_split(colnames(comm), "_"), function(x){return(x[[1]])})
colnames(comm) <- str_remove(colnames(comm), "\\.")

comm <- comm %>%
  rename("ASV" = "X")

#######################
#Assess control ASVs ####
#######################
#look at ASVs that were assigned to positive controls, so we can
#look at whether the community assigned any non-controls to these ASVs
positive <- comm %>%
  dplyr::select(ASV, CL12a, CL12b, CL12c, CL12d, CL42a, CL42b, CL42c, CL42d, QC1a, 
                QC1b, QC1c, QC1d)  %>%
  gather(sample, reads, CL12a:QC1d) %>%
  filter(reads > 0) %>%
  dplyr::select(ASV) %>%
  unique() %>%
  as_vector()

#this makes me think I need to be careful about the SME samples, but 
#we'll see
negative <- comm %>%
  dplyr::select(ASV, NEGa, NEGb, NEGc, SMEb) %>%
  gather(sample, reads, NEGa:SMEb) %>%
  filter(reads > 0)

#######################
#Subset controls ASVs from community####
#######################
#looks like positive control ASVs were not 
#assigned to any samples. Austen says this
#is sufficient correction check?

pos <- comm %>%
  filter(ASV %in% c(positive))

