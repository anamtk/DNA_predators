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

#######################
#OLD ANALYSIS: Subset runs to remove negative ####
#######################
colnames(comm)

a <- comm %>%
  dplyr::select("ASV", "NEGa", "CL12a", "CL42a", "HEV01a", "HEV02a", "HEV03a", "HEV04a",
                "HEV05a", "HEV07a", "HEV10a", "HEV11a", "HEV12a",
                "HEV13a", "HEV14a", "HEV15a", "HEV16a", "HEV17a", "HEV18a", 
                "HEV20a", "HEV21a", "HEV22a", "HEV23a", "HEV24a", "HEV25a", 
                "HEV26a", "HEV27a", "HEV29a", "HEV32a", "HEV34a", "HEV38a", "HEV39a", 
                "HEV40a", "HEV42a", "HEV45a", "HEV49a", "HEV50a", "HEV65a", "HEV66a",
                "HEV67a", "HEV68a", "HEV70a", "HEV71a", "HEV74a", "HEV76a", "HEV79a",
                "HEV81a", "HEV82a", "HEV83a", "HEV87a", "HEV88a", "HEV89a", 
                "NEO10a", "NEO11a", "NEO12a", "NEO13a", "NEO14a", "NEO15a", 
                "NEO16a", "NEO17a", "NEO18a", "NEO19a", "NEO1a",  
                "NEO20a", "NEO21a", "NEO22a", "NEO23a", "NEO24a", "NEO25a",
                "NEO26a", "NEO2a", "NEO3a", "NEO4a", "NEO5a", "NEO6a", "NEO7a",
                "NEO8a", "NEO9a", "QC1a", "SCY10a", "SCY12a", "SCY15a", "SCY16a", "SCY17a",
                "SCY3a", "SCY4a", "SCY7a", "SCY8a", "SCY9a") 

a[3:ncol(a)] <- a[3:ncol(a)]-a[,2]

a[a < 0] <- 0

write.csv(a, here("data", "outputs", "2_error_corrected", "run_a.csv"))
  
#this DF has two negatives, and the SMEb one I'm going to subtract from just the 
#SME samples? I don't know...
b <- comm %>%
  dplyr::select("ASV", "NEGb", "SMEb", "CEN10b", "CEN11b", "CEN12b", "CEN13b", "CEN14b", 
                "CEN15b", "CEN16b", "CEN17b", "CEN1b", "CEN2b", "CEN3b", "CEN4b", "CEN6b", 
                "CEN7b", "CEN8b", "CEN9b", "CL12b", "CL42b", "HEV07b", "HEV10b", "HEV11b", 
                "HEV12b", "HEV13b", "HEV14b", "HEV15b", "HEV16b", "HEV17b", "HEV18b", 
                "HEV20b","HEV21b", "HEV22b", "HEV23b", "HEV24b", "HEV25b", "HEV26b", 
                "HEV27b", "HEV29b", "PHH10b", "PHH16b", "PHH20b", "PHH21b", "PHH22b",
                "PHH23b", "PHH24b", "PHH25b", "PHH30b", "PHH32b", "PHH34b", "PHH35b", 
                "PHH38b", "PHH39b", "PHH40b", "PHH41b", "PHH42b", 
                "PHH43b", "PHH44b", "PHH45b", "PHH46b",
                "PHH47b", "PHH49b", "PHH4b", "PHH50b", "PHH51b", "PHH52b", "PHH53b",
                "PHH54b", "PHH55b", "PHH56b", "PHH57b", "PHH58b", "PHH59b", "PHH5b",
                "PHH60b", "PHH61b", "PHH62b", "PHH63b", "PHH65b", "PHH66b", "PHH67b", 
                "PHH68b", "PHH69b", "PHH6b", "PHH70b", "PHH71b", "PHH72b", "PHH73b",
                "PHH74b", "PHH75b", "PHH7b", "PHH8b", "PHH9b", "QC1b", "SME10b", 
                "SME11b", "SME12b", "SME13b", "SME14b", "SME1b", "SME2b", "SME3b", "SME5b",
                "SME6b", "SME7b", "SME8b", "SME9b") 

b1 <- b %>%
  dplyr::select("ASV", "NEGb", "CEN10b", "CEN11b", "CL12b", "CL42b", "QC1b", "CEN12b", 
                "CEN13b", "CEN14b", 
                "CEN15b", "CEN16b", "CEN17b", "CEN1b", "CEN2b", "CEN3b", "CEN4b", "CEN6b", 
                "CEN7b", "CEN8b", "CEN9b", "HEV07b", "HEV10b", "HEV11b", 
                "HEV12b", "HEV13b", "HEV14b", "HEV15b", "HEV16b", "HEV17b", "HEV18b", 
                "HEV20b","HEV21b", "HEV22b", "HEV23b", "HEV24b", "HEV25b", "HEV26b", 
                "HEV27b", "HEV29b", "PHH10b", "PHH16b", "PHH20b", "PHH21b", "PHH22b",
                "PHH23b", "PHH24b", "PHH25b", "PHH30b", "PHH32b", "PHH34b", "PHH35b", 
                "PHH38b", "PHH39b", "PHH40b", "PHH41b", "PHH42b", 
                "PHH43b", "PHH44b", "PHH45b", "PHH46b",
                "PHH47b", "PHH49b", "PHH4b", "PHH50b", "PHH51b", "PHH52b", "PHH53b",
                "PHH54b", "PHH55b", "PHH56b", "PHH57b", "PHH58b", "PHH59b", "PHH5b",
                "PHH60b", "PHH61b", "PHH62b", "PHH63b", "PHH65b", "PHH66b", "PHH67b", 
                "PHH68b", "PHH69b", "PHH6b", "PHH70b", "PHH71b", "PHH72b", "PHH73b",
                "PHH74b", "PHH75b", "PHH7b", "PHH8b", "PHH9b")

b1[3:ncol(b1)] <- b1[3:ncol(b1)]-b1[,2]

b1[b1 < 0] <- 0

write.csv(b1, here("data", "outputs", "2_error_corrected", "run_b1.csv"))

#this SME negative has relatively high reads... so I'm going to subset 
#just these since I think it might be contamination, but I'll have to wait until
#i get my lab notebook to verify.
b2 <- b %>%
  dplyr::select("ASV", "SMEb", "SME10b", 
                "SME11b", "SME12b", "SME13b", "SME14b", "SME1b", "SME2b", "SME3b", "SME5b",
                "SME6b", "SME7b", "SME8b", "SME9b")

b2[3:ncol(b2)] <- b2[3:ncol(b2)]-b2[,2]

b2[b2 < 0] <- 0

write.csv(b2, here("data", "outputs", "2_error_corrected", "run_b2.csv"))

b_all <- b1 %>%
  full_join(b2, by = "ASV")

write.csv(b_all, here("data", "outputs", "2_error_corrected", "run_b.csv"))

c <- comm %>%
  dplyr::select("ASV", "NEGc", "CL12c", "CL42c", "EUB10c", "EUB12c", "EUB13c", "EUB15c",
                "EUB16c", "EUB18c", "EUB19c", "EUB1c", "EUB20c", "EUB21c", "EUB23c", 
                "EUB24c", "EUB25c", "EUB26c", "EUB27c", "EUB28c", "EUB29c", "EUB2c",
                "EUB30c", "EUB31c", "EUB32c", "EUB34c", "EUB35c", "EUB36c", "EUB3c",
                "EUB4c", "EUB5c", "EUB6c", "EUB7c", "EUB8c", "EUB9c", "HEV07c", 
                "HEV10c", "HEV11c", "HEV12c", "HEV13c", "HEV14c", "HEV15c", "HEV16c", 
                "HEV17c", "HEV18c", "HEV20c", "HEV21c", "HEV22c", "HEV23c", "HEV24c", 
                "HEV25c", "HEV26c", "HEV27c", "HEV29c", "ISO10c", "ISO11c", "ISO12c",
                "ISO13c", "ISO14c", "ISO15c", "ISO16c", "ISO17c", "ISO18c", "ISO19c",
                "ISO1c", "ISO20c", "ISO21c", "ISO22c", "ISO23c", "ISO24c", "ISO25c",
                "ISO26c", "ISO27c", "ISO2c", "ISO30c", "ISO31c", "ISO32c", "ISO33c",
                "ISO34c", "ISO35c", "ISO36c", "ISO37c", "ISO4c", "ISO5c", "ISO6c",
                "ISO7c", "ISO8c", "ISO9c", "LRS2c", "LRS3c", "LRS5c", "LRS6c", "LRS7c",
                "LRS8c", "PAN10c", "PAN11c", "PAN12c", "PAN13c", "PAN14c",
                "PAN15c", "PAN1c", "PAN2c", "PAN3c", "PAN4c", "PAN5c", "PAN6c",  
                "PAN7c", "PAN8c", "PAN9c", "QC1c") 

c[3:ncol(c)] <- c[3:ncol(c)]-c[,2]

c[c < 0] <- 0

write.csv(c, here("data", "outputs", "2_error_corrected", "run_c.csv"))
#the negative disappeared from this set, so these are good to go!
d <- comm %>%
  dplyr::select("ASV", "CL12d", "CL42d", "HEV100d", "HEV101d", "HEV102d", "HEV103d",
                "HEV104d", "HEV105d", "HEV106d", "HEV107d", "HEV108d", "HEV109d",
                "HEV110d", "HEV111d", "HEV90d", "HEV91d", "HEV92d", "HEV93d", 
                "HEV94d", "HEV95d", "HEV96d", "HEV97d", "HEV98d", "HEV99d", "HEV07d",    
                "HEV10d", "HEV11d", "HEV12d", "HEV13d", "HEV14d", "HEV15d", "HEV16d",   
                "HEV17d", "HEV18d", "HEV20d", "HEV21d", "HEV22d", "HEV23d", "HEV24d",   
                "HEV25d", "HEV26d", "HEV27d", "HEV29d", "HEV65d", "HEV66d", "HEV67d",
                "HEV68d", "HEV70d", "HEV71d", "HEV74d", "HEV76d", "HEV79d", "HEV81d", 
                "HEV82d", "HEV83d", "HEV87d", "HEV88d", "HEV89d", "QC1d")

write.csv(d, here("data", "outputs", "2_error_corrected", "run_d.csv"))  

all_runs <- a %>%
  full_join(b_all, by = "ASV") %>%
  full_join(c, by = "ASV") %>%
  full_join(d, by = "ASV")

write.csv(all_runs, here("data", "outputs", "2_error_corrected", "all_runs.csv"))
