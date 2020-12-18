#this script lets you determine if you sampled each sample
#deeply enough. I did this for my test run to determine how many samples
#i could put on following runs. Specifically, I wanted to see how many reads minimum
#I would need per sample in order to get the full OTU diversity I observed
#This Requires that you convert something to a file format that vegan undersands
#and then you subsample the "diversity" (OTUs) and abundance (reads) of your dataset
#randomly to determine at what read number you reach the full diversity of OTUS

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", "gt", "vegan")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

#requires unlabeled rows that correspond to different samples, and columns that correspond to different species
#in my case, rows are samples, and columns correspond to different OTUs
data <- read.csv(here("data", 
                      "prelim_validation",
                      "all_raw_rarefaction.csv"), header = T)
str(data)
#remove unnecessary columns in the dataframe
data$OTU_ID <- NULL
data$X <- NULL
data <- t(data) #transpose data frame for analysis

S <- specnumber(data) #number of species per sample
S
raremax <- min(rowSums(data)) #minimum summed abundance of all species in a sample? 
raremax

#determines min sample size for each sample
#observing maximum species richness
rarecurve(data, step = 20, sample = raremax, col = "blue", cex = 0.6,
          xlab = "Number of reads", ylab = "OTU richness")

