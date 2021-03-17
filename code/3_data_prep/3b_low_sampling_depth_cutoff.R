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

package.list <- c("here", "tidyverse", "iNEXT")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

###########################
#Load data ####
###########################

#Dada2 combined runs ASVs, read counts by sample
comm <- read.csv(here("data", 
                      "outputs", 
                      "3a_remove_negatives", 
                      "ASVs_counts_all.csv"))

comm <- comm %>%
  dplyr::select(-X)

###########################
#Subset data for sequencing depth ####
###########################

#select samples minus controls and the ASV column
depth <- comm %>%
  column_to_rownames(var = "ASV") %>%
  dplyr::select(-CL12a, -CL12b, -CL12c, -CL12d, -CL42a, -CL42b, -CL42c,
                -CL42d, -SMEb, -QC1a, -QC1b, -QC1c,
                -QC1d, -NEGa, -NEGb, -NEGc)

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
#Quantiles to find change point in data below which sample removal is needed####
###########################
quants <- as.data.frame(quantile(sample_depth$n, c(.05, .06, .07, .08, 
                                                   .09, .10, .11, .12, .13, .14, 
                                                   .15, .16, .17, .18, .19, .2)))


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
#Subset samples for cross run comparisons####
###########################

#samples with sequencing depth below 11211 will be removed for community analyses
#I want DF to be exported in wide format, since this is how they will be rarefied
#in the next script
#subset samples for cross-run comparsions (none of these are below the
#sampling depth threshold)
cross <- depth %>%
  rownames_to_column(var = "ASV") %>% #take out row name
  dplyr::select("ASV", "HEV07a", "HEV10a", "HEV11a", "HEV12a", "HEV13a", 
                "HEV14a", "HEV15a", "HEV16a", "HEV17a", "HEV18a", "HEV20a", "HEV21a", "HEV22a",
                "HEV23a", "HEV24a", "HEV25a", "HEV26a", "HEV27a", "HEV29a",
                "HEV07b", "HEV10b", "HEV11b", "HEV12b", "HEV13b", "HEV14b", "HEV15b", "HEV16b",
                "HEV17b", "HEV18b", "HEV20b", "HEV21b", "HEV22b", "HEV23b", "HEV24b", "HEV25b",
                "HEV26b", "HEV27b", "HEV29b", 
                "HEV07c", "HEV10c", "HEV11c", "HEV12c", "HEV13c", "HEV14c", "HEV15c", "HEV16c",
                "HEV17c", "HEV18c", "HEV20c", "HEV21c", "HEV22c", "HEV23c", "HEV24c", "HEV25c",
                "HEV26c", "HEV27c", "HEV29c",
                "HEV07d", "HEV10d", "HEV11d", "HEV12d", "HEV13d", "HEV14d", "HEV15d", "HEV16d",
                "HEV17d", "HEV18d", "HEV20d", "HEV21d", "HEV22d", "HEV23d", "HEV24d", "HEV25d",
                "HEV26d", "HEV27d", "HEV29d") %>%
  gather(sample, reads, HEV07a:HEV29d) %>% #make long
  group_by(ASV) %>% #group by ASV so we can filter zero ASVs out from all
  filter(sum(reads) != 0) %>% #delete all ASVs that add to zero acros all samples
  pivot_wider(names_from = sample, #make wide again
              values_from = reads)
 
#Aside:
#these were on run D for the sterilization study, but I will not be including them
#duplicated for the community study. Currently they live no-where
#"HEV65d"  "HEV66d"  "HEV67d"  "HEV68d"  "HEV70d"  "HEV71d"  "HEV74d"  "HEV76d" 
#[341] "HEV79d"  "HEV81d"  "HEV82d"  "HEV83d"  "HEV87d"  "HEV88d"  "HEV89d" 

###########################
#Subset samples for community analyses####
###########################
#again, these will only include samples above our sequencing cut-off

#these are the samples which were sampled deeply enough based on
#our quantile cutoff
deep <- sample_depth %>%
  filter(n > 11211) %>%
  dplyr::select(site) #347 to 302

deep <- as.vector(deep$site)

#these are the data for community analyses, which still inculdes the 
#samples which were not sequenced deeply enough
all_data <- depth %>%
  rownames_to_column(var = "ASV") %>%
  dplyr::select("ASV", "HEV01a", "HEV02a", "HEV03a", "HEV04a", "HEV05a", "HEV07a", "HEV10a", "HEV11a",
                "HEV12a", "HEV13a", "HEV14a", "HEV15a", "HEV16a", "HEV17a", "HEV18a", "HEV20a",
                "HEV21a", "HEV22a", "HEV23a", "HEV24a", "HEV25a", "HEV26a", "HEV27a", "HEV29a",
                "HEV32a", "HEV34a", "HEV38a", "HEV39a", "HEV40a", "HEV42a", "HEV45a", "HEV49a",
                "HEV50a", "HEV65a", "HEV66a", "HEV67a", "HEV68a", "HEV70a", "HEV71a", "HEV74a", 
                "HEV76a", "HEV79a", "HEV81a", "HEV82a", "HEV83a", "HEV87a", "HEV88a", "HEV89a",
                "NEO10a", "NEO11a", "NEO12a", "NEO13a", "NEO14a", "NEO15a", "NEO16a", "NEO17a",
                "NEO18a", "NEO19a", "NEO1a", "NEO20a", "NEO21a", "NEO22a", "NEO23a", "NEO24a",
                "NEO25a", "NEO26a", "NEO2a", "NEO3a", "NEO4a", "NEO5a", "NEO6a", "NEO7a", "NEO8a",
                "NEO9a", "SCY10a", "SCY12a", "SCY15a", "SCY16a", "SCY17a", "SCY3a", "SCY4a", "SCY7a",
                "SCY8a", "SCY9a", "CEN10b", "CEN11b", "CEN12b", "CEN13b", "CEN14b", "CEN15b", 
                "CEN16b", "CEN17b", "CEN1b", "CEN2b", "CEN3b", "CEN4b", "CEN6b", "CEN7b", "CEN8b",
                "CEN9b", "PHH10b", "PHH16b", "PHH20b", "PHH21b", "PHH22b", "PHH23b", "PHH24b",
                "PHH25b", "PHH30b", "PHH32b", "PHH34b", "PHH35b", "PHH38b", "PHH39b", "PHH40b", 
                "PHH41b", "PHH42b", "PHH43b", "PHH44b", "PHH45b", "PHH46b", "PHH47b", "PHH49b",
                "PHH4b", "PHH50b", "PHH51b", "PHH52b", "PHH53b", "PHH54b", "PHH55b", "PHH56b",
                "PHH57b", "PHH58b", "PHH59b", "PHH5b", "PHH60b", "PHH61b", "PHH62b", "PHH63b", 
                "PHH65b", "PHH66b", "PHH67b", "PHH68b", "PHH69b", "PHH6b", "PHH70b", "PHH71b",
                "PHH72b", "PHH73b", "PHH74b", "PHH75b", "PHH7b", "PHH8b", "PHH9b", "SME10b", 
                "SME11b", "SME12b", "SME13b", "SME14b", "SME1b", "SME2b", "SME3b", "SME5b", 
                "SME6b", "SME7b", "SME8b", "SME9b", "EUB10c", "EUB12c", "EUB13c", "EUB15c",
                "EUB16c", "EUB18c", "EUB19c", "EUB1c", "EUB20c", "EUB21c", "EUB23c", "EUB24c",
                "EUB25c", "EUB26c", "EUB27c", "EUB28c", "EUB29c", "EUB2c", "EUB30c", "EUB31c",
                "EUB32c", "EUB34c", "EUB35c", "EUB36c", "EUB3c", "EUB4c", "EUB5c", "EUB6c",
                "EUB7c", "EUB8c", "EUB9c", "ISO10c", "ISO11c", "ISO12c", "ISO13c", "ISO14c",
                "ISO15c", "ISO16c", "ISO17c", "ISO18c", "ISO19c", "ISO1c", "ISO20c", "ISO21c",
                "ISO22c", "ISO23c", "ISO24c", "ISO25c", "ISO26c", "ISO27c", "ISO2c", "ISO30c",
                "ISO31c", "ISO32c", "ISO33c", "ISO34c", "ISO35c", "ISO36c", "ISO37c", "ISO4c",
                "ISO5c", "ISO6c", "ISO7c", "ISO8c", "ISO9c", "LRS2c", "LRS3c", "LRS5c", "LRS6c",
                "LRS7c", "LRS8c", "PAN10c", "PAN11c", "PAN12c", "PAN13c", "PAN14c", "PAN15c",
                "PAN1c", "PAN2c", "PAN3c", "PAN4c", "PAN5c", "PAN6c", "PAN7c", "PAN8c", 
                "PAN9c", "HEV100d", "HEV101d", "HEV102d", "HEV103d", "HEV104d", "HEV105d",
                "HEV106d", "HEV107d", "HEV108d","HEV109d", "HEV110d", "HEV111d", "HEV90d",
                "HEV91d", "HEV92d", "HEV93d", "HEV94d", "HEV95d", "HEV96d", "HEV97d", "HEV98d", "HEV99d") %>%
  gather(sample, reads, HEV01a:HEV99d) %>%
  group_by(ASV) %>% #group by ASV so we can filter zero ASVs out from all
  filter(sum(reads) != 0) #delete all ASVs that add to zero acros all samples
  
#40 samples had too-low sequencing, which included all ISO samples, so
#I'll be removing those samples from the dataframe for overall analyses
#I may still try to salvage data from the low samples, but we will seee

data <- all_data %>%
  filter(sample %in% deep) %>% #select only samples in the deep enough DF
  group_by(ASV) %>% #group by ASV so we can filter zero ASVs out from all
  filter(sum(reads) != 0)  %>%#delete all ASVs that add to zero acros all samples
  pivot_wider(names_from = sample, #make wide for rarefying
              values_from = reads)

#now we can make all_data wide as well for rarefying
all_data <- all_data %>%
  pivot_wider(names_from= sample,
              values_from = reads)

###########################
#Export for Rarefying####
###########################

#export wide format for rarefying 
write.csv(cross, here("data", 
                      "outputs", 
                      "3a_depth_corrected", 
                      "cross_run_samples.csv"))

write.csv(data, here("data", 
                     "outputs", 
                     "3a_depth_corrected", 
                     "depth_subset_samples.csv"))

