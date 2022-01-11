#############################
#Comparing DADA2 to UNOISE3 ASV Assignments
#Ana Miller-ter Kuile
#June 3, 2020
#############################

#This script looks at ASV assignment and abundances between DADA2 and UNOISE3 to 
#determine which denoising method is best suited for these data (e.g. gives the
#most data for most samples)

#############################
#Load packages ####
#############################

package.list <- c("here", "tidyverse")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

#############################
#Load dataframes ####
#############################

#UNOISE
ASV1 <- read.delim(here("data",
                        "raw_data",
                        "1_denoised_data", 
                        "unoise", 
                        "zotu_table_a.txt"), 
                   sep = '\t')
ASV1$X.OTU.ID <- as.character(ASV1$X.OTU.ID)

ASV2 <- read.delim(here("data",
                        "raw_data",
                        "1_denoised_data", 
                        "unoise", 
                        "zotu_table_b.txt"), sep = '\t')
ASV2$X.OTU.ID <- as.character(ASV2$X.OTU.ID)

ASV3 <- read.delim(here("data",
                        "raw_data",
                        "1_denoised_data", 
                        "unoise", 
                        "zotu_table_c.txt"), sep = '\t')
ASV3$X.OTU.ID <- as.character(ASV3$X.OTU.ID)

ASV4 <- read.delim(here("data",
                        "raw_data",
                        "1_denoised_data", 
                        "unoise", 
                        "zotu_table_d.txt"), sep = '\t')
ASV4$X.OTU.ID <- as.character(ASV4$X.OTU.ID)

ASV5 <- read.delim(here("data",
                        "raw_data",
                        "1_denoised_data", 
                        "unoise", 
                        "zotu_table_e.txt"), sep = '\t')
ASV5$X.OTU.ID <- as.character(ASV5$X.OTU.ID)

ASV6 <- read.delim(here("data",
                        "raw_data",
                        "1_denoised_data", 
                        "unoise", 
                        "zotu_table_f.txt"), sep = '\t')
ASV6$X.OTU.ID <- as.character(ASV6$X.OTU.ID)

ASV7 <- read.delim(here("data",
                        "raw_data",
                        "1_denoised_data", 
                        "unoise", 
                        "zotu_table_g.txt"), sep = '\t')
ASV7$X.OTU.ID <- as.character(ASV7$X.OTU.ID)

ASV8 <- read.delim(here("data",
                        "raw_data",
                        "1_denoised_data", 
                        "unoise", 
                        "zotu_table_h.txt"), sep = '\t')
ASV8$X.OTU.ID <- as.character(ASV8$X.OTU.ID)

ASV9 <- read.delim(here("data",
                        "raw_data",
                        "1_denoised_data", 
                        "unoise", 
                        "zotu_table_j.txt"), sep = '\t')
ASV9$X.OTU.ID <- as.character(ASV9$X.OTU.ID)

ASV10 <- read.delim(here("data",
                         "raw_data",
                         "1_denoised_data", 
                         "unoise", 
                         "zotu_table_k.txt"), sep = '\t')
ASV10$X.OTU.ID <- as.character(ASV10$X.OTU.ID)

ASV11 <- read.delim(here("data",
                         "raw_data",
                         "1_denoised_data", 
                         "unoise", 
                         "zotu_table_l.txt"), sep = '\t')
ASV11$X.OTU.ID <- as.character(ASV11$X.OTU.ID)

u3 <- ASV1 %>%
  left_join(ASV2, by = "X.OTU.ID") %>%
  left_join(ASV3, by = "X.OTU.ID") %>%
  left_join(ASV4, by = "X.OTU.ID") %>%
  left_join(ASV5, by = "X.OTU.ID") %>%
  left_join(ASV6, by = "X.OTU.ID") %>%
  left_join(ASV7, by = "X.OTU.ID") %>%
  left_join(ASV8, by = "X.OTU.ID") %>%
  left_join(ASV9, by = "X.OTU.ID") %>%
  left_join(ASV10, by = "X.OTU.ID") %>%
  left_join(ASV11, by = "X.OTU.ID") %>%
  rename("ASV" = "X.OTU.ID") %>%
  replace(., is.na(.), 0)

#DADA2
d2 <- read.csv(here("data", 
                    "raw_data",
                    "1_denoised_data",
                    "dada2",
                    "ASVs_counts_all.tsv"), sep = "\t")

#rename columns for simplicity
colnames(d2) <- sapply(str_split(colnames(d2), "_"), function(x){return(x[[1]])})
colnames(d2) <- str_remove(colnames(d2), "\\.")

d2 <- d2 %>%
  rename("ASV" = "X")

#############################
#Extract per sample read abundance from each####
#############################

d2_long <- d2 %>%
  gather(sample, reads, CEN10b:SMEb) %>% #make long
  group_by(sample) %>% #group by sample
  summarise(reads = sum(reads)) %>% #summarize total reads per sample
  mutate(approach = "together") #make a column that indicates these were D2 of all samples

u3_long <- u3 %>%
  gather(sample, reads, CEN01:SME14) %>% #make long
  group_by(sample) %>% #group by sample
  summarise(reads = sum(reads)) #summarize total reads per sample

reads_dada <- as.data.frame(d2_long$reads) #extract reads from D2 DF

reads_dada <- reads_dada %>% #create df of reads
  mutate(pipeline = "dada") %>% #add pipeline for comparison
  rename("reads" = "d2_long$reads") #rename for brevity and combining with U3

reads_unoise <- as.data.frame(u3_long$reads) #extract reads from U3 DF

reads_unoise <- reads_unoise %>% #create DF of reads
  mutate(pipeline = "unoise") %>% #add pipeline for comparison
  rename("reads" = "u3_long$reads") #rename for brevity and combining with D2

reads <- reads_dada %>% #combine reads per sample for dada2 and unoise3
  bind_rows(reads_unoise)

#compare these visaully
ggplot(reads, aes(x = reads, color = pipeline)) + 
  geom_freqpoly(size = 2) +
  theme_bw() +
  scale_x_log10() +
  labs(x = "Total reads", y = "Number of samples")

#DADA2 has fewer low-abundance samples than UNOISE 3, so I will be using the 
#DADA2 data for analyses in this project

