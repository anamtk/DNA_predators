###########################
#Controls check
#Ana Miller-ter Kuile
#March 3, 2021
###########################

#this is looks at controls in each run

###########################
#Load Packages ####
###########################

package.list <- c("here", "tidyverse")

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
                      "raw_data",
                      "1_denoised_data", 
                      "dada2", 
                      "b_combined_runs",
                      "ASVs_counts_all.tsv"), 
                 sep = "\t")

#rename columns for simplicity
colnames(comm) <- sapply(str_split(colnames(comm), "_"), function(x){return(x[[1]])})
colnames(comm) <- str_remove(colnames(comm), "\\.")

comm <- comm %>%
  rename("ASV" = "X")

###########################
#Subset controls and mutate ####
###########################

#select samples minus controls and the ASV column
controls <- comm %>%
  column_to_rownames(var = "ASV") %>%
  dplyr::select(CL12a, CL12b, CL12c, CL12d, CL42a, CL42b, CL42c,
                CL42d, NEGa, NEGb, NEGc, SMEb, QC1a, QC1b, QC1c,
                QC1d) 

controls <- controls %>%
  rownames_to_column(var = "ASV") %>%
  pivot_longer(CL12a:QC1d, 
               names_to = "sample",
               values_to = "reads") %>%
  filter(reads > 0)

negs <- controls %>%
  filter(sample %in% c("NEGa", "NEGb", "NEGc"))

negs %>%
  group_by(sample) %>%
  summarise(ASVs = n(),
            mean_reads = sum(reads),
            sd_reads = sd(reads)) %>%
  add_row(sample = "NEGd", ASVs = 0, mean_reads = 0, sd_reads = NA) %>%
  ggplot(aes(x = sample, y = mean_reads)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean_reads - sd_reads, ymax = mean_reads+sd_reads), 
                width = 0.2) +
  theme_bw()

negs %>%
    summarise(mean_reads = mean(reads),
              sd_reads = sd(reads),
              total = n(),
              se_reads=sd_reads/sqrt(total)) 

negs %>%
  group_by(sample) %>%
  summarise(sum_reads = sum(reads))

negs %>%
  summarise(mean_reads = sum(reads),
            sd_reads = sd(reads),
            total = n(),
            se_reads=sd_reads/sqrt(total)) 
  
pos <- controls %>%
  filter(!sample %in% c("NEGa", "NEGb", "NEGc", "SMEb"))

pos %>%
  group_by(sample) %>%
  summarise(ASVs = n(),
            mean_reads = sum(reads),
            sd_reads = sd(reads)) %>%
  ggplot(aes(x = sample, y = ASVs)) +
  geom_bar(stat = "identity") +
  theme_bw()

pos %>%
ggplot(aes(x = reads)) +
  geom_histogram() +
  facet_wrap(~sample) +
  theme_bw()

pos %>%
  filter(reads > 1656) %>%
  gt()
