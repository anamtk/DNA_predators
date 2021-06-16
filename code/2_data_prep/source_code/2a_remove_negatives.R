###########################
#Subtract negatives
#Ana Miller-ter Kuile
#March 17, 2021
###########################

# subtract negative ASVs from all other ASVs
# don't need to do for run D because the negative dropped out during filtering

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Load dataset ------------------------------------------------------------

#Dada2 combined runs ASVs, read counts by sample
comm <- read.csv(here("data", 
                      "raw_data",
                      "1_denoised_data", 
                      "dada2", 
                      "ASVs_counts_all.tsv"), 
                 sep = "\t")

#rename columns for simplicity
colnames(comm) <- sapply(str_split(colnames(comm), "_"), function(x){return(x[[1]])})
colnames(comm) <- str_remove(colnames(comm), "\\.")

comm <- comm %>%
  rename("ASV" = "X")

# Manipulate dataset ------------------------------------------------------
colnames(comm)

comm_long <- comm %>%
  pivot_longer(CEN10b:SMEb,
               names_to = "sample",
               values_to = "reads") %>%
  mutate(run = substr(sample, nchar(sample)-1+1, nchar(sample))) 

#NEGa
a <- comm_long %>%
  filter(run == "a") 

NEGa <- a %>%
  filter(sample == "NEGa") %>%
  rename(negative = reads) %>%
  dplyr::select(ASV, negative)

a <- a %>%
  left_join(NEGa, by = "ASV") %>%
  mutate(corrected_reads = reads-negative) %>%
  mutate(corrected_reads = ifelse(corrected_reads < 0, 0, corrected_reads))

#NEGb
b <- comm_long %>%
  filter(run == "b")

NEGb <- b %>%
  filter(sample == "NEGb") %>%
  rename(negative = reads) %>%
  dplyr::select(ASV, negative)

b <- b %>%
  left_join(NEGa, by = "ASV") %>%
  mutate(corrected_reads = reads-negative) %>%
  mutate(corrected_reads = ifelse(corrected_reads < 0, 0, corrected_reads))

#NEGc
c <- comm_long %>%
  filter(run == "c")

NEGc <- c %>%
  filter(sample == "NEGc") %>%
  rename(negative = reads) %>%
  dplyr::select(ASV, negative)

c <- c %>%
  left_join(NEGa, by = "ASV") %>%
  mutate(corrected_reads = reads-negative) %>%
  mutate(corrected_reads = ifelse(corrected_reads < 0, 0, corrected_reads))

#no negative because cleaned out
d <- comm_long %>%
  filter(run == "d")

# Re-combine --------------------------------------------------------------
 
cleaned <- a %>%
  bind_rows(b) %>%
  bind_rows(c) %>%
  dplyr::select(ASV, sample, corrected_reads) %>%
  rename(reads = corrected_reads) %>%
  bind_rows(d) %>%
  dplyr::select(-run) %>%
  pivot_wider(names_from = sample,
              values_from = reads)


# Export ------------------------------------------------------------------

write.csv(cleaned, (here("data", 
                       "outputs", 
                       "2a_remove_negatives", 
                       "ASVs_counts_all.csv")))

