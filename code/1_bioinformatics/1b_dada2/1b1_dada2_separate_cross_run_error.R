#Dada2 Error Rate Check Across Runs
#Ana Miller-ter Kuile
#February 3, 2020
#This was drawn from both Happy Belly Informatics and the github issue about 
#combining runs for dada2
#https://astrobiomike.github.io/amplicon/dada2_workflow_ex
#https://github.com/benjjneb/dada2/issues/257

#Installation Needs####
#First, make sure that Xcode (mac C++ functionality) was updated in this first
#in terminal: xcode-select --install

#installs BiocManager from BioConductor 
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.10")

#this code checks that all packages are up  to date
BiocManager::valid()

#install dada2
BiocManager::install("dada2")
BiocManager::install("tidybulk")

#load dada2 and get package version 
library(dada2); packageVersion("dada2")
library(tidybulk)
library(tidyverse) #for data cleaning and manipulation
#remotes::install_github("HuntsmanCancerInstitute/hciR")
#library(hciR) #as_matrix code
library(ShortRead)
library(here)
library(gridExtra)
library(patchwork)

#Checking Error of first run####
#from https://github.com/benjjneb/dada2/issues/257
#Additionally, you will want to check that the error models 
#from each run are giving basically the same result. 
#If the different runs have different error models, 
#we don't recommend pooling them together 
#(although you can as long as you use the highest 
#run-specific error model for the pool).

# Run A errors ---------------------------------------------------

#this is the samples object created by cutadapt - you may have to create a 
#new one in 
#terminal if you have moved folders
samples_a <- scan(here("data", "a_January_2019", 
                       "trimmed", "samples"), what = "character")

setwd(here("data", "a_January_2019", "trimmed"))

# one holding the file names of all the forward reads
forward_reads_1 <- paste0(samples_a, "_sub_R1_trimmed.fq.gz")
# and one with the reverse
reverse_reads_1 <- paste0(samples_a, "_sub_R2_trimmed.fq.gz")

# and variables holding file names for the forward and reverse
# filtered reads we're going to generate below
filtered_forward_reads_1 <- paste0(samples_a, "_sub_R1_filtered.fq.gz")
filtered_reverse_reads_1 <- paste0(samples_a, "_sub_R2_filtered.fq.gz")

filtered_out1 <- read.csv(here("data", "outputs", "filtered_out1.csv"))
filtered_out1

#error model of forward and reverse reads
#the truncate length above may alter this error structure, something to consider
err_forward_reads_1 <- learnErrors(filtered_forward_reads_1, randomize = TRUE,
                                   multithread = TRUE)
err_reverse_reads_1 <- learnErrors(filtered_reverse_reads_1, multithread = TRUE,
                                   randomize = TRUE)

error_1_f <- plotErrors(err_forward_reads_1, nominalQ=TRUE)
error_1_r <- plotErrors(err_reverse_reads_1, nominalQ=TRUE)

# Run B errors ---------------------------------------------------

setwd(here("data", "b_April_2019", "trimmed"))

samples_b <- scan("samples", what = "character")

# one holding the file names of all the forward reads
forward_reads_2 <- paste0(samples_b, "_sub_R1_trimmed.fq.gz")
# and one with the reverse
reverse_reads_2 <- paste0(samples_b, "_sub_R2_trimmed.fq.gz")

# and variables holding file names for the forward and reverse
# filtered reads we're going to generate below
filtered_forward_reads_2 <- paste0(samples_b, "_sub_R1_filtered.fq.gz")
filtered_reverse_reads_2 <- paste0(samples_b, "_sub_R2_filtered.fq.gz")

filtered_out2 <- read.csv(here("data", "outputs", "filtered_out2.csv"))

filtered_out2

#error model of forward and reverse reads
#the truncate length above may alter this error structure, something to consider
err_forward_reads_2 <- learnErrors(filtered_forward_reads_2, multithread = TRUE,
                                   randomize = TRUE)
err_reverse_reads_2 <- learnErrors(filtered_reverse_reads_2, multithread = TRUE,
                                   randomize = TRUE)

error_2_f <- plotErrors(err_forward_reads_2, nominalQ=TRUE)
error_2_r <- plotErrors(err_reverse_reads_2, nominalQ=TRUE)

# Run C errors ---------------------------------------------------

setwd(here("data", "c_September_2019", "trimmed"))

samples_c <- scan("samples", what = "character")

# one holding the file names of all the forward reads
forward_reads_3 <- paste0(samples_c, "_sub_R1_trimmed.fq.gz")
# and one with the reverse
reverse_reads_3 <- paste0(samples_c, "_sub_R2_trimmed.fq.gz")

# and variables holding file names for the forward and reverse
# filtered reads we're going to generate below
filtered_forward_reads_3 <- paste0(samples_c, "_sub_R1_filtered.fq.gz")
filtered_reverse_reads_3 <- paste0(samples_c, "_sub_R2_filtered.fq.gz")

filtered_out3 <- read.csv(here("data", "outputs", "filtered_out3.csv"))
filtered_out3

#error model of forward and reverse reads
#the truncate length above may alter this error structure, something to consider
err_forward_reads_3 <- learnErrors(filtered_forward_reads_3, multithread = TRUE,
                                   randomize = TRUE)
err_reverse_reads_3 <- learnErrors(filtered_reverse_reads_3, multithread = TRUE,
                                   randomize = TRUE)

error_3_f <- plotErrors(err_forward_reads_3, nominalQ=TRUE)
error_3_r <- plotErrors(err_reverse_reads_3, nominalQ=TRUE)

# Run D errors ---------------------------------------------------

setwd(here("data", "d_December_2019", "trimmed"))

samples_d <- scan("samples", what = "character")

# one holding the file names of all the forward reads
forward_reads_4 <- paste0(samples_d, "_sub_R1_trimmed.fq.gz")
# and one with the reverse
reverse_reads_4 <- paste0(samples_d, "_sub_R2_trimmed.fq.gz")

# and variables holding file names for the forward and reverse
# filtered reads we're going to generate below
filtered_forward_reads_4 <- paste0(samples_d, "_sub_R1_filtered.fq.gz")
filtered_reverse_reads_4 <- paste0(samples_d, "_sub_R2_filtered.fq.gz")


filtered_out4 <- read.csv(here("data", "outputs", "filtered_out4.csv"))

filtered_out4
#one of my negatives got filtered to zero so I need to delete that for the next
#step - it's throwing it off

#delete negative for the sequence table below
filtered_out4b <- filtered_out4 %>%
  filter(reads.out > 0) 

#filtered_out %>%
#  as_tibble(rownames = "sample")

#need to find the samples with zero values in reads.out. This pipe does that:
empty <- filtered_out4 %>%
  filter(reads.out <= 0) %>%
  select(X) %>%
  as_vector() %>%
  unname() %>%
  str_sub(end=-22)

empty #which samples are empty 

#now remove them from forward and reverse, first by creating a vector with the name of that sample
empty_forward <- paste0("NEG_S11_sub_R1_filtered.fq.gz")
empty_reverse <- paste0("NEG_S11_sub_R2_filtered.fq.gz")

#then removing from the filtered forward reads file
filtered_forward_reads_4b <- filtered_forward_reads_4[! filtered_forward_reads_4 %in% empty_forward]

#and reverse reads
filtered_reverse_reads_4b <- filtered_reverse_reads_4[! filtered_reverse_reads_4 %in% empty_reverse]

#and samples file
samples_d2 <- samples_d[! samples_d %in% empty]
samples_d2 <- samples_d[!samples_d == "NEG_S11"]
#this would probably look better if I included the truncLen argument in the filtering
#step, however, I don't want to do that just yet. Also - My bp region is ~363 bp long,
#so would have to truncate >= to 182 - I feel like truncating at all may be a bad idea
#plotQualityProfile(filtered_reverse_reads_4[17:20])

#error model of forward and reverse reads
#the truncate length above may alter this error structure, something to consider
err_forward_reads_4 <- learnErrors(filtered_forward_reads_4b, 
                                   multithread = TRUE,
                                   randomize = TRUE)
err_reverse_reads_4 <- learnErrors(filtered_reverse_reads_4b, 
                                   multithread = TRUE,
                                   randomize = TRUE)

error_4_f <- plotErrors(err_forward_reads_4, nominalQ=TRUE)
error_4_r <- plotErrors(err_reverse_reads_4, nominalQ=TRUE)

# Combine and export Error Data -------------------------------------------------------

er_1_f <- error_1_f$data
er_2_f <- error_2_f$data
er_3_f <- error_3_f$data
er_4_f <- error_4_f$data
er_1_r <- error_1_r$data
er_2_r <- error_2_r$data
er_3_r <- error_3_r$data
er_4_r <- error_4_r$data

er_1_f <- er_1_f %>%
  dplyr::select(Transition, Qual, count) %>%
  mutate(run = "One", reads = "Forward")

er_2_f <- er_2_f %>%
  dplyr::select(Transition, Qual, count) %>%
  mutate(run = "Two", reads = "Forward")

er_3_f <- er_3_f %>%
  dplyr::select(Transition, Qual, count) %>%
  mutate(run = "Three", reads = "Forward")

er_4_f <- er_4_f %>%
  dplyr::select(Transition, Qual, count) %>%
  mutate(run = "Four", reads = "Forward")

er_1_r <- er_1_r %>%
  dplyr::select(Transition, Qual, count) %>%
  mutate(run = "One", reads = "Reverse")

er_2_r <- er_2_r %>%
  dplyr::select(Transition, Qual, count) %>%
  mutate(run = "Two", reads = "Reverse")

er_3_r <- er_3_r %>%
  dplyr::select(Transition, Qual, count) %>%
  mutate(run = "Three", reads = "Reverse")

er_4_r <- er_4_r %>%
  dplyr::select(Transition, Qual, count) %>%
  mutate(run = "Four", reads = "Reverse")

errors_f <- er_1_f %>%
  bind_rows(er_2_f) %>%
  bind_rows(er_3_f) %>%
  bind_rows(er_4_f)
  
errors_r <- er_1_r %>%
  bind_rows(er_2_r) %>%
  bind_rows(er_3_r) %>%
  bind_rows(er_4_r)

write.csv(errors_f, here("data", "outputs", "forward_errors.csv"))
write.csv(errors_r, here("data", "outputs", "reverse_errors.csv"))


