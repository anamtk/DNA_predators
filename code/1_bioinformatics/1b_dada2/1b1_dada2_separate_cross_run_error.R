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


# Filter and Trim ---------------------------------------------------------


# Run A DADA2 ---------------------------------------------------

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

#plotQualityProfile(forward_reads)
#plotQualityProfile(reverse_reads)
# and just plotting the last 4 samples of the reverse reads
#plotQualityProfile(reverse_reads_1[17:20])
#plotQualityProfile(forward_reads_1[17:20])

#Happy Belly uses the info above to use the truncLen argument
#in the following filterAndTrim function. I am not going to do that just yet
#the original code looked like:
#filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
#                              reverse_reads, filtered_reverse_reads, maxEE=c(2,2),
#                              rm.phix=TRUE, minLen=175, truncLen=c(250,200))

#i'm going to do maxEE at 1, but may reconsider later 
#filtered_out1 <- filterAndTrim(forward_reads_1, 
#                               filtered_forward_reads_1,
#                               reverse_reads_1, 
#                               filtered_reverse_reads_1, 
#                               maxEE=c(1,1),
#                               rm.phix=TRUE, 
#                               minLen = 100, 
#                               multithread = TRUE,
#                               matchIDs = TRUE)

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

derep_forward_1 <- derepFastq(filtered_forward_reads_1, verbose=TRUE)
names(derep_forward_1) <- samples_a # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
derep_reverse_1 <- derepFastq(filtered_reverse_reads_1, verbose=TRUE)
names(derep_reverse_1) <- samples_a

dada_forward_1 <- dada(derep_forward_1, 
                       err=err_forward_reads_1, 
                       pool="pseudo", 
                       multithread = TRUE)
# dada_forward <- dada(derep_forward, err=err_forward_reads, pool="pseudo", multithread=TRUE) # problem running this way if on Binder
dada_reverse_1 <- dada(derep_reverse_1, 
                       err=err_reverse_reads_1, 
                       pool="pseudo", 
                       multithread = TRUE)
# dada_reverse <- dada(derep_reverse, err=err_reverse_reads, pool="pseudo", multithread=TRUE)

merged_amplicons_1 <- mergePairs(dada_forward_1, 
                                 derep_forward_1, 
                                 dada_reverse_1,
                                 derep_reverse_1, 
                                 trimOverhang=TRUE)

seqtab_1 <- makeSequenceTable(merged_amplicons_1)

seqtab.nochim_1 <- removeBimeraDenovo(seqtab_1, verbose=T) #273 bimeras of 844 input sequences

#loss of abundance from these 
sum(seqtab.nochim_1)/sum(seqtab_1) #0.9976155

getN <- function(x) sum(getUniques(x))

# making a little table
summary_tab_1 <- data.frame(row.names=samples_a, 
                            dada2_input=filtered_out1$reads.in,
                            filtered=filtered_out1$reads.out, 
                            dada_f=sapply(dada_forward_1, getN),
                            dada_r=sapply(dada_reverse_1, getN), 
                            merged=sapply(merged_amplicons_1, getN),
                            nonchim=rowSums(seqtab.nochim_1),
                            final_perc_reads_retained=round(rowSums(seqtab.nochim_1)/filtered_out1$reads.in*100, 
                                                            1))

summary_tab_1

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs_1 <- colnames(seqtab.nochim_1)
asv_headers_1 <- vector(dim(seqtab.nochim_1)[2], 
                        mode="character")

for (i in 1:dim(seqtab.nochim_1)[2]) {
  asv_headers_1[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta_1 <- c(rbind(asv_headers_1, asv_seqs_1))
write(asv_fasta_1, "ASVs_a.fa")

# count table:
asv_tab_1 <- t(seqtab.nochim_1)
row.names(asv_tab_1) <- sub(">", "", asv_headers_1)
write.table(asv_tab_1, "ASVs_counts_a.tsv", sep="\t", 
            quote=F, col.names=NA)


# Run B DADA2 ---------------------------------------------------

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

#plotQualityProfile(forward_reads)
#plotQualityProfile(reverse_reads)
# and just plotting the last 4 samples of the reverse reads
#plotQualityProfile(reverse_reads_2[17:20])
#plotQualityProfile(forward_reads_2[17:20])


#i'm going to do maxEE at 1, but may reconsider later 
#filtered_out2 <- filterAndTrim(forward_reads_2, 
#                               filtered_forward_reads_2,
#                               reverse_reads_2, 
#                               filtered_reverse_reads_2, 
#                               maxEE=c(1,1),
 #                              rm.phix=TRUE, 
  #                             minLen = 100, 
  #                             multithread = TRUE,
  #                             matchIDs = TRUE)

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

derep_forward_2 <- derepFastq(filtered_forward_reads_2, verbose=TRUE)
names(derep_forward_2) <- samples_b # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
derep_reverse_2 <- derepFastq(filtered_reverse_reads_2, verbose=TRUE)
names(derep_reverse_2) <- samples_b

dada_forward_2 <- dada(derep_forward_2, err=err_forward_reads_2, pool="pseudo", multithread = TRUE)
# dada_forward <- dada(derep_forward, err=err_forward_reads, pool="pseudo", multithread=TRUE) # problem running this way if on Binder
dada_reverse_2 <- dada(derep_reverse_2, err=err_reverse_reads_2, pool="pseudo", multithread = TRUE)
# dada_reverse <- dada(derep_reverse, err=err_reverse_reads, pool="pseudo", multithread=TRUE)

merged_amplicons_2 <- mergePairs(dada_forward_2, derep_forward_2, dada_reverse_2,
                                 derep_reverse_2, trimOverhang=TRUE)

seqtab_2 <- makeSequenceTable(merged_amplicons_2)

seqtab.nochim_2 <- removeBimeraDenovo(seqtab_2, verbose=T) #1023 bimeras of 1783 input sequences

#loss of abundance from these 
sum(seqtab.nochim_2)/sum(seqtab_2) #0.9735308

getN <- function(x) sum(getUniques(x))

# making a little table
summary_tab_2 <- data.frame(row.names=samples_b, 
                            dada2_input=filtered_out2$reads.in,
                            filtered=filtered_out2$reads.out, 
                            dada_f=sapply(dada_forward_2, getN),
                            dada_r=sapply(dada_reverse_2, getN), 
                            merged=sapply(merged_amplicons_2, getN),
                            nonchim=rowSums(seqtab.nochim_2),
                            final_perc_reads_retained=round(rowSums(seqtab.nochim_2)/filtered_out2$reads.in*100, 
                                                            1))

summary_tab_2

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs_2 <- colnames(seqtab.nochim_2)
asv_headers_2 <- vector(dim(seqtab.nochim_2)[2], mode="character")

for (i in 1:dim(seqtab.nochim_2)[2]) {
  asv_headers_2[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta_2 <- c(rbind(asv_headers_2, asv_seqs_2))
write(asv_fasta_2, "ASVs_b.fa")

# count table:
asv_tab_2 <- t(seqtab.nochim_2)
row.names(asv_tab_2) <- sub(">", "", asv_headers_2)
write.table(asv_tab_2, "ASVs_counts_b.tsv", sep="\t", quote=F, col.names=NA)

# Run C DADA2 ---------------------------------------------------

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

#plotQualityProfile(forward_reads)
#plotQualityProfile(reverse_reads)
# and just plotting the last 4 samples of the reverse reads
#plotQualityProfile(reverse_reads_2[17:20])
#plotQualityProfile(forward_reads_2[17:20])

#Happy Belly uses the info above to use the truncLen argument
#in the following filterAndTrim function. I am not going to do that just yet
#the original code looked like:
#filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
#                              reverse_reads, filtered_reverse_reads, maxEE=c(2,2),
#                              rm.phix=TRUE, minLen=175, truncLen=c(250,200))

#i'm going to do maxEE at 1, but may reconsider later 
#filtered_out3 <- filterAndTrim(forward_reads_3, 
#                               filtered_forward_reads_3,
#                               reverse_reads_3, 
#                               filtered_reverse_reads_3, 
#                               maxEE=c(1,1),
#                               rm.phix=TRUE, 
#                               minLen = 100, 
#                               multithread = TRUE,
#                               matchIDs = TRUE)


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

derep_forward_3 <- derepFastq(filtered_forward_reads_3, verbose=TRUE)
names(derep_forward_3) <- samples_c # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
derep_reverse_3 <- derepFastq(filtered_reverse_reads_3, verbose=TRUE)
names(derep_reverse_3) <- samples_c

dada_forward_3 <- dada(derep_forward_3, err=err_forward_reads_3, pool="pseudo", multithread = TRUE)
# dada_forward <- dada(derep_forward, err=err_forward_reads, pool="pseudo", multithread=TRUE) # problem running this way if on Binder

dada_reverse_3 <- dada(derep_reverse_3, err=err_reverse_reads_3, pool="pseudo", multithread = TRUE)
# dada_reverse <- dada(derep_reverse, err=err_reverse_reads, pool="pseudo", multithread=TRUE)

merged_amplicons_3 <- mergePairs(dada_forward_3, derep_forward_3, dada_reverse_3,
                                 derep_reverse_3, trimOverhang=TRUE)

seqtab_3 <- makeSequenceTable(merged_amplicons_3)

seqtab.nochim_3 <- removeBimeraDenovo(seqtab_3, verbose=T) #349 bimeras of 975 input sequences

#loss of abundance from these 
sum(seqtab.nochim_3)/sum(seqtab_3) #0.992098

getN <- function(x) sum(getUniques(x))

# making a little table
summary_tab_3 <- data.frame(row.names=samples_c, 
                            dada3_input=filtered_out3$reads.in,
                            filtered=filtered_out3$reads.out, 
                            dada_f=sapply(dada_forward_3, getN),
                            dada_r=sapply(dada_reverse_3, getN), 
                            merged=sapply(merged_amplicons_3, getN),
                            nonchim=rowSums(seqtab.nochim_3),
                            final_perc_reads_retained=round(rowSums(seqtab.nochim_3)/filtered_out3$reads.in*100, 
                                                            1))

summary_tab_3

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs_3 <- colnames(seqtab.nochim_3)
asv_headers_3 <- vector(dim(seqtab.nochim_3)[2], mode="character")

for (i in 1:dim(seqtab.nochim_3)[2]) {
  asv_headers_3[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta_3 <- c(rbind(asv_headers_3, asv_seqs_3))
write(asv_fasta_3, "ASVs_c.fa")

# count table:
asv_tab_3 <- t(seqtab.nochim_3)
row.names(asv_tab_3) <- sub(">", "", asv_headers_3)
write.table(asv_tab_3, "ASVs_counts_c.tsv", sep="\t", quote=F, col.names=NA)


# Run D DADA2 ---------------------------------------------------

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

#plotQualityProfile(forward_reads)
#plotQualityProfile(reverse_reads)
# and just plotting the last 4 samples of the reverse reads
#plotQualityProfile(reverse_reads_2[17:20])
#plotQualityProfile(forward_reads_2[17:20])

#Happy Belly uses the info above to use the truncLen argument
#in the following filterAndTrim function. I am not going to do that just yet
#the original code looked like:
#filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
#                              reverse_reads, filtered_reverse_reads, maxEE=c(2,2),
#                              rm.phix=TRUE, minLen=175, truncLen=c(250,200))

#i'm going to do maxEE at 1, but may reconsider later 
#filtered_out4 <- filterAndTrim(forward_reads_4, 
#                               filtered_forward_reads_4,
#                               reverse_reads_4, 
#                               filtered_reverse_reads_4,
#                               maxEE=c(1,1),
#                               rm.phix=TRUE, 
#                               minLen = 100, 
#                               multithread = TRUE,
#                               matchIDs = TRUE)

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

derep_forward_4 <- derepFastq(filtered_forward_reads_4b, verbose=TRUE)
names(derep_forward_4) <- samples_d2 # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
derep_reverse_4 <- derepFastq(filtered_reverse_reads_4b, verbose=TRUE)
names(derep_reverse_4) <- samples_d2

dada_forward_4 <- dada(derep_forward_4, 
                       err=err_forward_reads_4, 
                       pool="pseudo", 
                       multithread = TRUE)
# dada_forward <- dada(derep_forward, err=err_forward_reads, pool="pseudo", multithread=TRUE) # problem running this way if on Binder

dada_reverse_4 <- dada(derep_reverse_4, 
                       err=err_reverse_reads_4, 
                       pool="pseudo", 
                       multithread = TRUE)
# dada_reverse <- dada(derep_reverse, err=err_reverse_reads, pool="pseudo", multithread=TRUE)

merged_amplicons_4 <- mergePairs(dada_forward_4, 
                                 derep_forward_4, 
                                 dada_reverse_4,
                                 derep_reverse_4, 
                                 trimOverhang=TRUE)

seqtab_4 <- makeSequenceTable(merged_amplicons_4)

seqtab.nochim_4 <- removeBimeraDenovo(seqtab_4, verbose=T) #31 of 245 input sequences

#loss of abundance from these 
sum(seqtab.nochim_4)/sum(seqtab_4) #0.9994757

getN <- function(x) sum(getUniques(x))

# making a little table #not working! GAH
summary_tab_4 <- data.frame(row.names=samples_d2, 
                            dada3_input=filtered_out4b$reads.in,
                            filtered=filtered_out4b$reads.out, 
                            dada_f=sapply(dada_forward_4, getN),
                            dada_r=sapply(dada_reverse_4, getN), 
                            merged=sapply(merged_amplicons_4, getN),
                            nonchim=rowSums(seqtab.nochim_4),
                            final_perc_reads_retained=round(rowSums(seqtab.nochim_4)/filtered_out4b$reads.in*100, 1))

summary_tab_4

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs_4 <- colnames(seqtab.nochim_4)
asv_headers_4 <- vector(dim(seqtab.nochim_4)[2], mode="character")

for (i in 1:dim(seqtab.nochim_4)[2]) {
  asv_headers_4[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta_4 <- c(rbind(asv_headers_4, asv_seqs_4))
write(asv_fasta_4, "ASVs_d.fa")

# count table:
asv_tab_4 <- t(seqtab.nochim_4)
row.names(asv_tab_4) <- sub(">", "", asv_headers_4)
write.table(asv_tab_4, "ASVs_counts_d.tsv", sep="\t", quote=F, col.names=NA)


# Export Filtered Matrices -------------------------------------------------

write.csv(filtered_out4, "/Volumes/Passport/Dissertation/DNA_Stage_Structure_Chapter/data/outputs/filtered_out4.csv")
write.csv(filtered_out3, "/Volumes/Passport/Dissertation/DNA_Stage_Structure_Chapter/data/outputs/filtered_out3.csv")
write.csv(filtered_out2, "/Volumes/Passport/Dissertation/DNA_Stage_Structure_Chapter/data/outputs/filtered_out2.csv")
write.csv(filtered_out1, "/Volumes/Passport/Dissertation/DNA_Stage_Structure_Chapter/data/outputs/filtered_out1.csv")

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


