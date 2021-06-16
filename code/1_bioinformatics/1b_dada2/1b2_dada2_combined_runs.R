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

#All runs together ####

setwd(here("data",
           "raw_data",
           "0_raw_sequences",
           "e_combined_runs"))

samples_all <- scan("samples", what = "character")

# one holding the file names of all the forward reads
forward_reads_all <- paste0(samples_all, "_R1_trimmed.fq.gz")
# and one with the reverse
reverse_reads_all <- paste0(samples_all, "_R2_trimmed.fq.gz")

# and variables holding file names for the forward and reverse
# filtered reads we're going to generate below
filtered_forward_reads_all <- paste0(samples_all, "_R1_filtered.fq.gz")
filtered_reverse_reads_all <- paste0(samples_all, "_R2_filtered.fq.gz")

#i'm going to do maxEE at 1, but may reconsider later 
filtered_out_all <- filterAndTrim(forward_reads_all, filtered_forward_reads_all,
                                  reverse_reads_all, filtered_reverse_reads_all, maxEE=c(1,1),
                                  rm.phix=TRUE, minLen = 100, multithread = TRUE,
                                  matchIDs = TRUE)

filtered_out_all
#one of my negatives got filtered to zero so I need to delete that for the next
#step - it's throwing it off

#delete negative for the sequence table below
filtered_out_allb <- filtered_out_all %>%
  as_tibble(rownames = "sample") %>%
  filter(reads.out > 0) %>%
  as_matrix() #from the hciR package for this purpose!

#filtered_out %>%
#  as_tibble(rownames = "sample")

#need to find the samples with zero values in reads.out. This pipe does that:
empty <- filtered_out_all %>%
  as_tibble(rownames = "sample") %>%
  filter(reads.out <= 0) %>%
  select(sample) %>%
  as_vector() %>%
  unname() %>%
  str_sub(end=-18)

empty #which samples are empty 

#now remove them from forward and reverse, first by creating a vector with the name of that sample
empty_forward <- paste0(empty, "_R1_filtered.fq.gz")
empty_reverse <- paste0(empty, "_R2_filtered.fq.gz")
empty_reverse
#then removing from the filtered forward reads file
filtered_forward_reads_allb <- filtered_forward_reads_all[! filtered_forward_reads_all %in% empty_forward]

#and reverse reads
filtered_reverse_reads_allb <- filtered_reverse_reads_all[! filtered_reverse_reads_all %in% empty_reverse]

#and samples file
samples_all2 <- samples_all[! samples_all %in% empty]

#this would probably look better if I included the truncLen argument in the filtering
#step, however, I don't want to do that just yet. Also - My bp region is ~363 bp long,
#so would have to truncate >= to 182 - I feel like truncating at all may be a bad idea
#plotQualityProfile(filtered_reverse_reads_4[17:20])

#error model of forward and reverse reads
#the truncate length above may alter this error structure, something to consider
err_forward_reads_all <- learnErrors(filtered_forward_reads_allb, multithread = TRUE,
                                     randomize = TRUE)
err_reverse_reads_all <- learnErrors(filtered_reverse_reads_allb, multithread = TRUE,
                                     randomize = TRUE)

error_all_f <- plotErrors(err_forward_reads_all, nominalQ=TRUE)
error_all_r <- plotErrors(err_reverse_reads_all, nominalQ=TRUE)

derep_forward_all <- derepFastq(filtered_forward_reads_allb, verbose=TRUE)
names(derep_forward_all) <- samples_all2 # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
derep_reverse_all <- derepFastq(filtered_reverse_reads_allb, verbose=TRUE)
names(derep_reverse_all) <- samples_all2

#may need to use for the max integer warning. we'll see
#https://github.com/benjjneb/dada2/issues/432
#setDadaOpt( 
#  MIN_HAMMING=2,  
#  MIN_FOLD=2,
#  USE_QUALS=TRUE, 
#  BAND_SIZE=4,
#  VECTORIZED_ALIGNMENT=FALSE, 
#  USE_KMERS=TRUE,
#  MAX_CONSIST=200)      

dada_forward_all <- dada(derep_forward_all, err=err_forward_reads_all, pool="pseudo", multithread = TRUE)
# dada_forward <- dada(derep_forward, err=err_forward_reads, pool="pseudo", multithread=TRUE) # problem running this way if on Binder

dada_reverse_all <- dada(derep_reverse_all, err=err_reverse_reads_all, pool="pseudo", multithread = TRUE)
# dada_reverse <- dada(derep_reverse, err=err_reverse_reads, pool="pseudo", multithread=TRUE)

merged_amplicons_all <- mergePairs(dada_forward_all, derep_forward_all, dada_reverse_all,
                                   derep_reverse_all, trimOverhang=TRUE)

seqtab_all <- makeSequenceTable(merged_amplicons_all)

seqtab.nochim_all <- removeBimeraDenovo(seqtab_all, verbose=T) #1661 of 3399 input sequences

#loss of abundance from these 
sum(seqtab.nochim_all)/sum(seqtab_all) #0.990318

getN <- function(x) sum(getUniques(x))

# making a little table
summary_tab_all <- data.frame(row.names=samples_all2, dada3_input=filtered_out_allb[,1],
                              filtered=filtered_out_allb[,2], dada_f=sapply(dada_forward_all, getN),
                              dada_r=sapply(dada_reverse_all, getN), merged=sapply(merged_amplicons_all, getN),
                              nonchim=rowSums(seqtab.nochim_all),
                              final_perc_reads_retained=round(rowSums(seqtab.nochim_all)/filtered_out_allb[,1]*100, 1))

summary_tab_all

write_csv(summary_tab_all, "DADA2_sum_stats.csv")

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs_all <- colnames(seqtab.nochim_all)
asv_headers_all <- vector(dim(seqtab.nochim_all)[2], mode="character")

for (i in 1:dim(seqtab.nochim_all)[2]) {
  asv_headers_all[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta_all <- c(rbind(asv_headers_all, asv_seqs_all))
write(asv_fasta_all, here("raw_data",
                          "1_denoised_data",
                          "dada2",
                          "ASVs_all.fasta"))

# count table:
asv_tab_all <- t(seqtab.nochim_all)
row.names(asv_tab_all) <- sub(">", "", asv_headers_all)
write.table(asv_tab_all, here("raw_data",
                              "1_denoised_data",
                              "dada2",
                              "ASVs_counts_all.tsv"), sep="\t", quote=F, col.names=NA)

