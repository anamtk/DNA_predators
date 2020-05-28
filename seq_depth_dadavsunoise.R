#This script requires both the dada2 and the UNOISE community data, as well as tidyverse and 
#ggplot, and maybe some others. Will need to update later because it's probs important

comm_long <- comm %>%
  gather(sample, reads, CEN10b:SMEb) %>%
  group_by(sample) %>%
  summarise(reads = sum(reads)) %>%
  mutate(approach = "together")

data_long <- data %>%
  gather(sample, reads, CEN01:SME14) %>%
  group_by(sample) %>%
  summarise(reads = sum(reads))

reads_dada <- as.data.frame(comm_long$reads)

reads_dada <- reads_dada %>%
  mutate(pipeline = "dada") %>%
  rename("reads" = "comm_long$reads")

reads_unoise <- as.data.frame(data_long$reads)

reads_unoise <- reads_unoise %>%
  mutate(pipeline = "unoise") %>%
  rename("reads" = "data_long$reads")

reads <- reads_dada %>%
  bind_rows(reads_unoise)

ggplot(reads, aes(x = reads, color = pipeline)) +
  geom_freqpoly(size = 2) +
  theme_bw() +
  scale_x_log10() +
  labs(x = "Total reads", y = "Number of samples")

c_run <- read.csv(here("data", "dada_may", "separate", "ASVs_counts_c.tsv"), sep = "\t")

#rename columns for simplicity
colnames(c_run) <- sapply(str_split(colnames(c_run), "_"), function(x){return(x[[1]])})
colnames(c_run) <- str_remove(colnames(c_run), "\\.")

c_run <- c_run %>%
  rename("ASV" = "X")

colnames(c_run) <- paste(colnames(c_run), "c", sep = "")

c_run_long <- c_run %>%
  dplyr::select(-CL12c, -CL42c, -NEGc, -QC1c) %>%
  gather(sample, reads, EUB10c:PAN9c) %>%
  group_by(sample) %>%
  summarise(reads = sum(reads)) %>%
  mutate(approach = "separate")

c_run_samples <- c_run_long$sample

comm_long_c <- comm_long %>%
  filter(sample %in% c_run_samples)

c_run_an <- c_run_long %>%
  bind_rows(comm_long_c)

ggplot(c_run_an, aes(x = reads, color = approach)) +
  geom_freqpoly(size = 2, alpha = 0.6) +
  theme_bw() +
  labs(x = "Total reads", y = "Number of samples")
