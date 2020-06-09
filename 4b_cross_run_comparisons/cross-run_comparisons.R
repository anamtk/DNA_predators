###########################
#Cross-run comparisons
#Ana Miller-ter Kuile
#June 4, 2020
###########################

#This script looks at how ASV assignments and presence varied from run to
#run for a set of samples I ran on every single sequencing run that is part
#of this study. 
#I'll be treating this like a repeated measures ANOVA in mixed model format
#specifically setting up a model that asks, are these distributions different
#across different runs given a grouping by the ASVs that were assigned to
#each sample across each run

#Does each run vary in the presence of each ASV assigned to each sample?

#presence ~ run + (1|Sample/ASV)

###########################
#Load Packages####
###########################
library(here)
library(tidyverse)
library(ggplot2)
library(glmmTMB)
library(DHARMa)
library(effects)
library(emmeans)
library(MuMIn)
library(lattice)
library(performance)

###########################
#Load Data####
###########################

#rarefied cross-run samples
cross <- read.csv(here("data", "outputs", "4_rarefied", "cross_run_rare.csv"))

taxa <- read.csv(here("data", "outputs", "1_taxonomic_assignment", "ASV_taxonomies.csv"))
###########################
#Manipulate DF to long for analyses####
###########################

#Divide the samples up by sequencing run, give them a run category to distinguuish
#them, and then remove the a, b, c, d notation for GLMM below

#Run A
run_a <- cross %>%
  rename("ASV" = "X") %>% #rename X to ASV
  dplyr::select(ASV, HEV07a:HEV29a) %>% #select run 1
  gather(sample, reads, HEV07a:HEV29a) %>% #make this long
  mutate(run = "A") #give a run specific ID

run_a$sample <- str_sub(run_a$sample, start = 1, end =-2) #remove "a" at end

#Run B
run_b <- cross %>%
  rename("ASV" = "X") %>%
  dplyr::select(ASV, HEV07b:HEV29b) %>%
  gather(sample, reads, HEV07b:HEV29b) %>%
  mutate(run = "B")

run_b$sample <- str_sub(run_b$sample, start = 1, end =-2)

#Run C
run_c <- cross %>%
  rename("ASV" = "X") %>%
  dplyr::select(ASV, HEV07c:HEV29c) %>%
  gather(sample, reads, HEV07c:HEV29c) %>%
  mutate(run = "C")

run_c$sample <- str_sub(run_c$sample, start = 1, end =-2)

#Run D
run_d <- cross %>%
  rename("ASV" = "X") %>%
  dplyr::select(ASV, HEV07d:HEV29d) %>%
  gather(sample, reads, HEV07d:HEV29d) %>%
  mutate(run = "D")
  
run_d$sample <- str_sub(run_d$sample, start = 1, end =-2)

#bind them by rows for a long DF
all_run <- run_a %>%
  bind_rows(run_b) %>%
  bind_rows(run_c) %>%
  bind_rows(run_d) %>%
  group_by(ASV) %>%
  filter(sum(reads) > 0) %>%
  mutate(presence = ifelse(reads > 0,1,0))
# from 15580 - 13452 = 2128, or 28 ASVs?

#Factor varibles
all_run$run <- factor(all_run$run, levels = c("A", "B", "C", "D"))
all_run$ASV <- as.factor(all_run$ASV)
all_run$sample <- as.factor(all_run$sample)

###########################
#Cross-run comparison of ASV composition by sample####
###########################

#going to do presence-absence composition analysis

mod1 <- glmmTMB(presence ~ run + (1|sample) + (1|ASV),
                data = all_run,
                family = 'binomial')

mod_null <- glmmTMB(presence ~ 1 + (1|sample) + (1|ASV),
                data = all_run,
                family = 'binomial')

AICc(mod1, mod_null)

mod1 <- glmmTMB(presence ~ run + (1|sample) + (1|ASV),
                    data = all_run,
                    family = 'binomial')

binned_residuals(mod1) #looks not great, but not terrible
simulationOutput <- simulateResiduals(fittedModel = mod1) #ok
fit <- plot(simulationOutput, asFactor=TRUE) #ok

plot(allEffects(mod1)) #b is higher, C and D lower
em <- emmeans(mod1, "run")
pairs(em) #B-C and B-D different

###########################
#total ASVs per sample by run####
###########################
#non-zero ASVs per sample
num_ASV <- all_run %>%
  group_by(run, sample) %>%
  filter(reads > 0) %>%
  tally(name = "ASVs")

ASV_mod <- glmmTMB(ASVs ~ run + (1|sample),
                   data=num_ASV,
                   family = "genpois",
                   REML = FALSE)

ASV_null <- glmmTMB(ASVs ~ 1 + (1|sample),
                   data=num_ASV,
                   family = "genpois",
                   REML = FALSE)

AICc(ASV_mod, ASV_null)

ASV_mod <- glmmTMB(ASVs ~ run + (1|sample),
                   data=num_ASV,
                   family = "genpois")

simulationOutput <- simulateResiduals(fittedModel = ASV_mod) #ok
fit <- plot(simulationOutput, asFactor=TRUE) #ok
zi <- testZeroInflation(simulationOutput)
od <- testDispersion(simulationOutput)

plot(allEffects(ASV_mod)) 
em <- emmeans(ASV_mod, "run")
pairs(em) #A-B diff, A-D diff, B-C diff, B-D diff, C-D diff

###########################
#Concatenate by species and subset just prey ASVs ####
###########################

taxa_run <- all_run %>%
  left_join(taxa, by = "ASV") %>%
  filter(!unique_ID %in% c("Sparassidae", "Homo sapiens", "Hominidae", NA)) %>%
  group_by(run, sample, unique_ID) %>%
  summarise(reads = sum(reads))

prey_count <- taxa_run %>%
  group_by(run, sample) %>%
  filter(reads > 0) %>%
  tally(name = "species")

###########################
#Cross-run comparison of species composition by sample####
###########################

taxa_run <- taxa_run %>%
  mutate(presence = ifelse(reads >0, 1,0))

#going to do presence-absence composition analysis

sp1 <- glmmTMB(presence ~ run + (1|sample) + (1|unique_ID),
                data = taxa_run,
                family = 'binomial')

sp_null <- glmmTMB(presence ~ 1 + (1|sample) + (1|unique_ID),
                    data = taxa_run,
                    family = 'binomial')

AICc(sp1, sp_null)

binned_residuals(sp1) #looks not great, but not terrible
simulationOutput <- simulateResiduals(fittedModel = sp1) #ok
fit <- plot(simulationOutput, asFactor=TRUE) #ok

plot(allEffects(sp1)) #b is higher, C and D lower
em <- emmeans(sp1, "run")
pairs(em) #B-D different

###########################
#total species per sample by run####
###########################

spec_mod <- glmmTMB(species ~ run + (1|sample),
                   data=prey_count,
                   family = "genpois",
                   REML = FALSE)

spec_null <- glmmTMB(species ~ 1 + (1|sample),
                    data=prey_count,
                    family = "genpois",
                    REML = FALSE)

AICc(spec_mod, spec_null)

spec_mod <- glmmTMB(species ~ run + (1|sample),
                   data=prey_count,
                   family = "genpois")

simulationOutput <- simulateResiduals(fittedModel = spec_mod) #ok
fit <- plot(simulationOutput, asFactor=TRUE) #ok
zi <- testZeroInflation(simulationOutput)
od <- testDispersion(simulationOutput)

plot(allEffects(spec_mod)) 
em <- emmeans(spec_mod, "run")
pairs(em) #A-D diff, B-C diff, B-D diff, C-D diff


###########################
#Visualization and summary####
###########################
ggplot(num_ASV, aes(x = run, y = ASVs)) +
  geom_boxplot() + theme_bw() +
  labs(x = "Sequencing run", y = "ASVs per sample")

ggplot(prey_count, aes(x = run, y = species)) +
  geom_boxplot() + theme_bw() +
  labs(x = "Sequencing run", y = "Species per sample")
#summary, showing that on average, run B has 5 more ASVs per sample than
#run D, and 3 more ASVs per sample than run C
#not really sure how to correct this in the data... hmm...will think on it
num_ASV %>%
  group_by(run) %>%
  summarise(mean = mean(ASVs), total = n(), sd = sd(ASVs), se = sd/total)

###########################
#Supplement: By Sample visualization####
###########################
#by sample ASV presence
ggplot(all_run, aes(x = run, y = presence, fill = ASV)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  facet_wrap(~sample)

ggplot(taxa_run, aes(x = run, y = presence, fill = unique_ID)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  facet_wrap(~sample)
