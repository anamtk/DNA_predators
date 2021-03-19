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

package.list <- c("here", "tidyverse", "glmmTMB",
                  "DHARMa", "effects", "emmeans",
                  "MuMIn", "lattice", "performance")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

###########################
#Load Data####
###########################

#rarefied cross-run samples
cross <- read.csv(here("data", 
                       "outputs", 
                       "3b_depth_corrected", 
                       "cross_run_samples.csv"))

taxa <- read.csv(here("data", 
                      "outputs", 
                      "3a_remove_negatives", 
                      "ASVs_counts_all.csv"))

###########################
#Manipulate DF to long for analyses####
###########################

#Divide the samples up by sequencing run, give them a run category to distinguuish
#them, and then remove the a, b, c, d notation for GLMM below

#Run A
run_a <- cross %>%
  dplyr::select(ASV, HEV07a:HEV29a) %>% #select run 1
  gather(sample, reads, HEV07a:HEV29a) %>% #make this long
  mutate(run = "A") #give a run specific ID

run_a$sample <- str_sub(run_a$sample, start = 1, end =-2) #remove "a" at end

#Run B
run_b <- cross %>%
  dplyr::select(ASV, HEV07b:HEV29b) %>%
  gather(sample, reads, HEV07b:HEV29b) %>%
  mutate(run = "B")

run_b$sample <- str_sub(run_b$sample, start = 1, end =-2)

#Run C
run_c <- cross %>%
  dplyr::select(ASV, HEV07c:HEV29c) %>%
  gather(sample, reads, HEV07c:HEV29c) %>%
  mutate(run = "C")

run_c$sample <- str_sub(run_c$sample, start = 1, end =-2)

#Run D
run_d <- cross %>%
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

#Factor variables
all_run$run <- factor(all_run$run, levels = c("A", "B", "C", "D"))
all_run$ASV <- as.factor(all_run$ASV)
all_run$sample <- as.factor(all_run$sample)

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

simulationOutput <- simulateResiduals(fittedModel = ASV_mod, plot = T)
zi <- testZeroInflation(simulationOutput)
od <- testDispersion(simulationOutput)

plot(allEffects(ASV_mod)) 
em <- emmeans(ASV_mod, "run")
pairs(em) #A-B diff, A-D diff, B-C diff, B-D diff, C-D diff

###########################
#Concatenate by family and subset just prey ASVs ####
###########################

## AB note: lines 147-155 don't work for me
## error message:
# Error: Problem with `filter()` input `..1`.
# x object 'Family' not found
# ℹ Input `..1` is `!Family %in% c("Sparassidae", "Hominidae", "Muridae", NA, "Salmonidae")`.
# ℹ The error occurred in group 1: ASV = "ASV_1".
## I tried to look at `all_run` and `taxa` but neither of them have taxonomy
## I also couldn't find the taxonomy file!

taxa_run <- all_run %>%
  left_join(taxa, by = "ASV") %>%
  filter(!Family %in% c("Sparassidae", 
                        "Hominidae", 
                        "Muridae", 
                        NA,
                        "Salmonidae")) %>%
  group_by(run, sample, Family) %>%
  summarise(reads = sum(reads))

prey_count <- taxa_run %>%
  group_by(run, sample) %>%
  filter(reads > 0) %>%
  tally(name = "families")

###########################
#total families per sample by run####
###########################

fam_mod <- glmmTMB(families ~ run + (1|sample),
                   data=prey_count,
                   family = "genpois",
                   REML = FALSE)

fam_null <- glmmTMB(families ~ 1 + (1|sample),
                    data=prey_count,
                    family = "genpois",
                    REML = FALSE)

AICc(fam_mod, fam_null)

fam_mod <- glmmTMB(families ~ run + (1|sample),
                   data=prey_count,
                   family = "genpois")

simulationOutput <- simulateResiduals(fittedModel = fam_mod, plot =T) #ok
testZeroInflation(simulationOutput)
testDispersion(simulationOutput)

plot(allEffects(fam_mod)) 
em <- emmeans(fam_mod, "run")
em
pairs(em) #A-D diff, B-C diff, B-D diff

prey_count %>%
  group_by(run) %>%
  summarise(mean = mean(families),
            sd = sd(families),
            total = n(),
            se = sd/sqrt(total))
  
###########################
#Visualization and summary####
###########################
ggplot(num_ASV, aes(x = run, y = ASVs)) +
  geom_boxplot() + theme_bw() +
  labs(x = "Sequencing run", y = "Total ASVs per sample")

ggplot(prey_count, aes(x = run, y = families)) +
  geom_boxplot() + theme_bw() +
  labs(x = "Sequencing run", y = "Prey families per sample")

