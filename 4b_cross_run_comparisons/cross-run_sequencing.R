###########################
#Cross-run comparisons - sequencing depth
#Ana Miller-ter Kuile
#June 4, 2020
###########################

#This script looks at sequencing depth across runs to determine if some
#runs had higher sequencing depth than others. 

###########################
#Load Packages ####
###########################
library(tidyverse)
library(here)
library(ggplot2)
library(iNEXT)
library(cowplot)
library(glmmTMB)
library(MuMIn)
library(DHARMa)
library(emmeans)
library(effects)
###########################
#Load Data ####
###########################

samples <- read.csv(here("data", "outputs", "3_depth_corrected", "cross_run_samples.csv"))

###########################
#Clean data for sequencing depth ####
###########################

samples <- samples %>%
  dplyr::select(-X) %>%
  column_to_rownames(var = "ASV")


###########################
#Sequencing depth analysis ####
###########################

cross_depth <- iNEXT(samples, q=0, datatype="abundance") #this determines sequencing depth for each sample

#graph the interpolated and extrapolated sampling depth per sample
ggiNEXT(cross_depth, type=1, facet.var="none", grey = T, se = F) + 
  # can set se = F to remove shaded regions to see lines better 
  theme_bw() +
  labs(x = "Sequencing Depth", y = "ASV Richness", title = "DADA2 Sequencing Depth") +
  theme(legend.position = "none", axis.text = element_text(size = 20), 
        axis.title = element_text(size = 25))

###########################
#Sequencing depth analysis by run ####
###########################

run_a <- samples %>%
  dplyr::select(HEV07a:HEV29a)

a_depth <- iNEXT(run_a, q=0, datatype="abundance") #this determines sequencing depth for each sample

#graph the interpolated and extrapolated sampling depth per sample
a <- ggiNEXT(a_depth, type=1, facet.var="none", grey = T, se = F) + 
  # can set se = F to remove shaded regions to see lines better 
  theme_bw() +
  labs(x = "Sequencing Depth", y = "ASV Richness", title = "DADA2 Sequencing Depth") +
  theme(legend.position = "none", axis.text = element_text(size = 20), 
        axis.title = element_text(size = 25))

run_b <- samples %>%
  dplyr::select(HEV07b:HEV29b)

b_depth <- iNEXT(run_b, q=0, datatype="abundance") #this determines sequencing depth for each sample

#graph the interpolated and extrapolated sampling depth per sample
b <- ggiNEXT(b_depth, type=1, facet.var="none", grey = T, se = F) + 
  # can set se = F to remove shaded regions to see lines better 
  theme_bw() +
  labs(x = "Sequencing Depth", y = "ASV Richness", title = "DADA2 Sequencing Depth") +
  theme(legend.position = "none", axis.text = element_text(size = 20), 
        axis.title = element_text(size = 25))

run_c <- samples %>%
  dplyr::select(HEV07a:HEV29c)

c_depth <- iNEXT(run_c, q=0, datatype="abundance") #this determines sequencing depth for each sample

#graph the interpolated and extrapolated sampling depth per sample
c <- ggiNEXT(c_depth, type=1, facet.var="none", grey = T, se = F) + 
  # can set se = F to remove shaded regions to see lines better 
  theme_bw() +
  labs(x = "Sequencing Depth", y = "ASV Richness", title = "DADA2 Sequencing Depth") +
  theme(legend.position = "none", axis.text = element_text(size = 20), 
        axis.title = element_text(size = 25))

run_d <- samples %>%
  dplyr::select(HEV07a:HEV29d)

d_depth <- iNEXT(run_d, q=0, datatype="abundance") #this determines sequencing depth for each sample

#graph the interpolated and extrapolated sampling depth per sample
d <- ggiNEXT(d_depth, type=1, facet.var="none", grey = T, se = F) + 
  # can set se = F to remove shaded regions to see lines better 
  theme_bw() +
  labs(x = "Sequencing Depth", y = "ASV Richness", title = "DADA2 Sequencing Depth") +
  theme(legend.position = "none", axis.text = element_text(size = 20), 
        axis.title = element_text(size = 25))

plot_grid(a,b,c,d,ncol = 4, align = "vh")

###########################
#Sequencing depth values####
###########################

a_stats <- a_depth$DataInfo

a_stats <- a_stats %>%
  mutate(run = "A") #give a run specific ID

a_stats$site <- str_sub(a_stats$site, start = 1, end = -2) #remove "a" at end

b_stats <- b_depth$DataInfo

b_stats <- b_stats %>%
  mutate(run = "B")

b_stats$site <- str_sub(b_stats$site, start = 1, end = -2) #remove "b" at end

c_stats <- c_depth$DataInfo

c_stats <- c_stats %>%
  mutate(run = "C")

c_stats$site <- str_sub(c_stats$site, start = 1, end = -2) #remove "c" at end

d_stats <- d_depth$DataInfo

d_stats <- d_stats %>%
  mutate(run = "D")

d_stats$site <- str_sub(d_stats$site, start = 1, end = -2) #remove "d" at end

stats <- a_stats %>%
  bind_rows(b_stats) %>%
  bind_rows(c_stats) %>%
  bind_rows(d_stats)

###########################
#Sequencing depth mixed model####
###########################

stats$run <- factor(stats$run, levels = c("A", "B", "C", "D"))

stats$site <- as.factor(stats$site)

depth_mod <- glmmTMB(n ~ run + (1|site),
                     data=stats,
                     family = "genpois",
                     REML = FALSE)

depth_null <- glmmTMB(n ~ 1 + (1|site),
                     data=stats,
                     family = "genpois",
                     REML = FALSE)

AICc(depth_mod, depth_null)

depth_mod <- glmmTMB(n ~ run + (1|site),
                     data=stats,
                     family = "genpois")

simulation <- simulateResiduals(fittedModel = depth_mod)
fit <- plot(simulation, asFactor=TRUE)
zi <- testZeroInflation(simulation)
od <- testDispersion(simulation)

plot(allEffects(depth_mod))

em <- emmeans(depth_mod, "run")
pairs(em)
#only sig. difference is between A-D

ggplot(stats, aes(x = run, y = n)) +
  geom_boxplot() + theme_bw()
  