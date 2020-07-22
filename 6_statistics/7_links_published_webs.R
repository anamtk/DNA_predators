###########################
# Links per predator
# July 14, 2020
# Ana Miller-ter Kuile
###########################

#This code examines how many links per species for my 
#molecular web compared to those from published webs,
#all concatenated at the family level 

###########################
# Load packages
library(here)
library(tidyverse)
library(ggplot2)
library(glmmTMB)
library(DHARMa)
library(MuMIn)
library(effects)
library(emmeans)
###########################

###########################
# Load data
###########################
pal <- read.csv(here("data", "outputs",
                     "5_rarefied_taxonomic_sort",
                     "fam_prey_DNA.csv"))

pub_webs <- read.csv(here("data", "outputs", "6_pub_webs", "pub_webs_per_pred.csv"))

pub_webs <- pub_webs %>%
  dplyr::select(-X)

###########################
# Make Palmyra  link data
###########################

pal <- pal %>%
  group_by(pred_ID, Family) %>%
  summarise(reads = sum(reads)) %>%
  filter(reads > 0) %>%
  group_by(pred_ID) %>%
  tally(name = "links") %>%
  mutate(web = "Palmyra",
         species_richness = 409,
         coll_method = "HTS molecular",
         pub_year = 2020) %>%
  rename("consumer" = "pred_ID")

###########################
# Combine and visualize data
###########################

per_pred <- bind_rows(pal, pub_webs) 

per_pred$web <- as.factor(per_pred$web)
per_pred$coll_method <- as.factor(per_pred$coll_method)
per_pred$fct_yr <- as.factor(per_pred$pub_year)

###########################
# links per species by collection type
###########################

ggplot(per_pred, aes(x = coll_method, y = links)) +
  geom_boxplot() +
  theme_bw()

ggplot(per_pred, aes(x = coll_method, y = links/species_richness)) +
  geom_boxplot() +
  theme_bw()

ggplot(per_pred, aes(x = pub_year, y = species_richness)) +
  geom_point() +
  geom_smooth(method = 'lm', se = F) +
  theme_bw()

ggplot(per_pred, aes(x = pub_year, y = links)) +
  geom_point() +
  geom_smooth(method = 'lm', se = F) +
  theme_bw()

hist(per_pred$links)

pred_sp <- per_pred %>%
  group_by(web) %>%
  summarise(predators = n())

per_pred <- per_pred %>%
  left_join(pred_sp, by = "web")

#these stats are saying - correcting for web size
#what is the relationship between the number of links
#and the way in which links were assigned

m1 <- glmmTMB(links ~ coll_method + (1|species_richness),
              data = per_pred,
              family = "genpois")

summary(m1)
plot(allEffects(m1))

em <- emmeans(m1, "coll_method")
pairs(em)

fit <- simulateResiduals(m1, plot = T)
testDispersion(fit)


###########################
# stats with offset
###########################

m2 <- glmmTMB(links ~ coll_method + (1|species_richness),
              data = per_pred,
              offset = predators,
              family = "truncated_nbinom2")


summary(m2)
plot(allEffects(m2))

em <- emmeans(m2, "coll_method")
pairs(em)

fit <- simulateResiduals(m2, plot = T)
testDispersion(fit)
testUniformity(fit)
