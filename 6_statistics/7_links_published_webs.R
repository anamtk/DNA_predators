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
# Load and tidy family data
###########################
pal <- read.csv(here("data", "outputs",
                     "5_rarefied_taxonomic_sort",
                     "fam_prey_DNA.csv"))

pal_fam <- pal %>%
  distinct(Family) %>%
  tally(name = "fam_richness")

hines_fam <- read.csv(here("data", "outputs",
                           "6_pub_webs", 
                           "hines_predation.csv"))

hines_fam <- hines_fam %>%
  distinct(taxon_Family) %>%
  tally(name = "fam_richness")

laigle_fam <- read.csv(here("data", "outputs",
                            "6_pub_webs", 
                            "laigle_predation.csv"))

laigle_fam <- laigle_fam %>%
  distinct(Family) %>%
  tally(name = "fam_richness")

###########################
# Load and tidy predation link data
###########################

pal <- pal %>%
  group_by(pred_ID, Family) %>%
  summarise(reads = sum(reads)) %>%
  filter(reads > 0) %>%
  group_by(pred_ID) %>%
  tally(name = "links") %>%
  mutate(web = "Palmyra",
         fam_richness = 69)

hines <- read.csv(here("data", "outputs",
                       "6_pub_webs", 
                       "hines_per_pred.csv"))

hines <- hines %>%
  dplyr::select(-X) %>%
  rename("pred_ID" = "node_from") %>%
  mutate(web = "Hines",
         fam_richness = 84)
hines$pred_ID <- as.character(hines$pred_ID)

laigle <- read.csv(here("data", "outputs",
                        "6_pub_webs", 
                        "laigle_per_pred.csv"))

laigle <- laigle %>%
  dplyr::select(-X) %>%
  rename("pred_ID" = "consumer") %>%
  mutate(web = "Laigle",
         fam_richness = 22)
laigle$pred_ID <- as.character(laigle$pred_ID)

###########################
# Combine and visualize data
###########################
links <- pal %>%
  bind_rows(hines) %>%
  bind_rows(laigle)
links$web <- as.factor(links$web)
ggplot(links, aes(x = web, y = links/fam_richness)) +
  geom_boxplot() + theme_bw()

ggplot(links, aes(x = web, y = links)) +
  geom_boxplot() + theme_bw() +
  geom_segment(aes(x = 0.5, xend = 1.5,
                   y = 84, yend = 84)) +
  geom_segment(aes(x = 1.5, xend = 2.5,
                   y = 22, yend = 22)) +
  geom_segment(aes(x = 2.5, xend = 3.5,
                   y = 69, yend = 69))

ggplot(links, aes(x = web, y = fam_richness)) +
  geom_point() + theme_bw()

###########################
# stats with offset
###########################

m1 <- glmmTMB(links ~ web,
              data = links,
              offset = fam_richness,
              family = "poisson",
              REML = FALSE)

m2 <- glmmTMB(links ~ 1,
              data = links,
              family = "poisson",
              offset = fam_richness,
              REML = FALSE)

AICc(m1, m2)

m1 <- glmmTMB(links ~ web,
              data = links,
              offset = fam_richness,
              family = "poisson")

summary(m1)
plot(allEffects(m1))

em <- emmeans(m1, "web")
pairs(em)

###########################
# stats without offset
###########################

m3 <- glmmTMB(links ~ web,
              data = links,
              family = "poisson",
              REML = FALSE)

m4 <- glmmTMB(links ~ 1,
              data = links,
              family = "poisson",
              REML = FALSE)

AICc(m3, m4)

m3 <- glmmTMB(links ~ web,
              data = links,
              family = "poisson")

summary(m3)
plot(allEffects(m3))

em <- emmeans(m3, "web")
pairs(em)

