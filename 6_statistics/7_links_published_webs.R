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
library(ggeffects)
###########################

###########################
# Load data
###########################
pal <- read.csv(here("data", "outputs",
                     "5_rarefied_taxonomic_sort",
                     "fam_prey_DNA.csv"))

pub_webs <- read.csv(here("data", "outputs", "6_pub_webs", "pub_webs_per_pred.csv"))

pub_webs <- pub_webs %>%
  dplyr::select(-X) %>%
  mutate(web_type = "published")

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
         family_richness = 69,
         coll_method = "HTS molecular",
         pub_year = 2020,
         web_type = "HTS") %>%
  rename("consumer" = "pred_ID")

###########################
# Combine and clean data
###########################

per_pred <- bind_rows(pal, pub_webs) 

per_pred$web <- as.factor(per_pred$web)
per_pred$coll_method <- as.factor(per_pred$coll_method)

sizes <- per_pred %>%
  distinct(species_richness) 

quantile(sizes$species_richness)

per_pred$web_sz <- ifelse(per_pred$species_richness <= 8, "zero",
                          ifelse(per_pred$species_richness > 8 & per_pred$species_richness <= 21, "twenty-five",
                                 ifelse(per_pred$species_richness > 21 & per_pred$species_richness <= 28, "fifty",
                                        ifelse(per_pred$species_richness > 28 & per_pred$species_richness <= 55, "seventy-five", "one hundred"))))

per_pred$web_sz <- as.factor(per_pred$web_sz)

sizes_fam <- per_pred %>%
  distinct(family_richness) 

quantile(sizes_fam$family_richness)

per_pred$web_sz_fam <- ifelse(per_pred$family_richness <=1, "zero", 
                              ifelse(per_pred$family_richness > 0 & per_pred$species_richness <= 4.25, "twenty-five",
                                     ifelse(per_pred$family_richness > 4.25 & per_pred$family_richness <= 8, "fifty",
                                            ifelse(per_pred$family_richness > 8 & per_pred$family_richness <= 13.75, "seventy-five", "one hundred"))))

per_pred$web_sz_fam <- as.factor(per_pred$web_sz_fam)
###########################
# visualizations
###########################

ggplot(per_pred, aes(x = web_type, y = links, fill = coll_method)) +
  geom_boxplot() +
  geom_vline(xintercept = 1.5, color = "grey") +
  theme_bw()

ggplot(per_pred, aes(x = web_type, y = links)) +
  geom_boxplot() +
  theme_bw()

ggplot(per_pred, aes(x = pub_year, y = species_richness)) +
  geom_point() +
  geom_hline(yintercept = 100, linetype = "dashed") +
  theme_bw()

ggplot(per_pred, aes(x = pub_year, y = family_richness)) +
  geom_point() +
  geom_hline(yintercept = 20, linetype = "dashed") +
  theme_bw()

hist(per_pred$links)

###########################
# stats: links for HTS vs published methods with species_richness
###########################
#what is the relationship between web collection method
#on the number of links per predator species if we take
#into account that these come from webs of different sizes?

m1 <- glmmTMB(links ~ coll_method + web_sz + (1|web),
              data = per_pred,
              family = "genpois")

summary(m1)
plot(allEffects(m1))

em <- emmeans(m1, "coll_method")
pairs(em)

em <- emmeans(m1, "web_sz")
pairs(em)

fit <- simulateResiduals(m1, plot = T)
testDispersion(fit)

me <- ggpredict(m1, terms = c("coll_method"))
plot(me) +
  labs(x = "Link assignment method", y = "Number of links per predator species") +
  theme(axis.title = element_text(size = 15), axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))

###########################
# stats: links for HTS vs published methods with family_richness
###########################
#what is the relationship between web collection method
#on the number of links per predator species if we take
#into account that these come from webs of different sizes?

m2 <- glmmTMB(links ~ coll_method + web_sz_fam + (1|web),
              data = per_pred,
              family = "genpois")

summary(m2)
plot(allEffects(m2))

em <- emmeans(m2, "coll_method")
pairs(em)

em <- emmeans(m2, "web_sz_fam")
pairs(em)

fit <- simulateResiduals(m2, plot = T)
testDispersion(fit)

me <- ggpredict(m2, terms = c("coll_method"))
plot(me) +
  labs(x = "Link assignment method", y = "Predicted links per predator species") +
  theme(axis.title = element_text(size = 15), axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))


###########################
# stats: links for HTS vs published methods with log transformed family richness
###########################
#what is the relationship between web collection method
#on the number of links per predator species if we take
#into account that these come from webs of different sizes?
per_pred <- per_pred %>%
  mutate(log_fam_rich = log(family_richness))

m3 <- glmmTMB(links ~ coll_method + log_fam_rich + (1|web),
              data = per_pred,
              family = "genpois")

summary(m3)
plot(allEffects(m3))

em <- emmeans(m3, "coll_method")
pairs(em)

fit <- simulateResiduals(m3, plot = T)
testDispersion(fit)

me <- ggpredict(m3, terms = c("coll_method"))
plot(me) +
  labs(x = "Link assignment method", y = "Predicted links per predator species") +
  theme(axis.title = element_text(size = 15), axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))

