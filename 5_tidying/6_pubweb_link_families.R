###########################
# Links families
# July 22, 2020
# Ana Miller-ter Kuile
###########################

#going to be looking at how many families are in our
#links across webs. 

###########################
# Load packages
library(here)
library(tidyverse)
library(ggplot2)
library(glmmTMB)
library(emmeans)
library(DHARMa)
library(effects)
###########################

###########################
# Load data
###########################

pal <- read.csv(here("data", "outputs",
                     "5_rarefied_taxonomic_sort",
                     "fam_prey_DNA.csv"))


pub <- read.csv(here("data", "outputs", "6_pub_webs", "pub_webs_predation.csv"))

###########################
# Unique families
###########################

pal_fams <- pal %>%
  distinct(Family) %>%
  rename("taxon_Family" = "Family")

all_fams <- pub %>%
  distinct(taxon_Family) %>%
  bind_rows(pal_fams)

###########################
# OUtput to Google trophic groups
###########################

write.csv(all_fams, here("data", "outputs", "7_all_webs", 
                         "diet_families.csv"))

###########################
# Import the updated version
###########################

all_fams <- read.csv(here("data", "outputs", "7_all_webs", 
                          "diet_families_tp.csv"))

###########################
# Combine Pal and published web data
###########################
pub <- pub %>%
  dplyr::select(-X, -type, -taxonomy.name, -taxon_Kingdom,
                -taxon_Class, -taxon_Order, )

pal <- pal %>%
  group_by(pred_ID, Family) %>%
  summarise(reads = sum(reads)) %>%
  filter(reads > 0) %>%
  mutate(web = "Palmyra",
         species_richness = 409,
         coll_method = "HTS molecular",
         pub_year = 2020,
         original_name = Family,
         resource = Family) %>%
  rename("consumer" = "pred_ID",
         "taxon_Family" = "Family") %>%
  dplyr::select(-reads)

all_fam <- bind_rows(pal, pub) %>%
  left_join(all_fams, by = "taxon_Family") %>%
  dplyr::select(-X) %>%
  unite(tp, trophic_group, predatory_omnivory, sep = "", remove = FALSE)

###########################
# Create trophic position variables
###########################

all_fam$tp <- factor(all_fam$tp, levels = c("herbivore", "detritovore", "omnivoreno",
                                            "omnivoreyes", "parasite", "parasitoid", "predator"))

all_fam$tl <- ifelse(all_fam$tp == "herbivore", 1,
                     ifelse(all_fam$tp == "detritovore", 1,
                            ifelse(all_fam$tp == "omnivoreno", 1,
                                   ifelse(all_fam$tp == "omnivoreyes", 2,
                                          ifelse(all_fam$tp == "parasite", 3, 4)))))

all_fam$broad_tl <- ifelse(all_fam$broad_tp == "basal", 1,
                           ifelse(all_fam$broad_tp == "omnivorous", 2, 
                                  ifelse(all_fam$broad_tp =="", NA, 3)))
###########################
# Quick Vis
###########################

all_fam %>%
  group_by(coll_method) %>%
  summarise(mean = mean(tl, na.rm=T), 
            total = n(), 
            sd = sd(tl, na.rm=T), 
            se = sd/sqrt(total)) #%>%
  ggplot(aes(x = coll_method, y = mean)) +
  geom_point() +
  geom_errorbar(ymin = mean - se, ymax = mean + se) +
  theme_bw()

ggplot(all_fam, aes(x = coll_method, y = broad_tl)) +
  geom_boxplot() +
  theme_bw()

ggplot(all_fam, aes(x = web, y = broad_tl)) +
  geom_boxplot() +
  theme_bw()

all_fam %>%
  group_by(consumer, web, coll_method, species_richness, broad_tp) %>%
  summarise(proportion = n()/species_richness) %>%
  filter(broad_tp != "") %>%
  ggplot(aes(x = coll_method, y = proportion, fill = broad_tp)) +
  geom_boxplot() +
  theme_bw() +
  scale_y_log10()

###########################
# Quick stats
###########################

tp_df <- all_fam %>%
  group_by(consumer, web, coll_method, species_richness, broad_tp) %>%
  summarise(count = n()) %>%
  mutate(log_sp_rich = log(species_richness))

m1 <- glmmTMB(count ~ broad_tp + offset(log_sp_rich) + (1|coll_method),
              data = tp_df,
              family = "genpois")

summary(m1)
plot(allEffects(m1))

library(ggeffects)

me <- ggpredict(m1, terms = c("broad_tp", "coll_method"), type = "re")
plot(me)
