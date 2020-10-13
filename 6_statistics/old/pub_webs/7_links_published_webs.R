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
package.list <- c("here", "tidyverse", "ggplot2", "glmmTMB", "DHARMa", 
                  "MuMIn", "effects", "emmeans", "ggeffects")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}
#############################


###########################
# Load data
###########################
pal <- read.csv(here("data", "outputs",
                     "5_rarefied_taxonomic_sort",
                     "fam_prey_DNA.csv"))

pub_webs <- read.csv(here("data", "outputs", 
                          "6_pub_webs", 
                          "pub_webs_per_pred.csv"))

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

ggplot(per_pred, aes(x = web_type, y = links/family_richness, fill = coll_method)) +
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

hist(per_pred$links)

###########################
# stats: links for HTS vs published methods with family_richness
###########################
#what is the relationship between web collection method
#on the number of links per predator species if we take
#into account that these come from webs of different sizes?
hist(log(per_pred$links))

m2 <- glmmTMB(links ~ coll_method  + (1|web),
              data = per_pred,
              offset = log(family_richness),
              family = "nbinom2")

summary(m2)
plot(allEffects(m2))

em <- emmeans(m2, "coll_method")
pairs(em)

fit <- simulateResiduals(m2, plot = T) #KS test significant
testUniformity(fit) 
testDispersion(fit)

m2 <- glmmTMB(links ~ coll_method  + offset(log(family_richness)) + (1|web),
              data = per_pred,
              family = "nbinom2")

me <- ggpredict(m2, terms = c("coll_method"))
plot(me) +
  labs(x = "Link assignment method", y = "Predicted links per predator species") +
  theme(axis.title = element_text(size = 15), axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))

###########################
# stats: links for HTS vs published methods without family_richness
###########################

m2 <- glmmTMB(links ~ coll_method  + (1|web),
              data = per_pred,
              family = "truncated_nbinom2")


summary(m2)
plot(allEffects(m2))

em <- emmeans(m2, "coll_method")
pairs(em)


fit <- simulateResiduals(m2, plot = T) #KS test significant
testUniformity(fit) 
testDispersion(fit)

me <- ggpredict(m2, terms = c("coll_method"))
plot(me) +
  labs(x = "Link assignment method", y = "Predicted links per predator species") +
  theme(axis.title = element_text(size = 15), axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))

###########################
# stats: links for HTS vs published methods new webs only
###########################
new_pre_pred <- per_pred %>%
  filter(pub_year > 2000)

m2 <- glmmTMB(links ~ web_type + (1|web),
              data = new_per_pred,
              family = "truncated_genpois")


summary(m2)
plot(allEffects(m2))

fit <- simulateResiduals(m2, plot = T) #KS test significant
testUniformity(fit) 
testDispersion(fit)

me <- ggpredict(m2, terms = c("web_type"))
plot(me) +
  labs(x = "Link assignment method", y = "Predicted links per predator species") +
  theme(axis.title = element_text(size = 15), axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(new_per_pred, aes(x = web_type, y = links)) +
  geom_boxplot() +
  theme_bw()
