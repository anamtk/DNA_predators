#############################
#Trophic positions of published vs. HTS links
#August 18, 2020
#Ana Miller-ter Kuile
#############################

#This script examines whether the trophic positions of diet 
#items assigned via different published methods differ
#from those from the HTS data
#This is important because of the importance of intraguild
#predation in food web ecology

#############################
## Bring in Necessary Packages
package.list <- c("here", "tidyverse", "ggplot2", "glmmTMB", "DHARMa", 
                  "MuMIn", "effects", "emmeans", "ggeffects")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

#############################

#note, it may be important to break this analysis up
#by first: more specific trophic groups and then
#by courser scale - e.g. are they eating plants, vs. 
#other animals, or both.

#############################
#Load data
#############################

links <- read.csv(here("data", "outputs", 
                       "7_all_webs", 
                       "all_interactions_and_tp.csv"))

#############################
#Summarize by specific TP
#############################

links_tp <- links %>%
  unite(web_consumer, consumer, web, remove = F) %>%
  group_by(web_consumer, web, family_richness, coll_method,
           pub_year, tp) %>%
  tally(name = "link_number") 

#############################
#Summarize by broad TP
#############################

links_btp <- links %>%
  unite(web_consumer, consumer, web, remove = F) %>%
  group_by(web_consumer, web, family_richness, coll_method,
           pub_year, broad_tp) %>%
  tally(name = "link_number") 

###########################
# total link by broad TP, and total number of each tidying
###########################
links_btp <- links_btp %>%
  filter(broad_tp != "")

links_sp <- links %>%
  unite(web_consumer, consumer, web, remove = F) %>%
  group_by(web_consumer, web, family_richness, coll_method,
           pub_year) %>%
  tally(name = "total_link_number") %>%
  ungroup() %>%
  dplyr::select(web_consumer, total_link_number)

web_species_tp <- links %>%
  group_by(web, broad_tp) %>%
  distinct(resource) %>%
  summarise(total = n()) %>%
  filter(broad_tp != "") %>%
  pivot_wider(names_from = "broad_tp",
              values_from = "total")

web_species_tp[is.na(web_species_tp)] <- 0

links_btp <- links_btp %>%
  left_join(links_sp, by = "web_consumer") %>%
  left_join(web_species_tp, by = "web")

ggplot(links_btp, aes(x = coll_method, y = link_number/total_link_number, fill = broad_tp)) +
  geom_boxplot() + theme_bw()

###########################
# subset different TP for analyses
###########################
bt_basal <- links_btp %>%
  filter(broad_tp == "basal")

bt_omni <- links_btp %>%
  filter(broad_tp == "omnivorous")

bt_pred <- links_btp %>%
  filter(broad_tp == "predatory")

###########################
# basal TP analysis
###########################
bt_basal <- bt_basal %>%
  mutate(log_basal = log(basal),
         proportion = link_number/total_link_number)

m1 <- glmmTMB(link_number ~ coll_method + log_basal + (1|web),
              data = bt_basal,
              offset= log(total_link_number),
              family = "nbinom2")


summary(m1)
plot(allEffects(m1))

em <- emmeans(m1, "coll_method")
pairs(em)

fit <- simulateResiduals(m1, plot = T) 
testUniformity(fit)
testDispersion(fit)

m1 <- glmmTMB(link_number ~ coll_method + log_basal + offset(log(total_link_number))+ (1|web),
              data = bt_basal,
              family = "nbinom2")

me <- ggpredict(m1, terms = c("coll_method"))
plot(me) +
  labs(x = "Link assignment method", y = "Predicted basal links per predator species") +
  theme(axis.title = element_text(size = 15), axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))

#proportional model

m2 <- glmmTMB(proportion ~ coll_method + (1|web),
              data = bt_basal,
              weights = total_link_number,
              family = "binomial")

summary(m2)
plot(allEffects(m2))
em <- emmeans(m2, "coll_method")
pairs(em)
fit <- simulateResiduals(m2, plot = T) 

me <- ggpredict(m2, terms = c("coll_method"))
plot(me) +
  labs(x = "Link assignment method", y = "Predicted basal links per predator species") +
  theme(axis.title = element_text(size = 15), axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))

###########################
# omni TP analysis
###########################
bt_omni <- bt_omni %>%
  mutate(log_omni = log(omnivorous),
         proportion = link_number/total_link_number)

m1 <- glmmTMB(link_number ~ coll_method  + (1|web),
              data = bt_omni,
              offset = log(omnivorous),
              family = "genpois")

summary(m1)
plot(allEffects(m1))

em <- emmeans(m1, "coll_method")
pairs(em)

fit <- simulateResiduals(m1, plot = T) 
testUniformity(fit)
testDispersion(fit)

m1 <- glmmTMB(link_number ~ coll_method + offset(log(omnivorous)) + (1|web),
              data = bt_omni,
              family = "genpois")

me <- ggpredict(m1, terms = c("coll_method"))
plot(me) +
  labs(x = "Link assignment method", y = "Predicted omnivore links per predator species") +
  theme(axis.title = element_text(size = 15), axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))

#binomial propprtion model
m2 <- glmmTMB(proportion ~ coll_method + (1|web),
              data = bt_omni,
              weights = total_link_number,
              family = "binomial")

summary(m2)
plot(allEffects(m2))
em <- emmeans(m2, "coll_method")
pairs(em)
fit <- simulateResiduals(m2, plot = T) 

me <- ggpredict(m2, terms = c("coll_method"))
plot(me) +
  labs(x = "Link assignment method", y = "Predicted basal links per predator species") +
  theme(axis.title = element_text(size = 15), axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))


###########################
# pred TP analysis
###########################
bt_pred <- bt_pred %>%
  mutate(log_pred = log(predatory),
         proportion = link_number/total_link_number,
         web_type = ifelse(coll_method == "HTS molecular", "HTS", "literature")) 

m1 <- glmmTMB(link_number ~ web_type  + (1|web),
              data = bt_pred,
              offset = log(predatory),
              family = "nbinom2")

summary(m1)
plot(allEffects(m1))

em <- emmeans(m1, "web_type")
pairs(em)

fit <- simulateResiduals(m1, plot = T) 
testUniformity(fit)
testDispersion(fit)

m1 <- glmmTMB(link_number ~ coll_method + offset(log(predatory)) + (1|web),
              data = bt_pred,
              family = "genpois")

me <- ggpredict(m1, terms = c("coll_method"))
plot(me) +
  labs(x = "Link assignment method", y = "Predicted predatory links per predator species") +
  theme(axis.title = element_text(size = 15), axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))

#binomial propprtion model
m2 <- glmmTMB(proportion ~ web_type + (1|web) + (1|coll_method),
              data = bt_pred,
              weights = total_link_number,
              family = "betabinomial")

summary(m2)
plot(allEffects(m2))
em <- emmeans(m2, "coll_method")
pairs(em)
fit <- simulateResiduals(m2, plot = T) 

me <- ggpredict(m2, terms = c("web_type"))
plot(me) +
  labs(x = "Link assignment method", y = "Predicted basal links per predator species") +
  theme(axis.title = element_text(size = 15), axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(bt_pred, aes(x = web_type, y = proportion)) +
  geom_boxplot() +theme_bw()


###########################
# bayesian explorations
###########################
install.packages("MCMCglmm")
library(MCMCglmm)

m1 <- glmmTMB(link_number ~ coll_method + offset(log(predatory)) + (1|web),
              data = bt_pred,
              family = "genpois")

mod = MCMCglmm(log(Species+1) ~ log(Area) + log(Elev+1) + Temp +
                 log(number.islands) + log(distance) + log(age), random = ~Archipelago,
               data = islands.sub)

mb <- MCMCglmm(link_number ~ coll_method + log(total_link_number), 
               random = ~web,
               data = bt_pred)

plot(mb$Sol[,1])

acfplot(mb$VCV)

summary(mb)

?MCMCglmm


