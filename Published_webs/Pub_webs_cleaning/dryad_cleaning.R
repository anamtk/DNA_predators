###########################
# Add Dryad food webs to rmangal food webs
# July 13, 2020
# Ana Miller-ter Kuile
###########################

#need to add Dryad food webs to the rmangal food
#web dataframe, making sure they have consistent
#haming schemes.

#Given that I want to correct for total species richness
#going to only consider food webs in this analysis
#which then only includes the Laigle and Rohr datasets

###########################
# Load packages
library(here)
library(tidyverse)
library(ggplot2)
###########################

###########################
# Load data
###########################

laigle_ints <- read.csv(here("Published_webs", "Dryad", 
                             "Laigle_etal_2017", "interactions_Laigle.csv"))

laigle_nodes <- read.csv(here("Published_webs", "Dryad", 
                             "Laigle_etal_2017", "nodes_Laigle.csv"))

rohr_nodes <- read.csv(here("Published_webs", "Dryad", 
                          "Rohr_etal2014", "fw_broom.csv"))

rohr_ints <- read.csv(here("Published_webs", "Dryad", 
                           "Rohr_etal2014", "body_sizes.csv"))

rmangal_webs <- read.csv(here("data", "outputs", "6_pub_webs", "rmangal_webs_predation.csv"))

rmangal_webs <- rmangal_webs %>%
  mutate_at('consumer', as.character) %>%
  mutate_at('resource', as.character)
###########################
# Manipulation how-to
###########################

#need to manipulate both of these to include these variables
#resource, consumer, type, method, original_name, taxon_Kingdom
#taxon_Class, taxon_Order, taxon_Family, taxonomy.name,
#species_richness, web

###########################
# Laigle
###########################

laigle_nodes <- laigle_nodes %>%
  add_tally(name = "species_richness")

laigle_web <- laigle_ints %>%
  filter(interaction > 0) %>%
  dplyr::select(consumer, resource) %>%
  left_join(laigle_nodes, by = c("resource" = "X")) %>%
  filter(Kingdom == "Animalia") %>%
  mutate(web = "laigle") %>%
  dplyr::select(consumer, resource, Family, Order, 
                Class, Phylum, Kingdom, sp, species_richness, web) %>%
  rename("original_name" = "sp",
         "taxon_Kingdom" = "Kingdom",
         "taxon_Class" = "Class",
         "taxon_Order" = "Order",
         "taxon_Family" = "Family") %>%
  mutate(type = "predation",
         method = "various (observed, feeding trial, non-HTS molecular, inferred from co-occurence and body size rules)",
         taxonomy.name = original_name) %>%
  mutate_at("consumer", as.character) %>%
  mutate_at("resource", as.character)

###########################
# Rohr
###########################
rohr_nodes %>%
  tally()

rohr_web <- rohr_ints %>%
  filter(Type.of.feeding.interaction == "predacious") %>%
  dplyr::select(Link.methodology, Taxonomy.consumer, Type.of.feeding.interaction,
                Taxonomy.resource, taxon_Family, taxon_Order) %>%
  rename("method" = "Link.methodology",
         "consumer" = "Taxonomy.consumer",
         "type" = "Type.of.feeding.interaction",
         "resource" = "Taxonomy.resource") %>%
  mutate(species_richness = 55,
         web = "rohr",
         taxon_Class = NA,
         taxon_Kingdom= "Animalia",
         original_name = resource,
         taxonomy.name = resource)

###########################
# Bind all webs together
###########################

predation <- bind_rows(rmangal_webs, rohr_web, laigle_web) %>%
  dplyr::select(-X, -Phylum)


###########################
# Per predator links
###########################

per_pred <- predation %>%
  group_by(consumer, taxon_Family, species_richness, web) %>%
  tally(name = "links") %>%
  group_by(consumer, species_richness, web) %>%
  summarise(links = n())

ggplot(per_pred, aes(x = species_richness, y = links)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw()

ggplot(per_pred, aes(x = web, y = links)) +
  geom_boxplot() +
  theme_bw()

ggplot(per_pred, aes(x = web, y = links/species_richness)) +
  geom_boxplot() +
  theme_bw()

predation %>%
  distinct(web)

method <- predation %>%
  group_by(web, method) %>%
  summarise(sum = n())

web <- unique(predation$web)

coll_method <- c("observation and literature", "observation and literature", 
            "observation", "literature",
            "observation", "observation", "observation", "observation",
            "observation", "observation", "observation", "observation", 
            "observation", "observation", "observation", "observation",
            "observation", "observation", "observation", "observation",
            "observation and literature")

pub_year <- c(1983, 1983, 1974, 2019, 1978, 1978, 1978, 1978, 1943, 
              1926, 1953, 1939, 1939, 1939, 1939, 1983, 1974, 1974, 
              1984, 2000, 2017)

#author <- c("Choenly", "Choenly", "Cornaby", "Hines", "McKinnerney",
#            "McKinnerney", "McKinnerney", "McKinnerney", "Mohr",
#            "Richards", "Robinson", "Savely", "Savely", "Savely",
#            "Savely", "Schoenly", "Valeila", "Valeila", "Whittaker",
#             "Rohr", "Laigle")

methods <- as.data.frame(cbind(web, coll_method, pub_year))

predation <- predation %>%
  dplyr::select(-method) %>%
  left_join(methods, by = "web")

per_pred <- per_pred %>%
  left_join(methods, by = "web")

###########################
# Export data
###########################
write.csv(predation, here("data", "outputs", 
                       "6_pub_webs", "pub_webs_predation.csv"))

write.csv(per_pred, here("data", "outputs", 
                         "6_pub_webs", "pub_webs_per_pred.csv"))



