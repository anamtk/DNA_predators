##########################
# 1. Is there stage structure? -----
# Ana Miller-ter Kuile
# October 20, 2020
###########################

# this script does k-means clustering to
#determine if there are clusters in three
#predator species for which we have the largest
#sample sizes. maybe extend to others too? not sure
#then looks at whether these stages are predicted by 
#body size. 


# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 
                  "glmmTMB", "emmeans",
                  "MuMIn", "DHARMa",
                  "effects", "ggeffects",
                  "vegan")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Import data -------------------------------------------------------------

data <- read.csv(here("data", "outputs", "8_final_dataset",
                      "pred_prey_sizes_tp_DNAinteractions.csv"))

size <- data %>%
  dplyr::select(-X, -X.x, -X.y, -reads) %>%
  mutate(pred_mass_mg = exp(pred_log_mass_mg))

# Heteropoda venatoria ----------------------------------------------------
#species: 

sppData <- split(size, size$sample_str)

process <- function(x, sppSplit = sppData){
  species_ints <- sppSplit[[x]] %>%
    #filter(sample_str == species) %>%
    dplyr::select(sample, Family) %>%
    group_by(sample, Family) %>%
    mutate(presence = 1) %>%
    mutate(Family = factor(Family)) %>%
    mutate(presence = replace_na(presence, 0)) %>%
    pivot_wider(names_from = Family,
                values_from = presence,
                values_fill = 0) %>%
    column_to_rownames(var = "sample")
  
  maxClust <- min(ncol(species_ints), nrow(species_ints))
  
  cluster_model <- cascadeKM(species_ints, 1, maxClust, iter = 10000)
  plot(cluster_model, sortg=T)
  return(cluster_model$results)
}
lapply(1:length(sppData), FUN = process)

species <- size %>%
  dplyr::select(sample_str) %>%
  distinct(sample_str) %>%
  pull(sample_str)

for(i in length(species)){
  plots <- cluster_stages()
}


?as.vector
ints <- size %>%
  filter(sample_str == species)

SCY <- size %>%
  filter(sample_str == "SCY") %>%
  dplyr::select(sample, Family) %>%
  group_by(sample, Family) %>%
  mutate(presence = 1) %>%
  mutate(Family = factor(Family)) %>%
  mutate(presence = replace_na(presence, 0)) %>%
  pivot_wider(names_from = Family,
              values_from = presence,
              values_fill = 0) %>%
  column_to_rownames(var = "sample")

SCY_dist <- vegdist(SCY, method = "jaccard", binary = TRUE)
model_SCY <- cascadeKM(SCY_dist, 1, 6, iter = 10000)
plot(model_SCY, sortg=T)
model_SCY
SCY_clusts <- as.data.frame(model_SCY$partition)

SCY_clusts <- SCY_clusts %>%
  dplyr::select("6 groups") %>%
  rownames_to_column(var = "sample") %>%
  rename("cluster" = "6 groups")

SCY_size <- size %>%
  filter(sample_str == "SCY") %>%
  left_join(SCY_clusts, by = "sample") %>%
  mutate(cluster = as.factor(cluster))

ggplot(SCY_size, aes(x = cluster, y = pred_mass_mg)) +
  geom_boxplot()

SCY_size %>%
  group_by(cluster, broad_tp) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = broad_tp,
              values_from = count) %>%
  mutate(total = sum(basal + omnivorous + predatory, na.rm=T),
         intraguild = sum(omnivorous + predatory)) %>%
  ggplot(aes(x = cluster, y = intraguild/total)) +
  geom_bar(stat = "identity")

# NEO ---------------------------------------------------------------------

NEO <- size %>%
  filter(sample_str == "NEO") %>%
  dplyr::select(sample, Family) %>%
  group_by(sample, Family) %>%
  mutate(presence = 1) %>%
  mutate(Family = factor(Family)) %>%
  mutate(presence = replace_na(presence, 0)) %>%
  pivot_wider(names_from = Family,
              values_from = presence,
              values_fill = 0) %>%
  column_to_rownames(var = "sample")

model_NEO <- cascadeKM(NEO, 1, 5, iter = 10000)
plot(model_NEO, sortg=T)
model_NEO$results
NEO_clust <- as.data.frame(model_NEO$partition)

NEO_clust <- NEO_clust %>%
  dplyr::select("2 groups") %>%
  rownames_to_column(var = "sample") %>%
  rename("cluster" = "2 groups")

NEO_size <- size %>%
  filter(sample_str == "NEO") %>%
  left_join(NEO_clust, by = "sample") %>%
  mutate(cluster = as.factor(cluster)) 

ggplot(NEO_size, aes(x = cluster, y = pred_mass_mg)) +
  geom_boxplot()

NEO_size %>%
  group_by(cluster, broad_tp) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = broad_tp,
              values_from = count) %>%
  mutate(total = sum(basal + omnivorous + predatory),
         intraguild = sum(omnivorous + predatory)) %>%
  ggplot(aes(x = cluster, y = intraguild/total)) +
  geom_bar(stat = "identity")

# PHH ---------------------------------------------------------------------

PHH <- size %>%
  filter(sample_str == "PHH") %>%
  dplyr::select(sample, Family) %>%
  group_by(sample, Family) %>%
  mutate(presence = 1) %>%
  mutate(Family = factor(Family)) %>%
  mutate(presence = replace_na(presence, 0)) %>%
  pivot_wider(names_from = Family,
              values_from = presence,
              values_fill = 0) %>%
  column_to_rownames(var = "sample")

model_PHH <- cascadeKM(PHH, 1, 5, iter = 10000)
plot(model_PHH, sortg=T)
model_PHH$results
PHH_clust <- as.data.frame(model_PHH$partition)

PHH_clust <- PHH_clust %>%
  dplyr::select("2 groups") %>%
  rownames_to_column(var = "sample") %>%
  rename("cluster" = "2 groups")

PHH_size <- size %>%
  filter(sample_str == "PHH") %>%
  left_join(PHH_clust, by = "sample") %>%
  mutate(cluster = as.factor(cluster))

ggplot(PHH_size, aes(x = cluster, y = pred_mass_mg)) +
  geom_boxplot()

PHH_size %>%
  group_by(cluster, broad_tp) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = broad_tp,
              values_from = count) %>%
  mutate(total = sum(basal + omnivorous + predatory)) %>%
  pivot_longer(cols = basal:predatory,
               names_to = "broad_tp",
               values_to = "count") %>%
  ggplot(aes(x = cluster, y = count/total, fill = broad_tp)) +
  geom_bar(stat = "identity", position = "dodge")
