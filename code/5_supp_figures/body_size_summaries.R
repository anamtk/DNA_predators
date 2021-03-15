# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", "bipartite", "patchwork")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Import data -------------------------------------------------------------

data <- read.csv(here("data",
                      "outputs",
                      "3g_final_dataset",
                      "pred_prey_sizes_DNAinteractions.csv"))

all_bs <- read.csv(here("data",
                     "outputs",
                     "3e_master_body_size_lists",
                     "all_mass_length.csv"))

all_pal <- all_bs %>%
  filter(Source %in% c("Ana_UG", "Palmyra"))

nodes <- read.csv(here("data", 
                       "raw_data",
                       "4_body_size_data",
                       "palmyra_nodes.csv"))

prey <- read.csv(here("data",
                      "outputs",
                      "3e_master_body_size_lists",
                      "prey_mass_length.csv"))

# Body Size distributions -------------------------------------------------

predator_sizes <- data %>%
  mutate(pred_mass_mg = 10^pred_log_mass_mg) %>%
  distinct(sample, pred_mass_mg) %>%
  rename(mass = pred_mass_mg) %>%
  mutate(type = "predator samples") %>%
  dplyr::select(mass, type)
  
interaction_sizes <- data %>%
  rename(mass = mean_prey_mass_mg) %>%
  mutate(type = "prey frequency")

prey_sizes <- data %>%
  distinct(Family, mean_prey_mass_mg) %>%
  rename(mass = mean_prey_mass_mg) %>%
  mutate(type = "prey identity")

community_sizes <- all_pal %>%
  group_by(Family) %>%
  summarise(mass = mean(Mass_mg, na.rm = T)) %>%
  mutate(type = "community") %>%
  filter(Family != "") %>%
  dplyr::select(-Family)

sizes <- predator_sizes %>%
  bind_rows(prey_sizes) %>%
  bind_rows(interaction_sizes) %>%
  bind_rows(community_sizes)

ggplot(sizes, aes(x = mass, fill = type))  +
  geom_density(alpha = 0.6, color = "black") +
  scale_x_log10() +
  theme_bw() +
  facet_wrap(~type, nrow = 4)

# Prey Family Size Dists --------------------------------------------------
prey <- prey %>%
  filter(Family != "") %>%
  filter(!is.na(Mass_mg))

prey_size_graph <- prey %>%
  filter(Family != "") %>%
  mutate(Family = as.factor(Family)) %>%
  mutate(Family = fct_reorder(Family, Mass_mg, .fun='sd')) %>%
  ggplot(aes(x = Family, y = Mass_mg)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(height=0, width=0.2),
              alpha = 0.4) +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

prey_family_count <- prey %>%
  filter(Family != "") %>%
  mutate(Family = as.factor(Family)) %>%
  mutate(Family = fct_reorder(Family, Mass_mg, .fun = 'sd')) %>%
  distinct(Genus, Family) %>%
  group_by(Family) %>%
  tally() %>%
  ggplot(aes(x = Family, y = n)) +
    geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

prey_size_graph / prey_family_count


