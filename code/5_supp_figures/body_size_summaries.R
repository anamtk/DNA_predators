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

#Sizes of interactions observed:
interaction_sizes <- data %>%
  rename(mass = mean_prey_mass_mg) %>%
  mutate(type = "Consumed prey size distribution")

#Sizes of prey species averaged
prey_sizes <- data %>%
  distinct(Family, mean_prey_mass_mg) %>%
  rename(mass = mean_prey_mass_mg) %>%
  mutate(type = "Prey size range")

#Sizes of all species averaged
community_sizes <- all_pal %>%
  group_by(Family) %>%
  summarise(mass = mean(Mass_mg, na.rm = T)) %>%
  mutate(type = "Community size range") %>%
  filter(Family != "") %>%
  dplyr::select(-Family)

prey_dat <- prey_sizes %>%
  #bind_rows(interaction_sizes) %>%
  bind_rows(community_sizes)

#f0f0f0
#bdbdbd
#636363
ggplot(prey_dat, aes(x = mass, fill = type)) +
  geom_histogram(color = "black") +
  scale_fill_manual(values = c("#bdbdbd", "#bdbdbd")) +
  scale_x_log10() +
  labs(x = "Mass (mg)", y = "Count") +
  theme_bw() +
  facet_wrap(~type, nrow = 3) +
  theme(strip.background = element_rect(fill="white"),
        axis.text = element_text(size =20),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 20),
        legend.position = "none")

ggplot(prey_dat, aes(x = mass, fill = type)) +
  geom_histogram(color = "black", 
                 position = "identity", 
                 alpha = 0.7,
                 bins = 20) +
  scale_fill_manual(values = c("#bdbdbd", "#636363")) +
  scale_x_log10() +
  labs(x = "Mass (mg)", y = "Count") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        axis.text = element_text(size =20),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = "bottom",
        legend.title = element_blank())

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


