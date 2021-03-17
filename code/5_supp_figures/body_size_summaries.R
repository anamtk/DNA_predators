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

#average predator species sizes
pred_sp <- data %>%
  mutate(pred_mass = 10^pred_log_mass_mg) %>%
  distinct(sample, sample_str, pred_mass) %>%
  group_by(sample_str) %>%
  summarise(mass = mean(pred_mass)) %>%
  mutate(type = "Predator size range")

#Sizes of prey species averaged
prey_fams <- data %>%
  distinct(Family)

prey_sizes <- all_pal %>% 
  semi_join(prey_fams, by = "Family") %>%
  group_by(Family) %>%
  summarise(mass = mean(Mass_mg)) %>%
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
  bind_rows(community_sizes) %>%
  bind_rows(pred_sp)

#f0f0f0
#bdbdbd
#636363
prey_dat %>%
  mutate(type = factor(type, levels = c("Community size range",
                                        "Prey size range", 
                                        "Predator size range"))) %>%
ggplot(aes(x = mass, fill = type)) +
  geom_histogram(color = "black", 
                 position = "identity", 
                 alpha = 0.7,
                 bins = 20) +
  scale_fill_manual(values = c("#f0f0f0", "#bdbdbd", "#636363")) +
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

prey_fam_size <- prey %>%
  filter(Family != "") %>%
  group_by(Family) %>%
  summarise(mean_mass = mean(Mass_mg, na.rm =T),
            sd = sd(Mass_mg, na.rm=T),
            total = n(),
            se = sd/sqrt(total))

prey_family_sp <- prey %>%
  filter(Family != "") %>%
  distinct(Genus, Family) %>%
  group_by(Family) %>%
  tally(name = "number_of_species") %>%
  left_join(prey_fam_size, by = "Family")

prey_family_sp %>%
  mutate(number_of_species = as.factor(number_of_species)) %>%
ggplot(aes(x = number_of_species,
                           y = se)) +
  geom_boxplot() +
  scale_y_sqrt() +
  labs(x = "Number of species", 
       y = "Standard error of body mass within family") +
  theme_bw() 

