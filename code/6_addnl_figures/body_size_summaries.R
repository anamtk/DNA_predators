# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", "bipartite", "patchwork")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Import data -------------------------------------------------------------

prey <- read.csv(here("data",
                      "outputs",
                      "2g_master_body_size_lists",
                      "prey_mass_length.csv"))

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

