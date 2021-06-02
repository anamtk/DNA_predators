# Community body size distributions
# Ana Miller-ter Kuile
# June 2, 2021


# Load data -----------------------------------------------------------

library(here)

source(here("code", 
            "4_analyses",
            "4_source.R"))

all_bs <- read.csv(here("data",
                        "outputs",
                        "3g_master_body_size_lists",
                        "all_mass_length.csv"))

all_pal <- all_bs %>%
  filter(Source %in% c("Ana_UG", "Palmyra"))

nodes <- read.csv(here("data", 
                       "raw_data",
                       "4_body_size_data",
                       "palmyra_nodes.csv"))

prey <- read.csv(here("data",
                      "outputs",
                      "3g_master_body_size_lists",
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

# Figure ------------------------------------------------------------------

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

