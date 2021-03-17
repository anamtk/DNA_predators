##########################
# 2. Are prey:predator ratios explained by feeding mode? -----
# Ana Miller-ter Kuile
# October 29, 2020
###########################

# this script analyzes whether a set of predator traits including 
# hunting mode, web use and venom use, influence prey size

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 
                  "glmmTMB", "emmeans",
                  "MuMIn", "DHARMa",
                  "effects", "ggeffects",
                  "calecopal", "patchwork",
                  "emmeans", "gt")

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

size <- data %>%
  dplyr::select(-X, -reads) %>%
  mutate(pred_mass_mg = 10^(pred_log_mass_mg))

# Ratios by feeding interaction -------------------------------------------
size %>%
  mutate(ratio = pred_mass_mg/mean_prey_mass_mg) %>%
  ggplot(aes(x = ratio, fill = webs)) +
  geom_histogram() +
  theme_bw() +
  facet_wrap(~webs) +
  scale_x_log10() +
  geom_vline(xintercept = 1) 

size %>%
  mutate(ratio = pred_mass_mg/mean_prey_mass_mg) %>%
  ggplot(aes(x = ratio, fill = venom)) +
  geom_histogram() +
  theme_bw() +
  facet_wrap(~venom) +
  scale_x_log10() +
  geom_vline(xintercept = 1) 

size %>%
  mutate(ratio = pred_mass_mg/mean_prey_mass_mg) %>%
  ggplot(aes(x = ratio, fill = hunting_mode)) +
  geom_histogram() +
  theme_bw() +
  facet_wrap(~hunting_mode) +
  scale_x_log10() +
  geom_vline(xintercept = 1) 

ratios <- size %>%
  mutate(ratio = pred_mass_mg/mean_prey_mass_mg) %>%
  mutate(log_ratio = log10(ratio),
         log10_ratio = pred_log_mass_mg/mean_prey_log_mass_mg)

ratios %>%
  group_by(webs) %>%
  tally(name = "interactions")

ratios %>%
  distinct(sample_str, webs) %>%
  group_by(webs) %>%
  tally(name = "species")

ratios %>%
  group_by(venom) %>%
  tally(name = "interactions")

ratios %>%
  distinct(sample_str, venom) %>%
  group_by(venom) %>%
  tally(name = "species")

ratios %>%
  group_by(hunting_mode) %>%
  tally(name = "interactions")

ratios %>%
  distinct(sample_str, hunting_mode) %>%
  group_by(hunting_mode) %>%
  tally(name = "species")

ratios %>%
  distinct(sample, webs) %>%
  group_by(webs) %>%
  tally(name = "individuals")

ratios %>%
  distinct(sample, venom) %>%
  group_by(venom) %>%
  tally(name = "individuals")

ratios %>%
  distinct(sample, hunting_mode) %>%
  group_by(hunting_mode) %>%
  tally(name = "individuals")

ratios %>%
  group_by(webs) %>%
  summarise(mean_ratio = mean(log10_ratio, na.rm = T),
            median = median(ratio, na.rm = T),
            sd = sd(ratio, na.rm = T),
            total = n(),
            se = sd/sqrt(total))

ratios %>%
  group_by(venom) %>%
  summarise(mean_ratio = mean(log10_ratio, na.rm = T),
            median = median(ratio, na.rm = T),
            sd = sd(ratio, na.rm = T),
            total = n(),
            se = sd/sqrt(total))

ratios %>%
  group_by(hunting_mode) %>%
  summarise(mean_ratio = mean(log10_ratio, na.rm = T),
            median = median(ratio, na.rm = T),
            sd = sd(ratio, na.rm = T),
            total = n(),
            se = sd/sqrt(total))

hist(ratios$ratio)
hist(ratios$log_ratio)
hist(ratios$log10_ratio)

# traits
# number of speices
# number of individuals
# number of interactions

predator_traits <- ratios %>%
  group_by(sample_str, hunting_mode, venom, webs) %>%
  tally(name = "interactions") %>%
  mutate(
    species = case_when(
      sample_str == "CEN" ~ "Geophilomorpha sp",
      sample_str == "EUB" ~ "E. annulipes",
      sample_str == "HEV" ~ "H. venatoria",
      sample_str == "LRS" ~ "Oonopidae sp",
      sample_str == "NEO" ~ "N. theisi",
      sample_str == "PAN" ~ "P. flavescens",
      sample_str == "PHH" ~ "P. holdhausi", 
      sample_str == "SCY" ~ "S. longipes",
      sample_str == "SME" ~ "S. pallidus"
    )) %>%
  ungroup() %>%
  dplyr::select(species, interactions)

ratios %>%
  distinct(sample, sample_str, hunting_mode, venom, webs) %>%
  group_by(sample_str, hunting_mode, venom, webs) %>%
  tally(name = "samples") %>%
  mutate(
    species = case_when(
      sample_str == "CEN" ~ "Geophilomorpha sp",
      sample_str == "EUB" ~ "E. annulipes",
      sample_str == "HEV" ~ "H. venatoria",
      sample_str == "LRS" ~ "Oonopidae sp",
      sample_str == "NEO" ~ "N. theisi",
      sample_str == "PAN" ~ "P. flavescens",
      sample_str == "PHH" ~ "P. holdhausi", 
      sample_str == "SCY" ~ "S. longipes",
      sample_str == "SME" ~ "S. pallidus"
    )) %>%
  ungroup() %>%
  dplyr::select(species, hunting_mode, venom, webs, samples) %>%
  left_join(predator_traits, by = "species") %>%
  gt() %>%
  tab_header(
    title = "Number of samples and interactions per species and traits") 

# Models ------------------------------------------------------------------

m_hunting_mode <- glmmTMB(log_ratio ~ hunting_mode + (1|sample) + (1|sample_str),
             data = ratios)

m_webs <- glmmTMB(log_ratio ~ webs +  (1|sample) + (1|sample_str),
                data = ratios)

m_venom <- glmmTMB(log_ratio ~ venom + (1|sample) + (1|sample_str),
               data = ratios)

m_null <- glmmTMB(log_ratio ~ 1 + (1|sample) + (1|sample_str),
                  data = ratios)

AICc(m_hunting_mode, m_webs, m_venom, m_null)

simulateResiduals(m_webs, plot =T)

summary(m_webs)

#f0f0f0
#bdbdbd
#636363

ratios %>%
  summarise(max = min(ratio))
ratios %>%
  mutate(sample_str = fct_reorder(sample_str, pred_mass_mg, .fun='mean')) %>%
  ggplot(aes(x = sample_str, y = ratio, fill = webs)) +
  geom_boxplot(size = 0.75) +
  geom_jitter(width = 0.25, height = 0, shape = 1) +
  theme_bw() +
  scale_fill_manual(labels = c("No web use", "Web use"), values = c("#f0f0f0", "#636363")) +
  scale_y_log10(breaks = c(0.01, 1, 100, 10000)) +
  labs(x = "Web-using", y = "Predator:prey mass ratio") +
  theme(axis.text = element_text(size =20),
        axis.title = element_text(size = 25)) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 1) +
  theme(axis.text = element_text(size =20),
        axis.title = element_text(size = 25),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 20))

a <- AICc(m_webs, m_null, m_hunting_mode, m_venom)

a %>% 
  mutate(delta = AICc - 939.1664) %>% 
  rownames_to_column(var = "model") %>%
  gt() %>% 
  fmt_number(
    columns = vars("df", "AICc", "delta"),
    decimals = 2) %>% 
  tab_header(
    title = "Model selection of predator trait models") 
