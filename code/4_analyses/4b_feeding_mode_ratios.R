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
                      "2i_final_dataset", 
                      "pred_prey_sizes_DNAinteractions.csv"))

size <- data %>%
  dplyr::select(-X, -reads) %>%
  mutate(pred_mass_mg = 10^(pred_log_mass_mg))

# Ratios by feeding interaction -------------------------------------------

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
  group_by(pred_class) %>%
  tally(name = "interactions")

ratios %>%
  distinct(sample_str, venom) %>%
  group_by(venom) %>%
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

# traits
# number of speices
# number of individuals
# number of interactions

ratios %>%
  mutate(pred_sz = 
           case_when(ratio > 1  ~ "predator larger",
                     ratio == 1 ~ "same size",
                     ratio < 1 ~ "prey larger")) %>%
  group_by(pred_sz) %>%
  tally()

predator_traits <- ratios %>%
  group_by(sample_str, hunting_mode, venom, webs) %>%
  tally(name = "interactions") %>%
  mutate(
    species = case_when(
      sample_str == "CEN" ~ "Mecistocephalus sp.",
      sample_str == "EUB" ~ "E. annulipes",
      sample_str == "HEV" ~ "H. venatoria",
      sample_str == "LRS" ~ "Opopaea sp.",
      sample_str == "NEO" ~ "N. theisi",
      sample_str == "PAN" ~ "P. flavescens",
      sample_str == "PHH" ~ "P. holdhausi", 
      sample_str == "SCY" ~ "S. longipes",
      sample_str == "SME" ~ "S. pallidus"
    )) %>%
  ungroup() %>%
  dplyr::select(species, interactions)

ratios <- ratios %>%
  mutate(pred_class = 
           case_when(sample_str %in% c("HEV", "LRS", "NEO", "SCY", "SME") ~ "Arachnida",
                     sample_str %in% c("EUB", "PAN", "PHH") ~ "Insecta",
                     sample_str == "CEN" ~ "Chilopoda"))

ratios %>%
  distinct(sample, sample_str, hunting_mode, venom, webs, pred_class) %>%
  group_by(sample_str, hunting_mode, venom, webs, pred_class) %>%
  tally(name = "samples") %>%
  mutate(
    species = case_when(
      sample_str == "CEN" ~ "Mecistocephalus sp.",
      sample_str == "EUB" ~ "E. annulipes",
      sample_str == "HEV" ~ "H. venatoria",
      sample_str == "LRS" ~ "Opopaea sp.",
      sample_str == "NEO" ~ "N. theisi",
      sample_str == "PAN" ~ "P. flavescens",
      sample_str == "PHH" ~ "P. holdhausi", 
      sample_str == "SCY" ~ "S. longipes",
      sample_str == "SME" ~ "S. pallidus"
    )) %>%
  ungroup() %>%
  dplyr::select(species, venom, webs, pred_class, samples) %>%
  rename(Class = pred_class) %>%
  left_join(predator_traits, by = "species") %>%
  gt() %>%
  tab_header(
    title = "Number of samples and interactions per species and traits") 


# Models ------------------------------------------------------------------

m_webs <- glmmTMB(log_ratio ~ webs + (1|sample) + (1|sample_str),
                data = ratios)

m_venom <- glmmTMB(log_ratio ~ venom + (1|sample) + (1|sample_str),
               data = ratios)

m_class <- glmmTMB(log_ratio ~ pred_class + (1|sample) + (1|sample_str),
                   data = ratios)

m_null <- glmmTMB(log_ratio ~ 1 + (1|sample) + (1|sample_str),
                  data = ratios)

AICc(m_webs, m_class, m_venom, m_null)

simulateResiduals(m_class, plot =T)

summary(m_class)

plot(allEffects(m_class))

pairs(emmeans(m_class, "pred_class"))

#f0f0f0
#bdbdbd
#636363
#"LRS", "SCY", "NEO", "CEN",
#"SME", "EUB", "PHH", "PAN",
#"HEV"

#[1] "#C70000" "#EA7700" "#EEB00C" 
#"#C68A2C" 
#"#89742F" 
#"#496C3C" "#158D8E" "#067D8D"
#[9] "#114C54"

ratios %>%
  mutate(sample_str = fct_reorder(sample_str, pred_mass_mg, .fun='mean')) %>%
  ggplot(aes(x = sample_str, y = ratio, color = sample_str, fill= sample_str)) +
  geom_boxplot(size = 0.75, alpha = 0.6) +
  geom_jitter(width = 0.25, height = 0, shape = 1) +
  theme_bw() +
  scale_color_manual(
    values = c("#C70000", "#EA7700", "#EEB00C", "#89742F", "#114C54",
               "#C68A2C",
               "#496C3C", "#158D8E", "#067D8D")) +
  scale_fill_manual(
    values = c("#C70000", "#EA7700", "#EEB00C", "#89742F", "#114C54",
               "#C68A2C",
               "#496C3C", "#158D8E", "#067D8D")) +
  scale_y_log10(breaks = c(0.01, 1, 100, 10000)) +
  labs(x = "Predator class", y = "Predator:prey mass ratio") +
  theme(axis.text = element_text(size =20),
        axis.title = element_text(size = 25)) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 1) +
  theme(axis.text = element_text(size =20),
        axis.title = element_text(size = 25),
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 15)) +
  facet_grid(.~pred_class, scales = "free_x", space = "free")

a <- AICc(m_class, m_null, m_venom, m_webs)

a %>% 
  mutate(delta = AICc - 850.5784) %>% 
  arrange(delta) %>%
  rownames_to_column(var = "model") %>%
  gt() %>% 
  fmt_number(
    columns = vars("df", "AICc", "delta"),
    decimals = 2) %>% 
  tab_header(
    title = "Model selection of predator trait models") 
