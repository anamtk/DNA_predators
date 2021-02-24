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
                      "8_final_dataset", 
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
  tally()

ratios %>%
  group_by(venom) %>%
  tally()

ratios %>%
  group_by(hunting_mode) %>%
  tally()

ratios %>%
  distinct(sample, webs) %>%
  group_by(webs) %>%
  tally()

ratios %>%
  distinct(sample, venom) %>%
  group_by(venom) %>%
  tally()

ratios %>%
  distinct(sample, hunting_mode) %>%
  group_by(hunting_mode) %>%
  tally()

ratios %>%
  group_by(webs) %>%
  summarise(mean_ratio = mean(log10_ratio),
            median = median(ratio),
            sd = sd(ratio),
            total = n(),
            se = sd/sqrt(total))

ratios %>%
  group_by(venom) %>%
  summarise(mean_ratio = mean(log10_ratio),
            median = median(ratio),
            sd = sd(ratio),
            total = n(),
            se = sd/sqrt(total))

ratios %>%
  group_by(hunting_mode) %>%
  summarise(mean_ratio = mean(log10_ratio),
            median = median(ratio),
            sd = sd(ratio),
            total = n(),
            se = sd/sqrt(total))

hist(ratios$ratio)
hist(ratios$log_ratio)
hist(ratios$log10_ratio)

m_hunting_mode <- glmmTMB(log_ratio ~ hunting_mode + (1|sample_str),
             data = ratios)

m_webs <- glmmTMB(log_ratio ~ webs + (1|sample_str),
                data = ratios)

m_venom <- glmmTMB(log_ratio ~ venom + (1|sample_str),
               data = ratios)

AICc(m_hunting_mode, m_webs, m_venom)

summary(m_webs)

fit <- simulateResiduals(m_webs, plot =T)

plot(allEffects(m_webs))

ggplot(ratios, aes(x = webs, y = ratio)) +
  geom_boxplot(size = 0.75, fill = "#969696") +
  geom_jitter(width = 0.25, height = 0.25, shape = 1) +
  theme_bw() +
  scale_y_log10() +
  labs(x = "Web-using", y = "Predator:prey mass ratio") +
  theme(axis.text = element_text(size =20),
        axis.title = element_text(size = 25)) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 1) +
  theme(axis.text = element_text(size =20),
        axis.title = element_text(size = 25),
        axis.title.x = element_blank())


a <- AICc(m_webs, m_hunting_mode, m_venom)

a %>% 
  mutate(delta = AICc - 943.5464) %>% 
  rownames_to_column(var = "model") %>%
  gt() %>% 
  fmt_number(
    columns = vars("df", "AICc", "delta"),
    decimals = 2) %>% 
  tab_header(
    title = "Model selection of predator trait models") 
