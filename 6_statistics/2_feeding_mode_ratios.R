##########################
# 2. Are prey:predator ratios explained by feeding mode? -----
# Ana Miller-ter Kuile
# October 29, 2020
###########################

# this script analyzes whether predator feeding
#mode (active vs. non-active) influences the
#size of prey they eat beyond just a 
#species-specific size mode. OR, maybe 
#there is somethign that is varyign at the 
#species level that matters instead, regardless
#of hunting mode. 

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

data <- read.csv(here("data", "outputs", "8_final_dataset",
                      "pred_prey_sizes_tp_DNAinteractions.csv"))

size <- data %>%
  dplyr::select(-X, -X.x, -X.y, -reads) %>%
  mutate(pred_mass_mg = 10^(pred_log_mass_mg))

# Ratios by feeding interaction -------------------------------------------
size %>%
  mutate(ratio = pred_mass_mg/mean_prey_mass_mg) %>%
  mutate(webs = ifelse(sample_str %in% c("CEN", "EUB", "PAN", "PHH"), 
                       "no", "yes")) %>%
  ggplot(aes(x = ratio, fill = webs)) +
  geom_histogram() +
  theme_bw() +
  facet_wrap(~webs) +
  scale_x_log10() +
  geom_vline(xintercept = 1) 

ratios <- size %>%
  mutate(ratio = pred_mass_mg/mean_prey_mass_mg) %>%
  mutate(webs = ifelse(sample_str %in% c("CEN", "EUB", "PAN", "PHH"), 
                       "No web use", "Web use")) %>%
  mutate(log_ratio = log10(ratio),
         log10_ratio = pred_log_mass_mg/mean_prey_log_mass_mg)

ratios %>%
  group_by(webs) %>%
  tally()

ratios %>%
  distinct(sample, webs) %>%
  group_by(webs) %>%
  tally()

ratios %>%
  group_by(webs) %>%
  summarise(mean_ratio = mean(log10_ratio),
            median = median(ratio),
            sd = sd(ratio),
            total = n(),
            se = sd/sqrt(total))

hist(ratios$ratio)
hist(ratios$log_ratio)
hist(ratios$log10_ratio)

m <- glmmTMB(log_ratio ~ webs + (1|sample_str),
             data = ratios)

summary(m)

fit <- simulateResiduals(m, plot =T)

plot(allEffects(m))

dredge(m)

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
