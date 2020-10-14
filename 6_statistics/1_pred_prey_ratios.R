##########################
# 1. What are predator-prey ratios? -----
# Ana Miller-ter Kuile
# October 8, 2020
###########################

# this script analyzes predator-prey body size ratios, 
#answering the question:
# 1. What is pred:prey body size ratio and does it vary by species,
# body size, or feeding mode (tools, or more broad webs-no webs)
# note: may have to think about foraging mode a bit more in terms of
# active, web, sit-and-wait as I could see these having different
# metabolic cost-benefits that would determine ratios? maybe?

###########################
# Load packages-----
package.list <- c("here", "tidyverse", 
                  "glmmTMB", "emmeans",
                  "MuMIn", "DHARMa",
                  "effects", "ggeffects")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}
#############################

#############################
#Import Data -------
#############################
data <- read.csv(here("data", "outputs", "8_final_dataset",
                      "pred_prey_sizes_tp_DNAinteractions.csv"))

size <- data %>%
  dplyr::select(-X, -X.1, -X.x, -ASV, -ID_level, -X.y) %>%
  mutate(pred_mass_mg = exp(pred_log_mass_mg),
         mean_prey_pred_ratio = mean_prey_mass_mg/pred_mass_mg,
         min_prey_pred_ratio = min_prey_mass_mg/pred_mass_mg) 

#############################
#Data explorations of expected important variables -------
#############################
ggplot(size, aes(x = pred_mass_mg, y = mean_prey_mass_mg, color = sample_str)) +
  geom_abline(slope = 1, linetype = "dashed") +
  geom_point(size = 3) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  facet_wrap(~ Tools)

ggplot(size, aes(x = pred_mass_mg, y = mean_prey_mass_mg/pred_mass_mg, color = sample_str)) +
  geom_abline(slope = 1, linetype = "dashed") +
  geom_point(size = 3) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  facet_wrap(~ Feeding_mode)

ggplot(size, aes(x = pred_mass_mg, y = min_prey_mass_mg, color = sample_str)) +
  geom_abline(slope = 1, linetype = "dashed") +
  geom_point(size = 3) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  facet_wrap(~ Tools)

ggplot(size, aes(x = pred_mass_mg, y = min_prey_mass_mg, color = sample_str)) +
  geom_abline(slope = 1, linetype = "dashed") +
  geom_point(size = 3) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  facet_wrap(~ Feeding_mode)

ggplot(size, aes(x = pred_mass_mg, y = min_prey_mass_mg/pred_mass_mg, color = sample_str)) +
  geom_abline(slope = 1, linetype = "dashed") +
  geom_point(size = 3) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  facet_wrap(~ Feeding_mode)

#############################
#Basic model of mean -------
#############################

m1 <- glmmTMB(mean_prey_pred_ratio ~ pred_log_mass_mg*Feeding_mode + (1|sample_str),
              data = size)

summary(m1)

dredge(m1)

fit <- simulateResiduals(m1, plot = TRUE)

#############################
#log transformed model -------
#############################

m2 <- glmmTMB(min_prey_log_mass_mg ~ pred_log_mass_mg*sample_str,
              data = size)
summary(m2)
dredge(m2)

m3 <- glmmTMB(min_prey_log_mass_mg ~ pred_log_mass_mg + sample_str,
              data = size)

fit <- simulateResiduals(m3, plot = TRUE)

hist(log(size$mean_prey_pred_ratio))
hist(log(size$min_prey_pred_ratio))
hist(size$mean_prey_mass_mg)
hist(size$min_prey_mass_mg)
