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
  mutate(pred_mass_mg = exp(pred_log_mass_mg)) 

#############################
#Data explorations of expected important variables -------
#############################
ggplot(size, aes(x = pred_mass_mg, y = mean_prey_mass_mg, color = Tools)) +
  geom_abline(slope = 1, linetype = "dashed") +
  geom_point(size = 3) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  facet_wrap(~ sample_str)

ggplot(size, aes(x = pred_mass_mg, y = min_prey_mass_mg, color = Tools)) +
  geom_abline(slope = 1, linetype = "dashed") +
  geom_point(size = 3) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  facet_wrap(~ sample_str)
