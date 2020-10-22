##########################
# 1. What are PPMR slopes at individual vs. species level -----
# Ana Miller-ter Kuile
# October 8, 2020
###########################

#basically the findings i want to write up 
#(that i think are in there) are “larger predator individuals 
#within/across species do not eat larger prey than smaller 
#individuals and the ratio of prey size to predator size for 
#this dataset follows similar relationships to other 
#published datasets (using a huge dataset I just found on the internet)”


# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 
                  "glmmTMB", "emmeans",
                  "MuMIn", "DHARMa",
                  "effects", "ggeffects")

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
  mutate(pred_mass_mg = exp(pred_log_mass_mg)) 


# Size ratio visualizations -----------------------------------------------

#Does the ratio of prey size to predator size vary with
#predator size?
#mean
size %>%
  mutate(mean_ratio = pred_mass_mg/mean_prey_mass_mg) %>%
  ggplot(aes(x = pred_mass_mg, y = mean_ratio, color = sample_str)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se =F) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() 

#min
size %>%
  mutate(min_ratio = pred_mass_mg/min_prey_mass_mg) %>%
  ggplot(aes(x = pred_mass_mg, y = min_ratio, color = sample_str)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se =F) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() 

