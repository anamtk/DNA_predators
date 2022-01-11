################
# Source for statistical analyses
# Ana Miller-ter Kuile
# June 2, 2021

# this is the source for all the data cleaning in prep for the
# statistical analyses in this paper

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 
                  "glmmTMB", "emmeans",
                  "MuMIn", "DHARMa",
                  "effects", "ggeffects",
                  "calecopal", "patchwork",
                  "emmeans", "gt",
                  "glmm")

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

# dataframe for pred-prey body size models:
size <- data %>%
  dplyr::select(-X, -reads) %>%
  mutate(pred_mass_mg = 10^(pred_log_mass_mg))

# dataframe for ratio analyses
ratios <- size %>%
  mutate(ratio = pred_mass_mg/mean_prey_mass_mg) %>%
  mutate(log_ratio = log10(ratio),
         log10_ratio = pred_log_mass_mg/mean_prey_log_mass_mg) %>%
  mutate(pred_class = 
           case_when(sample_str %in% c("HEV", "LRS", "NEO", "SCY", "SME") ~ "Arachnida",
                     sample_str %in% c("EUB", "PAN", "PHH") ~ "Insecta",
                     sample_str == "CEN" ~ "Chilopoda"))
