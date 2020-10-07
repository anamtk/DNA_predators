##########################
# Predator mass-length relationships
# Ana Miller-ter Kuile
# October 7, 2020
###########################

# this script builds species-specific mass-length
#relationships for each predator species in my 
#samples so I can convert lengths to masses for analyses

###########################
# Load packages
package.list <- c("here", "tidyverse")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}
#############################

#############################
# Load and combine datasets --------
#############################

pred_size <- read.csv(here("data", "outputs", "8_prey_sizes", "pred_mass_length.csv"))

pred_size <- pred_size %>%
  mutate(Family = ifelse(Order == "Orthoptera", "Tettigonidae", Family)) %>%
  mutate(Family = ifelse(Family == "Tettigonidae", "Tettigoniidae", Family))

pred_id <- read.csv(here("data", "Predator_IDs.csv"))

pred_id <- pred_id %>%
  dplyr::select(pred_Family, sample_str)

pred_size <- pred_size %>%
  left_join(pred_id, by = c("Family" = "pred_Family")) %>%
  mutate(sample_str = ifelse(Order == "Geophilomorpha", "CEN", sample_str))


#############################
# Vis by species --------
#############################
ggplot(pred_size, aes(x = Length_mm, y = Mass_mg, color = sample_str)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~sample_str)

ggplot(pred_size, aes(x = Length_mm, y = Mass_mg, color = sample_str)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_x_log10() +
  scale_y_log10() 


