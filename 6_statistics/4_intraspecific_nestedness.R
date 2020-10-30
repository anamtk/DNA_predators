##########################
# 3. Is there a nested heirarchy within species predicted
#by body size? -----
# Ana Miller-ter Kuile
# October 29, 2020
###########################

# this script analyzes whether within predator
#species, there is nestedness
#specifically, with smaller individuals
#eating a subset of what 
#larger individauls eat

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 
                  "glmmTMB", "emmeans",
                  "MuMIn", "DHARMa",
                  "effects", "ggeffects",
                  "vegan")

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

meta <- read.csv(here("data", "Sample_metadata.csv"))

meta <- meta %>%
  dplyr::select(Method, Island, Habitat,
                Microhabitat, Year, Date.Collected,
                Extraction.ID, ID) %>%
  distinct(Method, Island, Habitat,
           Microhabitat, Year, Date.Collected,
           Extraction.ID, ID)

size <- size %>%
  left_join(meta, by = c("sample" = "Extraction.ID"))

# Matrices by species -----------------------------------------------
mat_hev <- size %>%
  filter(sample_str == "HEV") %>%
  dplyr::select(sample, Family) %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = Family,
              values_from = presence,
              values_fill = 0) %>%
  column_to_rownames(var = "sample")

mat_phh <- size %>%
  filter(sample_str == "PHH") %>%
  dplyr::select(sample, Family) %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = Family,
              values_from = presence,
              values_fill = 0) %>%
  column_to_rownames(var = "sample")


# Nestedness calculations -------------------------------------------------

#using this info
#https://rdrr.io/rforge/vegan/man/nestedtemp.html
#nad the jaccard nestedness (turnover + nestedness) 
#as well as the NODF since it seems really common

## Matrix
out <- nestedbetajac(mat_phh)
out
out2 <- nestednodf(mat_phh)
out2
plot(out2)
## Use oecosimu to assess the non-randomness
oecosimu(mat_hev, nestednodf, "quasiswap")

oecosimu(mat_phh, nestedbetajac, "quasiswap")

