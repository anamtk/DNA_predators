##########################
# Palmyra body size measurements progress -----
# Ana Miller-ter Kuile
# February 26, 2021
###########################

# this script looks at the most updated version of the palmyra body size 
# measurement dataset and subsets the families/species still needed

###########################
# Load packages-----
package.list <- c("here", "tidyverse")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}
#############################


# Import data -------------------------------------------------------------


data <- read.csv(here("data", 
                      "outputs", 
                      "3c_rarefied_taxonomic_sort", 
                      "fam_prey_DNA_conservative.csv"))

prey <- data %>%
  distinct(Domain, Phylum, Class, Order, Family)

prey_fams <- prey %>%
  dplyr::select(Family)

sizes <- read.csv(here("data", 
                       "raw_data",
                       "4_body_size_data", 
                       "Palmyra_BS_Feb2021.csv"))
str(sizes)

unique(sizes$Family)

sizes <- sizes %>%
  dplyr::select(Morphospecies, 
                Phylum, 
                Class, 
                Order.1, 
                Family,
                Body_Length_Mean_mm,
                Body_Mass_Mean_mg,
                Body_Mass_SD_mg,
                Body_Mass_N)

missing <- sizes %>%
  filter(is.na(Body_Mass_Mean_mg)) %>%
  distinct(Morphospecies, 
           Phylum, 
           Class, 
           Order.1, 
           Family,
           Body_Length_Mean_mm,
           Body_Mass_Mean_mg,
           Body_Mass_SD_mg,
           Body_Mass_N) %>%
  filter(Phylum != "")
str(missing)

missing_sp <- missing %>%
  semi_join(prey_fams, by = "Family")

write.csv(missing_sp, here("data", 
                           "outputs", 
                           "3g_palmyra_body_size_measurements",
                           "missing_palmyra_datafeb162021.csv"))





