###########################
# Palmyra food web data cleaning
# Ana Miller-ter Kuile
# March 15, 2021
###########################

# this script cleans the Palmyra size and node datasets
# for importing into the master body size lists. 

###########################
# Load packages
package.list <- c("here", "tidyverse")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}
#############################

#palmyra raw body size data
pal <- read.csv(here("data",
                     "raw_data",
                     "4_body_size_data", 
                     "uncleaned_Palmyra",
                     "Palmyra_BS_Feb2021.csv"))

#remove NA and zero values
pal <- pal %>%
  filter(!is.na(Total_Mass_mg)) %>%
  filter(Total_Mass_mg != 0)

#filter out rows that have no information
pal <- pal %>%
  dplyr::select(!c("Femur_Vol_mm3","Femur_Vol_ml",          
                        "Femur_Mass_g", "Femur_Mass_mg",         
                        "Tibia_Vol_mm3", "Tibia_Vol_ml",         
                        "Tibia_Mass_g","Tibia_Mass_mg",
                        "Difference_Mass_Wieght",
                        "Femur.length.mm",       
                        "Femur.width.mm",
                        "Femur.height.mm",       
                        "Tibia.length.mm",
                        "Tibia.width.mm",        
                        "Tibia.height.mm"))


#Node names
pal_IDs <- read.csv(here("data",
                         "raw_data",
                         "4_body_size_data",
                         "uncleaned_Palmyra",
                         "Palmyra_nodenames.csv"))

#filter just adult so there is just one row per ID
unique(pal_IDs$Stage_Name)
pal_IDs <- pal_IDs %>%
  filter(Stage_Name %in% c("Adult", "Worker", "L_Adult")) %>%
  dplyr::select(!c("Body_Length_Mean_mm",
                   "Body_Mass_Mean_mg",    
                   "Body_Mass_SD_mg",
                   "Body_Mass_N", 
                   "Body_Mass_Source",
                   "Body_Mass_Proxy", 
                   "Body_Mass_Notes",
                   "Stage_Name")) 
  
  
# Export ------------------------------------------------------------------

write.csv(pal, here("data", 
                    "raw_data",
                    "4_body_size_data",
                    "palmyra_food_web_sizes.csv"))

write.csv(pal_IDs, here("data", 
                        "raw_data",
                        "4_body_size_data",
                        "palmyra_nodes.csv"))











