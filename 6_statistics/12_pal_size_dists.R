###########################
# Load packages
package.list <- c("here", "tidyverse", "ggplot2", "bipartite", "glmmTMB",
                  "MuMIn", "emmeans")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}
#############################

#############################
prey <- read.csv(here("data", "prey_sizes_pal_fw.csv"))

ana_prey <- read.csv(here("data", "pred_prey_sizes_ana_fw.csv"))

nodes <- read.csv(here("data", "Palmyra_Terrestrial_Node_List_19_June.csv"))

#import the new DF with size info included
node_litsize <- read.csv(here("data", "outputs", "8_prey_sizes",
                              "pal_prey_fam_litsize.csv"))

pal <- read.csv(here("data", "outputs",
                     "5_rarefied_taxonomic_sort",
                     "fam_prey_DNA_conservative.csv"))

pal$sample <- str_sub(pal$sample, end=-2)

#############################

#prey size data from FW project subset with important factors
prey <- prey %>%
  dplyr::select(Morphospecies, Stage_Name,
                Order.1, Family,
                Body_Length_Mean_mm) %>%
  filter(!is.na(Body_Length_Mean_mm)) %>%
  rename("ID" = "Morphospecies", 
         "Order" = "Order.1")

#data taken from my predators summarised by species
ana_prey1 <- ana_prey %>%
  filter(!is.na(Length)) %>%
  group_by(Order, Family, ID) %>%
  summarise(Body_Length_Mean_mm = mean(Length, na.rm=T)) %>%
  mutate(Stage_Name = "unknown")

#combine the two
all_prey <- prey %>%
  bind_rows(ana_prey1)

#get the distinct prey families from palmyra DNA data
pal_fams <- pal %>%
  distinct(Family)

node_litsize <- node_litsize %>%
  rowwise() %>% 
  mutate(Body_Length_Mean_mm = (mean(c(min_length, max_length), na.rm=T))) %>%
  dplyr::select(Family, Body_Length_Mean_mm) %>%
  mutate(ID = "",
         Order = "",
         Stage_Name = "")

#add to all prey
all_prey <- all_prey %>%
  bind_rows(node_litsize)

#summarise size by family to bind to the interaction DF
prey_fams <- all_prey %>%
  group_by(Order, Family) %>%
  summarise(length_mean_mm = mean(Body_Length_Mean_mm), 
            sd_length = sd(Body_Length_Mean_mm),
            n = n(), 
            se_length = sd_length/sqrt(n))


ggplot(prey_fams, aes(x = length_mean_mm)) +
  geom_histogram() +
  theme_bw() +
  scale_x_log10()
