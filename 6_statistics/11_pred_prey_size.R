###########################
# Load packages
package.list <- c("here", "tidyverse", "ggplot2", "bipartite")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}
#############################

#############################

pal <- read.csv(here("data", "outputs",
                     "5_rarefied_taxonomic_sort",
                     "fam_prey_DNA.csv"))

pal$sample <- str_sub(pal$sample, end=-2)

meta <- read.csv(here("data",
                      "Sample_metadata.csv"))

prey <- read.csv(here("data", "prey_sizes_pal_fw.csv"))

ana_prey <- read.csv(here("data", "pred_prey_sizes_ana_fw.csv"))

nodes <- read.csv(here("data", "Palmyra_Terrestrial_Node_List_19_June.csv"))

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

#write those to CSV
write.csv(pal_fams, here("data", "outputs", "8_prey_sizes",
                         "pal_DNA_prey_fams.csv"))

#anti-join to get those for which size not available yet
no_sizes <- pal_fams %>%
  anti_join(all_prey, by = "Family")

#write so that can look up some sizes in literature
write.csv(no_sizes, here("data", "outputs", "8_prey_sizes", 
                         "pal_prey_fam_nosize.csv"))

#import the new DF with size info included
node_litsize <- read.csv(here("data", "outputs", "8_prey_sizes",
                              "pal_prey_fam_litsize.csv"))

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


#combine pred and prey BODY SIZE DATA TO interaction DF
preds <- meta %>%
  dplyr::select(Extraction.ID, Length_mm) %>%
  rename("Pred_Length" = "Length_mm")

prey <- prey_fams %>%
  dplyr::select(Family, length_mean_mm) %>%
  rename("Prey_Length" = "length_mean_mm")

pred_prey <- pal %>%
  filter(reads > 0) %>%
  left_join(preds, by = c("sample" = "Extraction.ID")) %>%
  left_join(prey, by = "Family")

ggplot(pred_prey, aes(x = Pred_Length, y = Prey_Length, color = pred_ID)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", 
                                "#33a02c", "#fb9a99", "#e31a1c", 
                                "#fdbf6f", "#ff7f00", "#cab2d6")) +
  geom_smooth(method = "lm", se=F) +
  theme_bw()
 
 
ggplot(pred_prey, aes(x = Pred_Length, y = Prey_Length/Pred_Length, color = pred_ID)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", 
                                "#33a02c", "#fb9a99", "#e31a1c", 
                                "#fdbf6f", "#ff7f00", "#cab2d6")) +
  geom_smooth(method = "lm", se=F) +
  theme_bw()
