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

pal <- read.csv(here("data", "outputs",
                     "5_rarefied_taxonomic_sort",
                     "fam_prey_DNA_conservative.csv"))

pal$sample <- str_sub(pal$sample, end=-2)

meta <- read.csv(here("data",
                      "Sample_metadata.csv"))

meta <- meta %>%
  group_by(ID, Extraction.ID) %>%
  summarise(Length_mm = mean(Length_mm))

prey <- read.csv(here("data", "prey_sizes_pal_fw.csv"))

ana_prey <- read.csv(here("data", "pred_prey_sizes_ana_fw.csv"))

nodes <- read.csv(here("data", "Palmyra_Terrestrial_Node_List_19_June.csv"))

pred_mode <- read.csv(here("data", "Predator_IDs.csv"))

pred_mode <- pred_mode %>%
  dplyr::select(pred_ID, Feeding_mode)

tp <- read.csv(here("data", "outputs", "7_all_webs",
                    "diet_families_tp.csv"))

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
  ungroup() %>%
  dplyr::select(Extraction.ID, Length_mm) %>%
  rename("Pred_Length" = "Length_mm")

prey <- prey_fams %>%
  dplyr::select(Family, length_mean_mm) %>%
  rename("Prey_Length" = "length_mean_mm")

pred_prey <- pal %>%
  filter(reads > 0) %>%
  left_join(preds, by = c("sample" = "Extraction.ID")) %>%
  left_join(prey, by = "Family")

pred_prey <- pred_prey %>%
  filter(pred_ID != "Euborellia annulipes") %>%
  filter(Prey_Length < 73)

pred_prey <- pred_prey %>%
  mutate(ratio = Pred_Length/Prey_Length)

m1 <- glmmTMB(ratio ~ Feeding_mode + (1|pred_ID),
              data=pred_prey)
summary(m1)

ggplot(pred_prey, aes(x = log(Pred_Length), y = log(Prey_Length), color = Feeding_mode)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", 
                                "#33a02c", "#fb9a99", "#e31a1c", 
                                "#fdbf6f", "#ff7f00")) +
  theme_bw() +
  facet_wrap(~pred_ID)

ggplot(pred_prey, aes(x = log(Pred_Length), y = log(Prey_Length), color = Feeding_mode)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#a6cee3", "#1f78b4")) +
  theme_bw() +
  facet_wrap(~Feeding_mode)

ggplot(pred_prey, aes(x = log(Pred_Length), y = log(Prey_Length), color = Feeding_mode)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#a6cee3", "#1f78b4")) +
  geom_smooth(method = "lm", se =F) +
  theme_bw()

ggplot(pred_prey, aes(x = Feeding_mode, y = Pred_Length/Prey_Length, fill = Feeding_mode)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = c("#a6cee3", "#1f78b4")) +
  theme_bw() +
  scale_y_log10()
 
ggplot(pred_prey, aes(x = Pred_Length, fill = Feeding_mode)) +
  geom_histogram(position = "dodge", binwidth = 1) +theme_bw()

overlap <- pred_prey %>%
  filter(Pred_Length >= 2 & Pred_Length <= 12)

ggplot(overlap, aes(x = Feeding_mode, y = Pred_Length/Prey_Length, fill = Feeding_mode)) +
  geom_boxplot() + 
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = c("#a6cee3", "#1f78b4")) +
  theme_bw() +
  scale_y_log10()

m2 <- glmmTMB(ratio ~ Feeding_mode + (1|pred_ID),
              data=overlap)
summary(m2)

#Species Level
pred_sp <- ana_prey1 %>%
  rename("Pred_Length" = "Body_Length_Mean_mm") %>%
  ungroup() %>%
  dplyr::select(ID, Pred_Length)

sp_pred <- pal %>%
  left_join(pred_sp, by = c("pred_ID" = "ID"))

sp_pred_prey <- sp_pred %>%
  left_join(prey, by = "Family") %>%
  filter(reads > 0) %>%
  filter(pred_ID != "Euborellia annulipes") %>%
  filter(Prey_Length < 73)

#interaction frequency
freq <- sp_pred_prey %>%
  group_by(pred_ID, Family, Pred_Length, Prey_Length) %>%
  summarise(frequency = n()) %>%
  left_join(pred_mode, by = "pred_ID")

ggplot(freq, aes(x = Pred_Length/Prey_Length, y = frequency, color = Feeding_mode)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_color_manual(values = c("#a6cee3", "#1f78b4")) +
  theme_bw() +
  scale_x_log10() +
  facet_wrap(~Feeding_mode)

freq_overlap <- freq %>%
  filter(Pred_Length >= 2 & Pred_Length <= 12)

ggplot(freq_overlap, aes(x = Pred_Length/Prey_Length, y = frequency, color = Feeding_mode)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_color_manual(values = c("#a6cee3", "#1f78b4")) +
  theme_bw() +
  scale_x_log10() +
  facet_wrap(~Feeding_mode)




