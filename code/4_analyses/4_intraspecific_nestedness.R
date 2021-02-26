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
                  "vegan", "iNEXT")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Import data -------------------------------------------------------------

data <- read.csv(here("data", 
                      "outputs",  
                      "8_final_dataset", 
                      "pred_prey_sizes_DNAinteractions.csv"))

size3 <- data %>%
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

size3 <- size3 %>%
  left_join(meta, by = c("sample" = "Extraction.ID"))


# Species accumulation curves ---------------------------------------------
#pal_kelp <- cal_palette("kelp1", n = 9, type = "continuous")
#pal_kelp
#LRS: "#C70000" 
# SCY: "#EA7700" 
#NEO: "#EEB00C" 
#CEN: "#C68A2C" 
#SME: "#89742F" 
#EUB: "#496C3C" 
#PHH: "#158D8E"
#PAN: "#067D8D" 
#HEV: "#114C54"

sample_sizes <- size3 %>%
  filter(sample_str %in% c("HEV", "PHH", "NEO")) %>%
  distinct(sample, sample_str) %>%
  group_by(sample_str) %>%
  summarise(frequency = n()) %>%
  mutate(Family = "Total")
  
frequencies <- size3 %>%
  filter(sample_str %in% c ("HEV", "PHH", "NEO")) %>%
  group_by(sample_str, Family) %>%
  summarise(frequency = n()) 

hev_freq <- frequencies %>%
  filter(sample_str == "HEV") %>%
  ungroup() %>%
  dplyr::select(-sample_str)
  
hev_freq <- rbind(data.frame(Family = "total", frequency = 53), hev_freq) 

hev_freq <- hev_freq %>%
  dplyr::select(-Family) %>%
  as_vector()

out_hev <- iNEXT(hev_freq, q=0, datatype = "incidence_freq")
out_hev$AsyEst
hev_acc <- ggiNEXT(out_hev, type=1, color.var="site") + 
  theme_bw() + 
  theme(legend.position="none",
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 25)) +
  scale_color_manual(values = c("#114C54")) +
  scale_fill_manual(values = c("#114C54")) +
  geom_ribbon(aes(ymin=42-11, ymax=42+11), color = NA, alpha=0.1) +
  geom_hline(yintercept = 42, linetype = "dashed", size =1) 

neo_freq <- frequencies %>%
  filter(sample_str == "NEO") %>%
  ungroup() %>%
  dplyr::select(-sample_str)

neo_freq <- rbind(data.frame(Family = "total", frequency = 24), neo_freq) 

neo_freq <- neo_freq %>%
  dplyr::select(-Family) %>%
  as_vector()

out_neo <- iNEXT(neo_freq, q=0, datatype = "incidence_freq")
out_neo$AsyEst
neo_acc <- ggiNEXT(out_neo, type=1, color.var="site") + 
  theme_bw() + 
  theme(legend.position="none",
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 25)) +
  scale_color_manual(values = c("#EEB00C")) +
  scale_fill_manual(values = c("#EEB00C")) +
  geom_ribbon(aes(ymin=44-11, ymax=44+11), color = NA, alpha=0.1) +
  geom_hline(yintercept = 44, linetype = "dashed", size =1) 

phh_freq <- frequencies %>%
  filter(sample_str == "PHH") %>%
  ungroup() %>%
  dplyr::select(-sample_str)

phh_freq <- rbind(data.frame(Family = "total", frequency = 42), phh_freq) 

phh_freq <- phh_freq %>%
  dplyr::select(-Family) %>%
  as_vector()

out_phh <- iNEXT(phh_freq, q=0, datatype = "incidence_freq")
out_phh$AsyEst
phh_acc <- ggiNEXT(out_phh, type=1, color.var="site") + 
  theme_bw() + 
  theme(legend.position="none",
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 25)) +
  scale_color_manual(values = c("#158D8E")) +
  scale_fill_manual(values = c("#158D8E")) +
  geom_ribbon(aes(ymin=36-17, ymax=36+17), color = NA, alpha=0.1) +
  geom_hline(yintercept = 36, linetype = "dashed", size =1) 
# Matrices by species -----------------------------------------------

mat_hev <- size3 %>%
  filter(sample_str == "HEV") %>%
  dplyr::select(sample, Family) %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = Family,
              values_from = presence,
              values_fill = 0) %>%
  column_to_rownames(var = "sample")

mat_neo <- size3 %>%
  filter(sample_str == "NEO") %>%
  dplyr::select(sample, Family) %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = Family,
              values_from = presence,
              values_fill = 0) %>%
  column_to_rownames(var = "sample")

mat_phh <- size3 %>%
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
out_hev <- nestedbetajac(mat_hev)
#out_hev
out_hev2 <- nestednodf(mat_hev)
#out_hev2
#plot(out_hev2)
## Use oecosimu to assess the non-randomness
oecosimu(mat_hev, nestednodf, "quasiswap")
oecosimu(mat_hev, nestedbetajac, "quasiswap")

## Matrix
out_neo <- nestedbetajac(mat_neo)
#out_neo
out_neo2 <- nestednodf(mat_neo)
#out_neo2
#plot(out_neo2)
## Use oecosimu to assess the non-randomness
oecosimu(mat_neo, nestednodf, "quasiswap")
oecosimu(mat_neo, nestedbetajac, "quasiswap")

## Matrix
out_phh <- nestedbetajac(mat_phh)
#out_phh
out_phh2 <- nestednodf(mat_phh)
#out_phh2
#plot(out_phh2)
## Use oecosimu to assess the non-randomness
oecosimu(mat_phh, nestednodf, "quasiswap")
oecosimu(mat_phh, nestedbetajac, "quasiswap")

# Plot lack of nestedness -------------------------------------------------

hev_matrix <- size3 %>%
  filter(sample_str == "HEV") %>%
  mutate(pred_mass_mg = as.factor(pred_mass_mg)) %>%
  group_by(pred_mass_mg, Family) %>%
  summarise(total = n()) %>%
  pivot_wider(names_from = pred_mass_mg,
              values_from = total,
              values_fill = 0) %>%
  mutate(total = rowSums(.[-1])) %>%
  arrange(total) %>%
  mutate(Family = factor(Family, levels = Family)) %>%
  dplyr::select(-total) %>%
  pivot_longer(!Family,
               names_to = "pred_mass_mg",
               values_to = "total") %>%
  mutate(total = as.factor(total)) %>%
  mutate(pred_mass_mg = as.numeric(pred_mass_mg),
         pred_mass_mg = as.factor(pred_mass_mg)) %>%
  ggplot(aes(x = pred_mass_mg, y = Family, fill = total)) +
  geom_tile() +
  theme_bw() +
  scale_fill_manual(values = c("#FFFFFF", "#114C54", "#114C54")) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 25)) +
  labs(x = "Increasing predator size", y = "Prey family")

neo_matrix <- size3 %>%
  filter(sample_str == "NEO") %>%
  mutate(pred_mass_mg = as.factor(pred_mass_mg)) %>%
  group_by(pred_mass_mg, Family) %>%
  summarise(total = n()) %>%
  pivot_wider(names_from = pred_mass_mg,
              values_from = total,
              values_fill = 0) %>%
  mutate(total = rowSums(.[-1])) %>%
  arrange(total) %>%
  mutate(Family = factor(Family, levels = Family)) %>%
  dplyr::select(-total) %>%
  pivot_longer(!Family,
               names_to = "pred_mass_mg",
               values_to = "total") %>%
  mutate(total = as.factor(total)) %>%
  mutate(pred_mass_mg = as.numeric(pred_mass_mg),
         pred_mass_mg = as.factor(pred_mass_mg)) %>%
  ggplot(aes(x = pred_mass_mg, y = Family, fill = total)) +
  geom_tile() +
  theme_bw() +
  scale_fill_manual(values = c("#FFFFFF", "#EEB00C", "#EEB00C")) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 25)) +
  labs(x = "Increaseing predator size", y = "Prey family")

phh_matrix <- size3 %>%
  filter(sample_str == "PHH") %>%
  mutate(pred_mass_mg = as.factor(pred_mass_mg)) %>%
  group_by(pred_mass_mg, Family) %>%
  summarise(total = n()) %>%
  pivot_wider(names_from = pred_mass_mg,
              values_from = total,
              values_fill = 0) %>%
  mutate(total = rowSums(.[-1])) %>%
  arrange(total) %>%
  mutate(Family = factor(Family, levels = Family)) %>%
  dplyr::select(-total) %>%
  pivot_longer(!Family,
               names_to = "pred_mass_mg",
               values_to = "total") %>%
  mutate(total = as.factor(total)) %>%
  mutate(pred_mass_mg = as.numeric(pred_mass_mg),
         pred_mass_mg = as.factor(pred_mass_mg)) %>%
  ggplot(aes(x = pred_mass_mg, y = Family, fill = total)) +
  geom_tile() +
  theme_bw() +
  scale_fill_manual(values = c("#FFFFFF", "#158D8E", "#158D8E")) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 25)) +
  labs(x = "Increasing predator size", y = "Prey family")



