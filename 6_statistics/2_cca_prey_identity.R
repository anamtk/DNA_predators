##########################
# 3. Is diet composition determined by size or predator identity? -----
# Ana Miller-ter Kuile
# October 29, 2020
###########################

# this script analyzes whether predator identity
#size, or both determines the composition of
#prey in a given community, Sand PG Canopy

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 
                  "glmmTMB", "emmeans",
                  "MuMIn", "DHARMa",
                  "effects", "ggeffects",
                  "vegan", "remotes",
                  "eulerr", "dummies",
                  "calecopal")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Import data -------------------------------------------------------------

data2 <- read.csv(here("data", "outputs", "8_final_dataset",
                      "pred_prey_sizes_tp_DNAinteractions.csv"))

#select variables to remove
size2 <- data2 %>%
  dplyr::select(-X, -X.x, -X.y, -reads) %>%
  mutate(pred_mass_mg = exp(pred_log_mass_mg)) 

#metadata for samples
meta <- read.csv(here("data", "Sample_metadata.csv"))

#select variables of interest
meta <- meta %>%
  dplyr::select(Method, Island, Habitat,
                Microhabitat, Year, Date.Collected,
                Extraction.ID, ID) %>%
  distinct(Method, Island, Habitat,
           Microhabitat, Year, Date.Collected,
           Extraction.ID, ID)

#join to the interactions DF
size2 <- size2 %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  unite(environment, c("Habitat", "Microhabitat"), sep = "_", remove = F) 

# Data summary -------------------------------------------------------

size2 %>%
  distinct(sample, sample_str) %>%
  group_by(sample_str) %>%
  summarise(n = n())

# Matrix for CCA -----------------------------------------------

#matrix with diet as columns, samples as rows
mat <- size2 %>%
  dplyr::select(sample, Family) %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = Family,
              values_from = presence,
              values_fill = 0) %>% #fill missing values with 0's
  column_to_rownames(var = "sample")

#metadata
meta <- size2 %>%
  dplyr::select(sample, sample_str, Habitat, Microhabitat, pred_log_mass_mg) %>%
  distinct(sample, sample_str, Habitat, Microhabitat, pred_log_mass_mg) %>%
  column_to_rownames(var = "sample")

#check that rownames are the same:
all.equal(rownames(mat), rownames(meta))

# CCA -----------------------------------------------------------
#working off this tutorial
#http://dmcglinn.github.io/quant_methods/lessons/multivariate_models.html
# vegan requires that we write out each term if we are not going to 
# convert the factor to a dummy matrix 
# alternatively we could use a shorthand approach
cca <-  cca(mat ~ . , data=meta)

#view the CCA loadings
cca
#screeplot(cca_pg)
summary(cca)
#how much variation explained in this RDA
RsquareAdj(cca)
#ANOVA of whole model
anova(cca, permutations=10000)
#ANOVA of model terms
anova(cca, by='margin', permutations=10000)

#pretty ggplot plot 
#site metadata
sites <- meta %>%
  rownames_to_column("site")
#get the site (sample) scores out and attach to site metadata
CCAscores <- scores(cca, display = "sites") %>% 
  as.data.frame() %>% 
  rownames_to_column("site") %>%
  left_join(sites, by = "site")

#get the vectors out representing the loadings by species
CCAvect <- scores(cca, display = c("cn")) %>% 
  as.data.frame() %>%
  rownames_to_column(var = "ID") %>%
  filter(str_detect(ID, "sample_str")) %>%
  mutate(sample_str = str_sub(ID, -3)) %>%
  filter(sample_str %in% c("CEN", "HEV", "PAN", "PHH")) %>%
  mutate(pred_sp = ifelse(sample_str == "CEN", "Geophilomorpha sp.",
                          ifelse(sample_str == "HEV", "H. venatoria",
                                 ifelse(sample_str == "PAN", "P. flavescens",
                                        "P. holdhausi"))))

#get vector out representing the loading by body size
CCAsizevect <- scores(cca, display = "bp") %>% 
  as.data.frame() %>%
  rownames_to_column(var = "ID") %>%
  filter(ID == "pred_log_mass_mg")

pal_kelp <- cal_palette("kelp1", n = 9, type = "continuous")

pred_labels <- c("CEN" = "Geophilomorpha sp.", "EUB" = "E. annulipes", 
                 "HEV" = "H. venatoria", "LRS" = "Oonopidae sp.", 
                 "NEO" = "N. theisi", "PAN" = "P. flavescens",
                 "PHH" = "P. holdhausi", "SCY" = "S. longipes",
                 "SME" = "S. pallidus")

CCA_plot <- ggplot() +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_point(data = CCAscores, aes(x = CCA1, y = CCA2, color = sample_str), size = 3) +
  geom_segment(data = CCAvect, 
               aes(x = 0, y = 0, xend = CCA1, yend = CCA2), 
               arrow = arrow(length = unit(0.2, "cm")),
               size = 1) +
  geom_text(data = CCAvect, aes(x = CCA1, y = CCA2, label = pred_sp), 
            nudge_y = -0.1, nudge_x = -0.15, size = 5) +
  theme_bw() +
  scale_color_manual(values = pal_kelp,
                     labels = pred_labels) +
  labs(x = "CCA1 (14.4%)",
       y = "CCA2 (12.9%)") +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25))

# Euler plot --------------------------------------------------------------

#get variance by variable for Euler graph
#create environmental variables
species <- dummy(meta$sample_str)
hab <- dummy(meta$Habitat)
micro <- dummy(meta$Microhabitat)
bs <- meta$pred_log_mass_mg

#get variation by each
plot(varpart(mat, species, bs, hab, micro))
varpart(mat, species, bs, hab, micro)
#create vector of these variations

#Explanatory tables:
#X1:  species
#X2:  bs
#X3:  hab
#X4:  micro 

var <-  c("Predator_species" = 12, 
                    "Predator_mass" = 1, 
                    "Habitat" = 4, 
                    "Microhabitat" = 1, 
                    "Predator_species&Predator_mass" = 1,
                    "Predator_species&Habitat" = 1,
                    "Predator_species&Microhabitat" = 6,
                    "Predator_mass&Habitat" = 0,
                    "Predator_mass&Microhabitat" = 0,
                    "Habitat&Microhabitat" = 0,
                    "Predator_species&Predator_mass&Habitat" = 1,
                    "Predator_species&Predator_mass&Microhabitat" = 0,
                    "Predator_species&Habitat&Microhabitat" = 0,
                    "Predator_mass&Habitat&Microhabitat" = 0)
fit3 <- euler(var, 
              shape = "ellipse")

euler <- plot(fit3, 
     quantities = TRUE,
     fills = c("#f0f0f0", "#d9d9d9", "#bdbdbd", "#969696"))