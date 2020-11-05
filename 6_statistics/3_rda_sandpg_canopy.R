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
                  "eulerr", "dummies")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}
#remotes::install_github("gavinsimpson/ggvegan")
library(ggvegan)
# Import data -------------------------------------------------------------

data2 <- read.csv(here("data", "outputs", "8_final_dataset",
                      "pred_prey_sizes_tp_DNAinteractions.csv"))

size2 <- data %>%
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

size2 <- size2 %>%
  left_join(meta, by = c("sample" = "Extraction.ID"))


# Data explorations -------------------------------------------------------

size2 %>%
  filter(Island == "Sand" & 
           Habitat == "PG" & 
           Microhabitat == "Canopy") %>%
  distinct(sample, sample_str, pred_log_mass_mg) %>%
  ggplot(aes(x = pred_log_mass_mg, fill = sample_str)) +
  geom_histogram(color = "black", position = "identity", alpha = 0.8) +
  theme_bw()

# Matrix for Sand PG Canopy -----------------------------------------------

#matrix
mat_pg <- size2 %>%
  filter(Island == "Sand" & 
           Habitat == "PG" & 
           Microhabitat == "Canopy") %>%
  dplyr::select(sample, Family) %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = Family,
              values_from = presence,
              values_fill = 0) %>%
  column_to_rownames(var = "sample")

#metadata for RDA
meta_pg <- size2 %>%
  filter(Island == "Sand" & 
           Habitat == "PG" & 
           Microhabitat == "Canopy") %>%
  dplyr::select(sample, sample_str, pred_log_mass_mg) %>%
  distinct(sample, sample_str, pred_log_mass_mg) %>%
  column_to_rownames(var = "sample")

#check that rownames are the same:
all.equal(rownames(mat_pg), rownames(meta_pg))

# RDA PG Canopy -----------------------------------------------------------
#working off this tutorial
#http://dmcglinn.github.io/quant_methods/lessons/multivariate_models.html
# vegan requires that we write out each term if we are not going to 
# convert the factor to a dummy matrix 
# alternatively we could use a shorthand approach
rda_pg <-  rda(mat_pg ~ . , data=meta_pg)

rda_pg

#how much variation explained in this RDA
RsquareAdj(rda_pg)
#ANOVA of whole model
anova(rda_pg, permutations=10000)
#ANOVA of model terms
anova(rda_pg, by='margin', permutations=10000)

#pretty ggplot plot (need to play with this)
autoplot <- autoplot(rda_pg, 
         layers = c("sites", "biplot", "centroids"),
         arrows = FALSE) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed")

#get variance by variable for Euler graph
#create environmental variables
species <- dummy(meta_pg$sample_str)
bs <- meta_pg$pred_log_mass_mg

#get variation by each
varpart(mat_pg, species, bs)

#create vector of these variations
var_pg <-  c("species" = 16, 
             "size&species" = 3, 
             "size" = 0)

#make a Euler and plot it
fit3 <- euler(var_pg, shape = "ellipse")
euler <- plot(fit3,
     quantities = TRUE)





