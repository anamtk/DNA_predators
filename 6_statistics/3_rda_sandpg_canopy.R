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

# Matrix for Sand PG Canopy -----------------------------------------------

#matrix
mat_pg <- size %>%
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
meta_pg <- size %>%
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

RsquareAdj(rda_pg)

plot(rda_pg, type='n', scaling=1)
orditorp(rda_pg, display='sp', cex=0.5, scaling=1, col='blue')
text(rda_pg, display='cn', col='red')

anova(rda_pg, permutations=10000)
anova(rda_pg, by='margin', permutations=10000)





pca <- rda(dune)
autoplot(rda_pg) +
  theme_bw()

autoplot.rda(rda_pg)
## Just the species scores
autoplot(rda_pg, 
         layers = c("sites", "biplot", "centroids"),
         arrows = FALSE) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed")
?autoplot

#get variance by 
species <- dummy(meta_pg$sample_str)
bs <- meta_pg$pred_log_mass_mg

varpart(mat_pg, species, bs)

var_pg <-  c("species" = 16, 
             "size&species" = 3, 
             "size" = 0)

fit3 <- euler(var_pg, shape = "ellipse")
plot(fit3,
     quantities = TRUE)




