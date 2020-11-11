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

size2 %>%
  filter(Island == "Sand" & 
           Habitat == "PG" & 
           Microhabitat == "Canopy") %>%
  filter(sample_str != "PHH") %>%
  distinct(sample, sample_str, pred_log_mass_mg) %>%
  ggplot(aes(x = pred_log_mass_mg, fill = sample_str)) +
  geom_histogram(color = "black", position = "identity", alpha = 0.8) +
  theme_bw()

size2 %>%
  distinct(sample, sample_str) %>%
  group_by(sample_str) %>%
  summarise(n = n())

# Matrix for Sand PG Canopy -----------------------------------------------

#matrix for RDA with diet as columns, samples as rows
mat_pg <- size2 %>%
  filter(Island == "Sand" & 
           Habitat == "PG" & 
           Microhabitat == "Canopy") %>%
  dplyr::select(sample, Family) %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = Family,
              values_from = presence,
              values_fill = 0) %>% #fill missing values with 0's
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

# CCA PG Canopy -----------------------------------------------------------
#working off this tutorial
#http://dmcglinn.github.io/quant_methods/lessons/multivariate_models.html
# vegan requires that we write out each term if we are not going to 
# convert the factor to a dummy matrix 
# alternatively we could use a shorthand approach
cca_pg <-  cca(mat_pg ~ . , data=meta_pg)

#view the RDA loadings
cca_pg
#screeplot(cca_pg)
summary(cca_pg)
#how much variation explained in this RDA
RsquareAdj(cca_pg)
#ANOVA of whole model
anova(cca_pg, permutations=10000)
#ANOVA of model terms
anova(cca_pg, by='margin', permutations=10000)

#pretty ggplot plot 
#site metadata
sites <- meta_pg %>%
  rownames_to_column("site")
#get the site (sample) scores out and attach to site metadata
CCAscores <- scores(cca_pg, display = "sites") %>% 
  as.data.frame() %>% 
  rownames_to_column("site") %>%
  left_join(sites, by = "site")

#get the vectors out representing the loadings by species
CCAvect <- scores(cca_pg, display = c("cn")) %>% 
  as.data.frame() %>%
  rownames_to_column(var = "ID") %>%
  mutate(sample_str = str_sub(ID, -3))

#get vector out representing the loading by body size
CCAsizevect <- scores(cca_pg, display = "bp") %>% 
  as.data.frame() %>%
  rownames_to_column(var = "ID") %>%
  filter(ID == "pred_log_mass_mg")

cols <-  c("EUB" = "#EA7700", "HEV" = "#EEB00C","NEO" = "#89742F", "PHH" = "#158D8E")

CCA_plot <- ggplot() +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_point(data = CCAscores, aes(x = CCA1, y = CCA2, color = sample_str), size = 3) +
  geom_segment(data = CCAvect, 
               aes(x = 0, y = 0, xend = CCA1, yend = CCA2), 
               arrow = arrow(length = unit(0.2, "cm")),
               size = 1) +
  geom_segment(data = CCAsizevect,
               aes(x = 0, y =0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.2, "cm")),
               size = 1) +
  geom_text(data = CCAvect, aes(x = CCA1, y = CCA2, label = sample_str), 
            nudge_x = -0.15, size = 5) +
  geom_text(data = CCAsizevect, aes(x = CCA1, y = CCA2), label = "Predator mass",
            size = 5, nudge_y = 0.2, nudge_x = 0.1) +
  theme_bw() +
  scale_color_manual(values = cols) +
  labs(x = "CCA1 (39.8%)",
       y = "CCA2 (32.0%)") +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25))

# Euler plot --------------------------------------------------------------

#get variance by variable for Euler graph
#create environmental variables
species <- dummy(meta_pg$sample_str)
bs <- meta_pg$pred_log_mass_mg
#get variation by each
varpart(mat_pg, species, bs)
#create vector of these variations
var_pg <-  c("Predator species" = 16.2, 
             "Predator size" = 0.1,
             "Predator size&Predator species" = 3.4)

#make a Euler and plot it

fit3 <- euler(var_pg, shape = "ellipse")

euler <- plot(fit3,
              quantities = TRUE)

# Matrix minus PHH -----------------------------------------------

#matrix for RDA with diet as columns, samples as rows
mat_sm <- size2 %>%
  filter(Island == "Sand" & 
           Habitat == "PG" & 
           Microhabitat == "Canopy") %>%
  filter(sample_str != "PHH") %>%
  dplyr::select(sample, Family) %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = Family,
              values_from = presence,
              values_fill = 0) %>% #fill missing values with 0's
  column_to_rownames(var = "sample")

#metadata for RDA
meta_sm <- size2 %>%
  filter(Island == "Sand" & 
           Habitat == "PG" & 
           Microhabitat == "Canopy") %>%
  filter(sample_str != "PHH") %>%
  dplyr::select(sample, sample_str, pred_log_mass_mg) %>%
  distinct(sample, sample_str, pred_log_mass_mg) %>%
  column_to_rownames(var = "sample")

#check that rownames are the same:
all.equal(rownames(mat_sm), rownames(meta_sm))

# CCA PG Canopy w/o PHH -----------------------------------------------------------
#working off this tutorial
#http://dmcglinn.github.io/quant_methods/lessons/multivariate_models.html
# vegan requires that we write out each term if we are not going to 
# convert the factor to a dummy matrix 
# alternatively we could use a shorthand approach
cca_sm <-  cca(mat_sm ~ . , data=meta_sm)

#view the RDA loadings
cca_sm
#screeplot(cca_sm)
summary(cca_sm)
#how much variation explained in this RDA
RsquareAdj(cca_sm)
#ANOVA of whole model
anova(cca_sm, permutations=10000)
#ANOVA of model terms
anova(cca_sm, by='margin', permutations=10000)

#pretty ggplot plot 
#site metadata
sites_sm <- meta_sm %>%
  rownames_to_column("site")
#get the site (sample) scores out and attach to site metadata
CCAscores_sm <- scores(cca_sm, display = "sites") %>% 
  as.data.frame() %>% 
  rownames_to_column("site") %>%
  left_join(sites, by = "site")

#get the vectors out representing the loadings by species
CCAvect_sm <- scores(cca_sm, display = c("cn")) %>% 
  as.data.frame() %>%
  rownames_to_column(var = "ID") %>%
  mutate(sample_str = str_sub(ID, -3))

#get vector out representing the loading by body size
CCAsizevect_sm <- scores(cca_sm, display = "bp") %>% 
  as.data.frame() %>%
  rownames_to_column(var = "ID") %>%
  filter(ID == "pred_log_mass_mg")

cols_sm <-  c("EUB" = "#EA7700", "HEV" = "#EEB00C","NEO" = "#89742F")

CCA_plot_sm <- ggplot() +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_point(data = CCAscores_sm, aes(x = CCA1, y = CCA2, color = sample_str), size = 3) +
  geom_segment(data = CCAvect_sm, 
               aes(x = 0, y = 0, xend = CCA1, yend = CCA2), 
               arrow = arrow(length = unit(0.2, "cm")),
               size = 1) +
  geom_segment(data = CCAsizevect_sm,
               aes(x = 0, y =0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.2, "cm")),
               size = 1) +
  geom_text(data = CCAvect_sm, aes(x = CCA1, y = CCA2, label = sample_str), 
            nudge_x = -0.15, size = 5) +
  geom_text(data = CCAsizevect_sm, aes(x = CCA1, y = CCA2), label = "Predator mass",
            size = 5, nudge_y = 0.2, nudge_x = 0.1) +
  theme_bw() +
  scale_color_manual(values = cols_sm) +
  labs(x = "CCA1 (56.5%)",
       y = "CCA2 (27.2%)") +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25))
CCA_plot_sm

# Euler plot w/o PHH--------------------------------------------------------------

#get variance by variable for Euler graph
#create environmental variables
species <- dummy(meta_sm$sample_str)
bs <- meta_sm$pred_log_mass_mg
#get variation by each
varpart(mat_sm, species, bs)
#create vector of these variations
var_sm <-  c("Predator species" = 28.1, 
             "Predator size" = 0.7,
             "Predator size&Predator species" = 0)

#make a Euler and plot it

fit <- euler(var_sm, shape = "ellipse")

euler2 <- plot(fit,
              quantities = TRUE)
euler2
