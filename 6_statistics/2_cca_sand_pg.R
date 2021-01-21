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
  mutate(pred_mass_mg = 10^(pred_log_mass_mg)) 

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

# Matrix for SAND PG -----------------------------------------------
sand_pg <- size2 %>%
  filter(Island == "Sand" & Habitat == "PG") %>%
  filter(Microhabitat != "Soil")

sand_pg %>%
  distinct(sample) %>%
  tally()

# Sand PG CCA -----------------------------------------------------------------

#matrix with diet as columns, samples as rows
mat_pg <- sand_pg %>%
  dplyr::select(sample, Family) %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = Family,
              values_from = presence,
              values_fill = 0) %>% #fill missing values with 0's
  column_to_rownames(var = "sample")

#metadata
meta_pg <- sand_pg %>%
  dplyr::select(sample, sample_str, pred_log_mass_mg) %>%
  distinct(sample, sample_str, pred_log_mass_mg) %>%
  column_to_rownames(var = "sample")

#check that rownames are the same:
all.equal(rownames(mat_pg), rownames(meta_pg))

# CCA -----------------------------------------------------------
#working off this tutorial
#http://dmcglinn.github.io/quant_methods/lessons/multivariate_models.html
# vegan requires that we write out each term if we are not going to 
# convert the factor to a dummy matrix 
# alternatively we could use a shorthand approach
cca_pg <-  cca(mat_pg ~ . , data=meta_pg)

#view the CCA loadings
cca_pg
#screeplot(cca_pg)
summary(cca_pg)
#how much variation explained in this RDA
RsquareAdj(cca_pg)
#ANOVA of whole model
anova(cca_pg, permutations=10000)
#ANOVA of model terms
anova(cca_pg, by='margin', permutations=10000)

cca_pg_simple <- update(cca_pg, . ~ . - pred_log_mass_mg)
anova(cca_pg_simple, cca_pg)

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
  filter(str_detect(ID, "sample_str")) %>%
  mutate(sample_str = str_sub(ID, -3)) %>%
  filter(sample_str %in% c("EUB", "HEV", "NEO", "PHH")) %>%
  mutate(pred_sp = ifelse(sample_str == "EUB", "E. annulipes",
                          ifelse(sample_str == "HEV", "H. venatoria",
                                 ifelse(sample_str == "NEO", "N. theisi",
                                        "P. holdhausi"))))

#get vector out representing the loading by body size
CCAsizevect <- scores(cca_pg, display = "bp") %>% 
  as.data.frame() %>%
  rownames_to_column(var = "ID") %>%
  filter(ID == "pred_log_mass_mg")

pal_kelp <- cal_palette("kelp1", n = 9, type = "continuous")

pred_labels <- c("CEN" = "Geophilomorpha sp.", "EUB" = "E. annulipes", 
                 "HEV" = "H. venatoria", "LRS" = "Oonopidae sp.", 
                 "NEO" = "N. theisi", "PAN" = "P. flavescens",
                 "PHH" = "P. holdhausi", "SCY" = "S. longipes",
                 "SME" = "S. pallidus")

#LRS: "#C70000" 
# SCY: "#EA7700" 
#NEO: "#EEB00C" 
#CEN: "#C68A2C" 
#SME: "#89742F" 
#EUB: "#496C3C" 
#PHH: "#158D8E"
#PAN: "#067D8D" 
#HEV: "#114C54"

ggplot() +
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
  scale_color_manual(values = c("#496C3C", "#114C54", "#EEB00C", "#158D8E"),
                     labels = pred_labels) +
  labs(x = "CCA1 (5.7%)",
       y = "CCA2 (5.0%)",
       colour = "Predator species") +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25))

sand_pg %>%
  distinct(sample, pred_mass_mg, sample_str) %>%
  ggplot(aes(x = pred_mass_mg, fill = sample_str)) +
  geom_histogram(position = "identity", 
                 color = "black", alpha = 0.6) +
  scale_fill_manual(values = c("#496C3C", "#114C54", "#EEB00C", "#158D8E"),
                     labels = pred_labels) +
  scale_x_log10() +
  labs(x = "Predator mass (mg)",
       y = "Count", 
       fill = "Predator species") +
  theme_bw()
