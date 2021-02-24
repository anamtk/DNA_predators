##########################
# bipartite vis of web use vs no web use -----
# Ana Miller-ter Kuile
# Feb 24, 2021
###########################

# make a bipartite vis of prey diet by web building

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 
                  "bipartite", "calecopal")

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

#select variables to remove
size <- data %>%
  dplyr::select(-X, -reads) %>%
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
size <- size %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  unite(environment, c("Habitat", "Microhabitat"), sep = "_", remove = F) 

# Data summary -------------------------------------------------------

# Matrix for vis -----------------------------------------------
prey <- data %>%
  group_by(sample_str, Order, Family) %>%
  tally(name = "frequency")

prey <- prey %>%
  arrange(Order) %>%
  mutate(Family = factor(Family, levels = c(Family)))

prey_IDs <- prey %>%
  ungroup() %>%
  dplyr::select(Order, Family) %>%
  distinct(Order, Family)

prey <- prey %>%
  mutate(predator = "all")

prey_mat <- prey %>%
  pivot_wider(names_from = sample_str,
              values_from = frequency,
              values_fill = 0) %>%
  column_to_rownames(var = "Family") %>%
  dplyr::select(-Order)

plotweb(prey_mat,
        method = "cca", empty = TRUE, labsize = 1, ybig = 1,  y.width.low = 0.1, 
        y.width.high = 0.1, low.spacing = 0.04, high.spacing = NULL,
        arrow="no",  col.interaction="grey80", col.high = "grey10", 
        col.low="grey10",  bor.col.interaction =NA, bor.col.high=NA, 
        bor.col.low=NA, high.lablength = NULL, low.lablength = NULL,
        sequence=NULL, low.abun=NULL, low.abun.col="green", 
        bor.low.abun.col ="black", high.abun=NULL, high.abun.col="red", 
        bor.high.abun.col="black", text.rot=90, text.high.col="black", 
        text.low.col="black", adj.high=NULL, adj.low=NULL, plot.axes = FALSE,
        low.y=1, high.y=1.5, add=FALSE, y.lim=NULL, x.lim=NULL, low.plot=TRUE, 
        high.plot=TRUE, high.xoff = 0, low.xoff = 0, high.lab.dis = NULL, 
        low.lab.dis = NULL, abuns.type="additional")

plotweb(prey_mat,
        method = "normal", empty = TRUE, labsize = 1, ybig = 1,  y.width.low = 0.1, 
        y.width.high = 0.1, low.spacing = 0.04, high.spacing = NULL,
        arrow="no",  col.interaction="grey80", col.high = "grey10", 
        col.low="grey10",  bor.col.interaction =NA, bor.col.high=NA, 
        bor.col.low=NA, high.lablength = NULL, low.lablength = NULL,
        sequence=NULL, low.abun=NULL, low.abun.col="green", 
        bor.low.abun.col ="black", high.abun=NULL, high.abun.col="red", 
        bor.high.abun.col="black", text.rot=90, text.high.col="black", 
        text.low.col="black", adj.high=NULL, adj.low=NULL, plot.axes = FALSE,
        low.y=1, high.y=1.5, add=FALSE, y.lim=NULL, x.lim=NULL, low.plot=TRUE, 
        high.plot=TRUE, high.xoff = 0, low.xoff = 0, high.lab.dis = NULL, 
        low.lab.dis = NULL, abuns.type="additional")

cal_palette("kelp2", n = 12, type = "continuous")
