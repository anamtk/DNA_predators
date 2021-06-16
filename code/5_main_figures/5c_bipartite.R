##########################
# bipartite vis of web use vs no web use -----
# Ana Miller-ter Kuile
# Feb 24, 2021
###########################

# make a bipartite vis of prey diet 

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
                       "2i_final_dataset", 
                       "pred_prey_sizes_DNAinteractions.csv"))

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

prey_mat <- prey %>%
  pivot_wider(names_from = sample_str,
              values_from = frequency,
              values_fill = 0) %>%
  column_to_rownames(var = "Family") %>%
  dplyr::select(-Order) %>%
  dplyr::select(HEV, PAN, PHH, EUB, SME, CEN,
                NEO, SCY, LRS)

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

plotweb(prey_mat,
        method = "normal", empty = TRUE, labsize = 1, ybig = 1,  y.width.low = 0.1, 
        y.width.high = 0.1, low.spacing = 0.02, high.spacing = 0.1,
        arrow="no",  col.interaction="grey80", col.high = "grey10", 
        col.low="grey10",  bor.col.interaction =NA, bor.col.high=NA, 
        bor.col.low=NA, high.lablength = NULL, low.lablength = NULL,
        sequence=NULL, low.abun=NULL, low.abun.col="green", 
        bor.low.abun.col ="black", high.abun=NULL, high.abun.col="red", 
        bor.high.abun.col="black", text.rot=90, text.high.col="black", 
        text.low.col="black", adj.high=NULL, adj.low=NULL, plot.axes = FALSE,
        low.y=1, high.y=1.5, add=FALSE, y.lim=NULL, x.lim=NULL, low.plot=TRUE, 
        high.plot=TRUE, high.xoff = 0.15, low.xoff = 0, high.lab.dis = NULL, 
        low.lab.dis = NULL, abuns.type="additional")

# Frequency sample size corrected -----------------------------------------
prednum <- data %>%
  distinct(sample, sample_str) %>%
  group_by(sample_str) %>%
  tally(name = "samp_size")

preycorr <- data %>%
  group_by(sample_str, Order, Family) %>%
  tally(name = "frequency") %>%
  left_join(prednum, by = "sample_str") %>%
  mutate(proportion = frequency/samp_size) %>%
  dplyr::select(-frequency, -samp_size) %>%
  arrange(Order) %>%
  mutate(Family = factor(Family, levels = c(Family)))

preycorr_mat <- preycorr %>%
  pivot_wider(names_from = sample_str,
              values_from = proportion,
              values_fill = 0) %>%
  column_to_rownames(var = "Family") %>%
  dplyr::select(-Order) %>%
  dplyr::select(HEV, PAN, PHH, EUB, SME, CEN,
                NEO, SCY, LRS)

plotweb(preycorr_mat,
        method = "normal", empty = TRUE, labsize = 1, ybig = 1,  y.width.low = 0.1, 
        y.width.high = 0.1, low.spacing = 0.02, high.spacing = 0.1,
        arrow="no",  col.interaction="grey80", col.high = "grey10", 
        col.low="grey10",  bor.col.interaction =NA, bor.col.high=NA, 
        bor.col.low=NA, high.lablength = NULL, low.lablength = NULL,
        sequence=NULL, low.abun=NULL, low.abun.col="green", 
        bor.low.abun.col ="black", high.abun=NULL, high.abun.col="red", 
        bor.high.abun.col="black", text.rot=90, text.high.col="black", 
        text.low.col="black", adj.high=NULL, adj.low=NULL, plot.axes = FALSE,
        low.y=1, high.y=1.5, add=FALSE, y.lim=NULL, x.lim=NULL, low.plot=TRUE, 
        high.plot=TRUE, high.xoff = 0.15, low.xoff = 0, high.lab.dis = NULL, 
        low.lab.dis = NULL, abuns.type="additional")

