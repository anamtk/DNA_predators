###########################
# rMangal web cleaning
# July 13, 2020
# Ana Miller-ter Kuile
###########################

#This script/function needs to take an interaction list from rmangal
#and extract only predation links from it, attached to the ID for
#each prey node

#then it needs to summarise those predatin link by predator species

###########################
# Load packages
library(here)
library(tidyverse)
library(ggplot2)
###########################

#need to break this into two for the lapply

#function 1. subset predation interactions

#Function 2. bind to the nodes file

#Function 3. compute per predator species links

pred_func <- function(interactions, nodes) {
  predation <- interactions %>%
    dplyr::select(node_from, node_to, type, method) %>%
    filter(type == "predation") %>%
    left_join(nodes, by = c("node_to" = "node_id"))
  return(predation)
}

per_pred_func <- function(predation) {
  per_pred <- predation %>%
    group_by(node_from, taxon_Family) %>%
    tally(name = "links") %>%
    group_by(node_from) %>%
    summarise(links = n())
  return(per_pred)
}


files <- list.files(path="path/to/dir", pattern="*.txt", full.names=TRUE, recursive=FALSE)
lapply(files, function(x) {
  t <- read.table(x, header=TRUE) # load file
  # apply function
  out <- function(t)
    # write to file
    write.table(out, "path/to/output", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
})

hines_int <- read.csv(here("Published_webs", "final_cut", "hines_2019", "interactions_hines.csv"))

hines_node <- read.csv(here("Published_webs", "final_cut", "hines_2019", "nodes_hines.csv"))

hines_pred <- pred_func(hines_int, hines_node)

hines_per_pred <- per_pred_func(hines_pred)
###########################
# Export data
###########################
write.csv(predation, here("data", "outputs", 
                          "6_pub_webs", "hines_predation.csv"))

write.csv(per_pred, here("data", "outputs", 
                         "6_pub_webs", "hines_per_pred.csv"))


