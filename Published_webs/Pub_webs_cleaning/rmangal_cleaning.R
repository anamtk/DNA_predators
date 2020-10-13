###########################
# rMangal web cleaning
# July 13, 2020
# Ana Miller-ter Kuile
###########################

#This script/function needs to take an interaction list from rmangal
#and extract only predation links from it, attached to the ID for
#each prey node

#then it needs to summarise those predation link by predator species

###########################
# Load packages
library(here)
library(tidyverse)
library(ggplot2)
###########################

###########################
# All predation links
###########################

    # path to folder that holds multiple .csv files
folder <- here("Published_webs", "rmangal")
nodes <- list.files(path=folder, pattern="nodes*", recursive = T, full.names = T) # create list of all .csv files in folder
interactions <- list.files(path=folder, pattern="interactions*", recursive = T, full.names = T)
#webs <- list(ints, nodes)

for(i in 1:length(nodes)){
  ## Get only the relevant .csvs in each pair
  node <- read.csv(nodes[i])
  interaction <- read.csv(interactions[i])
  
  ## Manipulate and merge
  node <- node %>%
    add_tally(name = "species_richness")
  
  predation <- interaction %>%
     dplyr::select(node_from, node_to, type, method) %>%
     filter(type == "predation") %>%
     left_join(node, by = c("node_from" = "node_id"))
  
  ## Get web ID
  underscores <- str_locate_all(nodes[i], "s_")[[1]]
  predation$web <- as.character(str_sub(nodes[i], start = underscores[nrow(underscores),2]+1, end = -5))

  ## Add to output tbl
  if(i == 1){
    outTbl <- predation
  }else{
    outTbl <- bind_rows(outTbl, predation)
  }
}

outTbl <- outTbl %>%
  dplyr::select(node_from, node_to, type, method, original_name,
                taxon_Kingdom, taxon_Class, taxon_Order, taxon_Family, 
                taxonomy.name, species_richness, web) %>%
  rename("consumer" = "node_to",
         "resource" = "node_from")

###########################
# Export data
###########################
write.csv(outTbl, here("data", "outputs", 
                          "6_pub_webs", "rmangal_webs_predation.csv"))


