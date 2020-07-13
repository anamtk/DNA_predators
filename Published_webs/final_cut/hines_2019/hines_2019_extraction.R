library(rmangal)
library(here)
library(tidyverse)


germany <- search_datasets(list(id = 178)) %>%
  get_collection()

nodes <- (germany$nodes)
interactions <- (germany$interactions)

write.csv(nodes, here("Published_webs", "final_cut", "hines_2019", "nodes_hines.csv"))

write.csv(interactions, here("Published_webs", "final_cut", "hines_2019", "interactions_hines.csv"))
