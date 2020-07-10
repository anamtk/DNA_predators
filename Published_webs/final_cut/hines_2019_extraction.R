library(rmangal)
library(here)
library(tidyverse)


germany <- search_datasets(list(id = 178)) %>%
  get_collection()

nodes <- (germany$nodes)
interactions <- (germany$interactions)
