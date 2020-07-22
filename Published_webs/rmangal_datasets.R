##########################
#rmangal data download
#July 15, 2020
#Ana Miller-ter Kuile
##########################

#This code downloads datasets from rmangal that include invertebrate
#terrestrial predation AND have nodes resolved at family level or
#better. I'll be using these to ask how well the molecular data do
#compared to published invertebrate interactions and interaction
#networks

#I already downloadaed the hines 2019 dataset earlier in its own
#code, so this is everything but that one. 

##########################
#Load packages
##########################
library(here)
library(tidyverse)
library(rmangal)

##########################
#Download and write each dataset
##########################

##########################
#Cornaby
##########################
cornaby <- search_datasets(list(id = 91)) %>%
  get_collection()
#2nd one has predation:
corn_node <- cornaby[[2]]$nodes
corn_int <- cornaby[[2]]$interactions

write.csv(corn_node, here("Published_webs", "final_cut", "Cornaby_1974", "nodes_cornaby.csv"))

write.csv(corn_int, here("Published_webs", "final_cut", "Cornaby_1974", "interactions_cornaby.csv"))

##########################
#Mohr
##########################
mohr <- search_datasets(list(id =107)) %>%
  get_collection()
mohr_node <- mohr$nodes
mohr_inter <- mohr$interactions

write.csv(mohr_node, here("Published_webs", "final_cut", "Mohr_1943", "nodes_mohr.csv"))

write.csv(mohr_inter, here("Published_webs", "final_cut", "Mohr_1943", "interactions_mohr.csv"))

##########################
#Robinson
##########################
robinson <- search_datasets(list(id =114)) %>%
  get_collection()
rob_int <- robinson$interactions
rob_node <- robinson$nodes

write.csv(rob_node, here("Published_webs", "final_cut", "Robinson_1953", "nodes_robinson.csv"))

write.csv(rob_int, here("Published_webs", "final_cut", "Robinson_1953", "interactions_robinson.csv"))

##########################
#Savely
##########################

savely <- search_datasets(list(id =115)) %>%
  get_collection()
save1_node <- savely[[1]]$nodes
save2_node <- savely[[2]]$nodes
save3_node <- savely[[3]]$nodes
save4_node <- savely[[4]]$nodes

save1_int <- savely[[1]]$interactions
save2_int <- savely[[2]]$interactions
save3_int <- savely[[3]]$interactions
save4_int <- savely[[4]]$interactions

write.csv(save1_node, here("Published_webs", "final_cut", "Savely_1939", "nodes_save1.csv"))
write.csv(save2_node, here("Published_webs", "final_cut", "Savely_1939", "nodes_save2.csv"))
write.csv(save3_node, here("Published_webs", "final_cut", "Savely_1939", "nodes_save3.csv"))
write.csv(save4_node, here("Published_webs", "final_cut", "Savely_1939", "nodes_save4.csv"))

write.csv(save1_int, here("Published_webs", "final_cut", "Savely_1939", "interactions_save1.csv"))
write.csv(save2_int, here("Published_webs", "final_cut", "Savely_1939", "interactions_save2.csv"))
write.csv(save3_int, here("Published_webs", "final_cut", "Savely_1939", "interactions_save3.csv"))
write.csv(save4_int, here("Published_webs", "final_cut", "Savely_1939", "interactions_save4.csv"))

##########################
#Schoenly
##########################

schoenly <- search_datasets(list(id =117)) %>%
  get_collection()
scho_node <- schoenly$nodes
scho_int <- schoenly$interactions

write.csv(scho_node, here("Published_webs", "final_cut", "Schoenly_1983", "nodes_schoenly.csv"))
write.csv(scho_int, here("Published_webs", "final_cut", "Schoenly_1983", "interactions_schoenly.csv"))

##########################
#Valeila
##########################

valeila <- search_datasets(list(id =124)) %>%
  get_collection()
val1_node <- valeila[[1]]$nodes
val2_node <- valeila[[2]]$nodes

val1_int <- valeila[[1]]$interactions
val2_int <- valeila[[2]]$interactions

write.csv(val1_node, here("Published_webs", "final_cut", "Valeila_1974", "nodes_val1.csv"))
write.csv(val2_node, here("Published_webs", "final_cut", "Valeila_1974", "nodes_val2.csv"))
write.csv(val1_int, here("Published_webs", "final_cut", "Valeila_1974", "interactions_val1.csv"))
write.csv(val2_int, here("Published_webs", "final_cut", "Valeila_1974", "interactions_val2.csv"))

##########################
#Whittaker
##########################

whitt <- search_datasets(list(id =126)) %>%
  get_collection()
whit_node <- whitt$nodes
whit_int <- whitt$interactions

write.csv(whit_node, here("Published_webs", "final_cut", "Whittaker_1984", "nodes_whittaker.csv"))
write.csv(whit_int, here("Published_webs", "final_cut", "Whittaker_1984", "interactions_whittaker.csv"))

##########################
#McKinnerney
##########################
mckin <- search_datasets(list(id =151)) %>%
  get_collection()

mckin2_int <- mckin[[2]]$interactions #2,4,5, 6
mckin2_nodes <- mckin[[2]]$nodes

mckin4_int <- mckin[[4]]$interactions #2,4,5, 6
mckin4_nodes <- mckin[[4]]$nodes

mckin5_int <- mckin[[5]]$interactions #2,4,5, 6
mckin5_nodes <- mckin[[5]]$nodes

mckin6_int <- mckin[[6]]$interactions #2,4,5, 6
mckin6_nodes <- mckin[[6]]$nodes

write.csv(mckin2_nodes, here("Published_webs", "final_cut", "Mckinnerney_1978", "nodes_mckin2.csv"))
write.csv(mckin2_int, here("Published_webs", "final_cut", "Mckinnerney_1978", "interactions_mckin2.csv"))
write.csv(mckin4_nodes, here("Published_webs", "final_cut", "Mckinnerney_1978", "nodes_mckin4.csv"))
write.csv(mckin4_int, here("Published_webs", "final_cut", "Mckinnerney_1978", "interactions_mckin4.csv"))
write.csv(mckin5_nodes, here("Published_webs", "final_cut", "Mckinnerney_1978", "nodes_mckin5.csv"))
write.csv(mckin5_int, here("Published_webs", "final_cut", "Mckinnerney_1978", "interactions_mckin5.csv"))
write.csv(mckin6_nodes, here("Published_webs", "final_cut", "Mckinnerney_1978", "nodes_mckin6.csv"))
write.csv(mckin6_int, here("Published_webs", "final_cut", "Mckinnerney_1978", "interactions_mckin6.csv"))

##########################
#McKinnerney
##########################

choen <- search_datasets(list(id =152)) %>%
  get_collection()

choen1_node <- choen[[1]]$nodes
choen1_int <- choen[[1]]$interactions

choen2_node <- choen[[2]]$nodes
choen2_int <- choen[[2]]$interactions

write.csv(choen1_node, here("Published_webs", "final_cut", "Choenly_reid_1983", "nodes_choen1.csv"))
write.csv(choen1_int, here("Published_webs", "final_cut", "Choenly_reid_1983", "interactions_choen1.csv"))

write.csv(choen2_node, here("Published_webs", "final_cut", "Choenly_reid_1983", "nodes_choen2.csv"))
write.csv(choen2_int, here("Published_webs", "final_cut", "Choenly_reid_1983", "interactions_choen2.csv"))

##########################
#Richards
##########################
rich <- search_datasets(list(id =154)) %>%
  get_collection()
rich_node <- rich$nodes
rich_int <- rich$interactions

write.csv(rich_node, here("Published_webs", "final_cut", "Richards_1926", "nodes_rich.csv"))
write.csv(rich_int, here("Published_webs", "final_cut", "Richards_1926", "interactions_rich.csv"))


