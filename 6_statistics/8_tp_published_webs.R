#############################
#Trophic positions of published vs. HTS links
#August 18, 2020
#Ana Miller-ter Kuile
#############################

#This script examines whether the trophic positions of diet 
#items assigned via different published methods differ
#from those from the HTS data
#This is important because of the importance of intraguild
#predation in food web ecology

#############################
#Load packages
library(here)
library(tidyverse)
library(ggplot2)
library(glmmTMB)
library(DHARMa)
library(MuMIn)
library(effects)
library(emmeans)
library(ggeffects)
library(vegan)
#############################

#note, it may be important to break this analysis up
#by first: more specific trophic groups and then
#by courser scale - e.g. are they eating plants, vs. 
#other animals, or both.

#############################
#Load data
#############################

links <- read.csv(here("data", "outputs", 
                       "7_all_webs", 
                       "all_interactions_and_tp.csv"))

#############################
#Summarize by specific TP
#############################

links_tp <- links %>%
  unite(web_consumer, consumer, web, remove = F) %>%
  group_by(web_consumer, web, family_richness, coll_method,
           pub_year, tp) %>%
  tally(name = "link_number") 

meta_tp <- links_tp %>%
  filter(!is.na(tp)) %>%
  dplyr::select(web_consumer, web, family_richness, coll_method,
                pub_year) %>%
  distinct()

tp_matrix <- links_tp %>%
  ungroup() %>%
  dplyr::select(web_consumer, tp, link_number) %>%
  filter(!is.na(tp)) %>%
  pivot_wider(names_from = "tp",
              values_from = "link_number") 

tp_matrix[is.na(tp_matrix)] <- 0

rownames(tp_matrix) <- tp_matrix$web_consumer
tp_matrix <- tp_matrix %>%
  dplyr::select(-web_consumer)

tp_matrix <- tp_matrix[rowSums(tp_matrix) > 0,] #1291


#############################
#NMDS by specific TP
#############################

#I’m using the metaMDS function from the vegan package
nmds1 <- metaMDS(tp_matrix, distance = "bray", k=2)
#then look at a stress plot. This is to evaluate whether the NMDS is actually a good representation of the structure of your data.
#I don’t know a ton about evaluating these, but it should be roughly a linear stair step
stressplot(nmds1)

#making a data frame of the points from your NMDS
nmds1_df<-data.frame(MDS1=nmds1$points[,1], MDS2=nmds1$points[,2])


#combining with metadata for plotting
nmds1_df_meta<-nmds1_df %>%
  bind_cols(meta_tp)

#I’m gong to plot my points in different colors based on the “Site” that they came from, so I need to make sure that my “Site” is a factor, not numeric.
nmds1_df_meta$web<-as.factor(nmds1_df_meta$web)
nmds1_df_meta$coll_method <- as.factor(nmds1_df_meta$coll_method)

nmds1_df_meta_1 <- nmds1_df_meta %>%
  filter(MDS1 > -34)
#using ggplot to plot it
#only points
ggplot(nmds1_df_meta_1, aes(x=MDS1,y=MDS2, color=web))+
  geom_point(mapping = aes(x=MDS1, y=MDS2),
             size = 6)+
  coord_fixed()+
  theme_bw()

#ellipse
ggplot(nmds1_df_meta_1, aes(x=MDS1,y=MDS2, color=coll_method))+
  geom_point(mapping = aes(x=MDS1, y=MDS2),
             size = 6)+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()


#############################
#Summarize by broad TP
#############################

links_btp <- links %>%
  unite(web_consumer, consumer, web, remove = F) %>%
  group_by(web_consumer, web, family_richness, coll_method,
           pub_year, broad_tp) %>%
  tally(name = "link_number") 

meta_btp <- links_btp %>%
  filter(broad_tp != "") %>%
  dplyr::select(web_consumer, web, family_richness, coll_method,
                pub_year) %>%
  distinct()

btp_matrix <- links_btp %>%
  ungroup() %>%
  dplyr::select(web_consumer, broad_tp, link_number) %>%
  filter(broad_tp != "") %>%
  pivot_wider(names_from = "broad_tp",
              values_from = "link_number") 

btp_matrix[is.na(btp_matrix)] <- 0

rownames(btp_matrix) <- btp_matrix$web_consumer
btp_matrix <- btp_matrix %>%
  dplyr::select(-web_consumer)

btp_matrix <- btp_matrix[rowSums(btp_matrix) > 0,] #1291

#############################
#NMDS by broad TP
#############################

#I’m using the metaMDS function from the vegan package
nmds1 <- metaMDS(btp_matrix, distance = "bray", k=2)
#then look at a stress plot. This is to evaluate whether the NMDS is actually a good representation of the structure of your data.
#I don’t know a ton about evaluating these, but it should be roughly a linear stair step
stressplot(nmds1)

#making a data frame of the points from your NMDS
nmds1_df<-data.frame(MDS1=nmds1$points[,1], MDS2=nmds1$points[,2])


#combining with metadata for plotting
nmds1_df_meta<-nmds1_df %>%
  bind_cols(meta_btp)

#I’m gong to plot my points in different colors based on the “Site” that they came from, so I need to make sure that my “Site” is a factor, not numeric.
nmds1_df_meta$web<-as.factor(nmds1_df_meta$web)
nmds1_df_meta$coll_method <- as.factor(nmds1_df_meta$coll_method)

#using ggplot to plot it
#only points
ggplot(nmds1_df_meta_1, aes(x=MDS1,y=MDS2, color=web))+
  geom_point(mapping = aes(x=MDS1, y=MDS2),
             size = 6)+
  coord_fixed()+
  theme_bw()

#ellipse
ggplot(nmds1_df_meta_1, aes(x=MDS1,y=MDS2, color=coll_method))+
  geom_point(mapping = aes(x=MDS1, y=MDS2),
             size = 2)+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()

links_btp <- links_btp %>%
  filter(broad_tp != "")

links_sp <- links %>%
  unite(web_consumer, consumer, web, remove = F) %>%
  group_by(web_consumer, web, family_richness, coll_method,
           pub_year) %>%
  tally(name = "total_link_number") %>%
  ungroup() %>%
  dplyr::select(web_consumer, total_link_number)

links_btp <- links_btp %>%
  left_join(links_sp, by = "web_consumer")

ggplot(links_btp, aes(x = coll_method, y = link_number/total_link_number, fill = broad_tp)) +
  geom_boxplot() + theme_bw()
