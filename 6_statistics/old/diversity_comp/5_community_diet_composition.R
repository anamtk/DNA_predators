###########################
#Diet Composition by community
#Ana Miller-ter Kuile
#June 9, 2020
###########################

#This code will look at diet composition within a community of predators
#from the same sampling area - exploring how species and body size
#determines diet in these environments. 

#How does species identity or body size influence diet composition
#within a predator community?
##########################
#Load packages####
###########################
library(here)
library(tidyverse)
library(ggplot2)
library(vegan)

#########################
#Load data####
###########################
#data from the taxonomic sort, which is already concatenated
#by unique ID within a predator individual
dna <- read.csv(here("data", "outputs", "5_rarefied_taxonomic_sort", "prey_DNA.csv"))

#sample metadata for sample sizes
meta <- read.csv(here("data", "Sample_metadata.csv"))

#########################
#assess possible communities####
###########################
dna$sample <- str_sub(dna$sample, start = 1, end =-2) #remove "a" at end

test <- dna %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  filter(Year == "2017") %>%
  filter(Habitat != "Lab area") %>%
  distinct(sample, pred_ID, Method, Island, Habitat, 
           Microhabitat, Date.Collected, No..Individuals) %>%
  dplyr::select(sample, pred_ID, Method, Island, Habitat, 
                Microhabitat, Date.Collected, No..Individuals) %>%
  group_by(Island, Habitat, Microhabitat) %>%
  summarise(final_indiv = sum(No..Individuals), final_samples = n(), final_species = n_distinct(pred_ID))

test2 <- meta %>%
  filter(Year == "2017") %>%
  filter(Habitat != "Lab area") %>%
  distinct(Island, Habitat, Microhabitat, ID, Extraction.ID, No..Individuals) %>%
  group_by(Island, Habitat, Microhabitat) %>%
  summarise(tot_Indivs = sum(No..Individuals), tot_samples = n(), tot_species = n_distinct(ID))

comm_samples <- test %>%
  left_join(test2, by =c("Island", "Habitat", "Microhabitat"))


test3 <- dna %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  distinct(sample, pred_ID, Method, Island, Habitat, 
           Microhabitat, Date.Collected, No..Individuals) %>%
  dplyr::select(sample, pred_ID, Method, Island, Habitat, 
                Microhabitat, Date.Collected, No..Individuals) %>%
  group_by(pred_ID) %>%
  summarise(Indivs = sum(No..Individuals), samples = n())
  
#########################
#PF Cooper canopy - experimental community####
###########################

san_pg_can <- dna %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  filter(Island == "Sand" & Habitat == "PG" & Microhabitat == "Canopy") %>%
  group_by(sample, unique_ID, pred_ID, Class_prey, Order_prey, Family_prey,
           Genus_prey, Species_prey, ID_level, reads, Method, Island, 
           Habitat, Microhabitat, Year, No..Individuals) %>%
  summarise(mean_mm_length = mean(Length_mm), 
            total = n(), 
            sd = sd(Length_mm), 
            se=sd/sqrt(total))

san_pg_meta <- san_pg_can %>%
  ungroup() %>%
  dplyr::select(sample, pred_ID, mean_mm_length) %>%
  distinct(sample, pred_ID, mean_mm_length)

san_pg_dat <- san_pg_can %>%
  ungroup() %>%
  dplyr::select(sample, unique_ID, reads) %>%
  group_by(sample, unique_ID) %>%
  summarise(reads = sum(reads)) %>%
  pivot_wider(names_from = unique_ID,
              values_from = reads) %>%
  column_to_rownames(var = "sample")


san_pg_dat <- san_pg_dat[, which(colSums(san_pg_dat) != 0)] #108 to 45
#I’m using the metaMDS function from the vegan package
nmds1 <- metaMDS(san_pg_dat, distance = "jaccard", binary = TRUE, k=2)
#then look at a stress plot. This is to evaluate whether the NMDS is actually a good representation of the structure of your data.
#I don’t know a ton about evaluating these, but it should be roughly a linear stair step
stressplot(nmds1)

#making a data frame of the points from your NMDS
nmds1_df<-data.frame(MDS1=nmds1$points[,1], MDS2=nmds1$points[,2])
#combining with metadata for plotting
nmds1_df_meta<-cbind(san_pg_meta, nmds1_df)

#I’m gong to plot my points in different colors based on the “Site” that they came from, so I need to make sure that my “Site” is a factor, not numeric.
nmds1_df_meta$pred_ID<-as.factor(nmds1_df_meta$pred_ID)

#using ggplot to plot it
#only points
plotpoints<-ggplot(nmds1_df_meta, aes(x=MDS1,y=MDS2, color=pred_ID, size = mean_mm_length))+
  geom_point()+
  coord_fixed()+
  theme_bw()
plotpoints
#ellipse
ellipse1<-ggplot(nmds1_df_meta, aes(x=MDS1,y=MDS2, color=pred_ID, size = mean_mm_length))+
  geom_point()+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()
ellipse1

ellipse2<-ggplot(nmds1_df_meta, aes(x=MDS1,y=MDS2, color=mean_mm_length))+
  geom_point(size = 6)+
  #stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()
ellipse2


