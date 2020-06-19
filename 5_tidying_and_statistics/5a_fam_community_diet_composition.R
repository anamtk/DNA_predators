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
dna <- read.csv(here("data", "outputs", "5_rarefied_taxonomic_sort", "fam_prey_DNA.csv"))

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
#PG sand canopy - experimental community####
###########################

san_pg_can <- dna %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  filter(Island == "Sand" & Habitat == "PG" & Microhabitat == "Canopy") %>%
  group_by(sample, Family, pred_ID, reads, Method, Island, 
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
  dplyr::select(sample, Family, reads) %>%
  group_by(sample, Family) %>%
  summarise(reads = sum(reads)) %>%
  filter(!Family %in% c("Muridae", "Hominidae")) %>%
  pivot_wider(names_from = Family,
              values_from = reads) %>%
  column_to_rownames(var = "sample") %>%
  replace(., is.na(.), 0)

zeros <- san_pg_dat[which(rowSums(san_pg_dat) == 0),]

zeros$sample <- rownames(zeros)


san_pg_dat <- san_pg_dat[, which(colSums(san_pg_dat) != 0)] #67 to 28
san_pg_dat <- san_pg_dat[which(rowSums(san_pg_dat) != 0),] #72 to 66
#I’m using the metaMDS function from the vegan package
nmds1 <- metaMDS(san_pg_dat, distance = "jaccard", binary = TRUE, k=2, trymax = 1000)
#then look at a stress plot. This is to evaluate whether the NMDS is actually a good representation of the structure of your data.
#I don’t know a ton about evaluating these, but it should be roughly a linear stair step
stressplot(nmds1)

san_pg_meta <- san_pg_meta %>%
  anti_join(zeros, by = "sample")
  

#making a data frame of the points from your NMDS
nmds1_df<-data.frame(MDS1=nmds1$points[,1], MDS2=nmds1$points[,2])
#combining with metadata for plotting
nmds1_df<-cbind(san_pg_meta, nmds1_df)

#I’m gong to plot my points in different colors based on the “Site” that they came from, so I need to make sure that my “Site” is a factor, not numeric.
nmds1_df$pred_ID<-as.factor(nmds1_df$pred_ID)

#using ggplot to plot it
#only points
plotpoints<-ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=pred_ID, size = mean_mm_length))+
  geom_point()+
  coord_fixed()+
  theme_bw()
plotpoints
#ellipse
ellipse1<-ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=pred_ID, size = mean_mm_length))+
  geom_point()+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()
ellipse1

ellipse2<-ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=mean_mm_length))+
  geom_point(size = 6)+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()
ellipse2

nmds1_df$quantiles <- ifelse(nmds1_df$mean_mm_length > 0 & nmds1_df$mean_mm_length <= 1.7, 0,
                                  ifelse(nmds1_df$mean_mm_length > 1.7 & nmds1_df$mean_mm_length <= 7.5625, 1,
                                         ifelse(nmds1_df$mean_mm_length > 7.5625 & nmds1_df$mean_mm_length <= 11.2, 2,
                                                ifelse(nmds1_df$mean_mm_length > 11.2 & nmds1_df$mean_mm_length <= 13.775, 3, 4
                                                       ))))


nmds1_df$quantiles <- factor(nmds1_df$quantiles, levels = c("0", "1", "2", "3", "4"))
ellipse3<-ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=quantiles))+
  geom_point(size = 6)+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()
ellipse3


#########################
#PF canopy ALL - experimental community####
###########################

pf_can <- dna %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  filter(Habitat == "PF" & Microhabitat == "Canopy") %>%
  group_by(sample, Family, pred_ID, reads, Method, Island, 
           Habitat, Microhabitat, Year, No..Individuals) %>%
  summarise(mean_mm_length = mean(Length_mm), 
            total = n(), 
            sd = sd(Length_mm), 
            se=sd/sqrt(total))

pf_meta <- pf_can %>%
  ungroup() %>%
  dplyr::select(sample, pred_ID, mean_mm_length) %>%
  distinct(sample, pred_ID, mean_mm_length)

pf_dat <- pf_can %>%
  ungroup() %>%
  dplyr::select(sample, Family, reads) %>%
  group_by(sample, Family) %>%
  summarise(reads = sum(reads)) %>%
  filter(!Family %in% c("Muridae", "Hominidae")) %>%
  pivot_wider(names_from = Family,
              values_from = reads) %>%
  column_to_rownames(var = "sample") %>%
  replace(., is.na(.), 0)



pf_dat <- pf_dat[, which(colSums(pf_dat) != 0)] #67 to 24
pf_dat <- pf_dat[which(rowSums(pf_dat) != 0),] #none
#I’m using the metaMDS function from the vegan package
nmds2 <- metaMDS(pf_dat, distance = "jaccard", binary = TRUE, k=2, trymax = 1000)
#then look at a stress plot. This is to evaluate whether the NMDS is actually a good representation of the structure of your data.
#I don’t know a ton about evaluating these, but it should be roughly a linear stair step
stressplot(nmds2)
#making a data frame of the points from your NMDS
nmds2_df<-data.frame(MDS1=nmds2$points[,1], MDS2=nmds2$points[,2])
#combining with metadata for plotting
nmds2_df<-cbind(pf_meta, nmds2_df)

#I’m gong to plot my points in different colors based on the “Site” that they came from, so I need to make sure that my “Site” is a factor, not numeric.
nmds2_df$pred_ID<-as.factor(nmds2_df$pred_ID)

#using ggplot to plot it
#only points
plotpoints<-ggplot(nmds2_df, aes(x=MDS1,y=MDS2, color=pred_ID, size = mean_mm_length))+
  geom_point()+
  coord_fixed()+
  theme_bw()
plotpoints
#ellipse
ellipse1<-ggplot(nmds2_df, aes(x=MDS1,y=MDS2, color=pred_ID))+
  geom_point(size = 6)+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()
ellipse1

ellipse2<-ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=mean_mm_length))+
  geom_point(size = 6)+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()
ellipse2

nmds1_df$quantiles <- ifelse(nmds1_df$mean_mm_length > 0 & nmds1_df$mean_mm_length <= 1.7, 0,
                             ifelse(nmds1_df$mean_mm_length > 1.7 & nmds1_df$mean_mm_length <= 7.5625, 1,
                                    ifelse(nmds1_df$mean_mm_length > 7.5625 & nmds1_df$mean_mm_length <= 11.2, 2,
                                           ifelse(nmds1_df$mean_mm_length > 11.2 & nmds1_df$mean_mm_length <= 13.775, 3, 4
                                           ))))


nmds1_df$quantiles <- factor(nmds1_df$quantiles, levels = c("0", "1", "2", "3", "4"))
ellipse3<-ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=quantiles))+
  geom_point(size = 6)+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()
ellipse3
#########################
#Size within species experiment####
###########################

HEV <- dna %>%
  filter(sample_str == "HEV") %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  group_by(sample, Family, pred_ID, reads, Method, Island, 
           Habitat, Microhabitat, Year, No..Individuals) %>%
  summarise(mean_mm_length = mean(Length_mm), 
            total = n(), 
            sd = sd(Length_mm), 
            se=sd/sqrt(total))

HEV_meta <- HEV %>%
  ungroup() %>%
  dplyr::select(sample, pred_ID, mean_mm_length) %>%
  distinct(sample, pred_ID, mean_mm_length)


HEV_dat <- HEV %>%
  ungroup() %>%
  dplyr::select(sample, Family, reads) %>%
  group_by(sample, Family) %>%
  summarise(reads = sum(reads)) %>%
  filter(!Family %in% c("Muridae", "Hominidae")) %>%
  pivot_wider(names_from = Family,
              values_from = reads) %>%
  column_to_rownames(var = "sample") %>%
  replace(., is.na(.), 0)

zeros_HEV <- HEV_dat[which(rowSums(HEV_dat) == 0),] #9

zeros_HEV$sample <- rownames(zeros_HEV)


HEV_dat <- HEV_dat[, which(colSums(HEV_dat) != 0)] #66 to 30
HEV_dat <- HEV_dat[which(rowSums(HEV_dat) != 0),] #70 to 61
#I’m using the metaMDS function from the vegan package
nmds_hev <- metaMDS(HEV_dat, distance = "jaccard", binary = TRUE, k=2, trymax = 1000)
#then look at a stress plot. This is to evaluate whether the NMDS is actually a good representation of the structure of your data.
#I don’t know a ton about evaluating these, but it should be roughly a linear stair step
stressplot(nmds_hev)

HEV_meta <- HEV_meta %>%
  anti_join(zeros_HEV, by = "sample")


#making a data frame of the points from your NMDS
nmds_df<-data.frame(MDS1=nmds_hev$points[,1], MDS2=nmds_hev$points[,2])
#combining with metadata for plotting
nmds_df_meta<-cbind(HEV_meta, nmds_df)

#I’m gong to plot my points in different colors based on the “Site” that they came from, so I need to make sure that my “Site” is a factor, not numeric.
nmds_df_meta$pred_ID<-as.factor(nmds_df_meta$pred_ID)

#using ggplot to plot it
#only points
plotpoints<-ggplot(nmds_df_meta, aes(x=MDS1,y=MDS2, color = sample, size = mean_mm_length))+
  geom_point()+
  coord_fixed()+
  theme_bw()
plotpoints
