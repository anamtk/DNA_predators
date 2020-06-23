###########################
#Species Stage Structure
#Ana Miller-ter Kuile
#June 22, 2020
###########################

#This code will look at how diet changes/doesn't change with
#species body size. 

#For this script, it is probably a good idea to start explorations
#with the species for which we have large sample sizes and/or
#large ranges of body size

#In this logic - Euborellia, Phisis, and Heteropoda are probably first

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
dna_fam <- read.csv(here("data", "outputs", "5_rarefied_taxonomic_sort", "fam_prey_DNA.csv"))

#sample metadata for sample sizes
meta <- read.csv(here("data", "Sample_metadata.csv"))


totals <- meta %>%
  distinct(ID, Extraction.ID, No..Individuals) %>%
  group_by(ID) %>%
  summarise(total = sum(No..Individuals))

meta <- meta %>%
  left_join(totals, by = "ID")

ggplot(meta, aes(x = Length_mm)) +
  geom_histogram() + theme_bw() +
  facet_wrap(~ID) 

meta <- meta %>%
  group_by(Method, Island, Habitat, Microhabitat, 
           Year, ID, Extraction.ID, No..Individuals,
           total) %>%
  summarise(avg_mm = mean(Length_mm))

#########################
#Heteropoda venatoria####
###########################
HEV <- dna_fam %>%
  filter(sample_str == "HEV")

HEV$sample <- str_sub(HEV$sample, start = 1, end =-2) #remove "a" at end

HEV <- HEV %>%
  left_join(meta, by= c("sample" = "Extraction.ID"))

HEV_meta <- HEV %>%
  dplyr::select(sample, Method, Island, Habitat, Microhabitat, Year, avg_mm)

HEV_dat <- HEV %>%
  dplyr::select(sample, reads, Family) %>%
  group_by(sample, Family) %>%
  summarise(reads = sum(reads)) %>%
  group_by(sample, Family) %>%
  summarise(reads = ifelse(reads > 0, 1, 0)) %>%
  pivot_wider(names_from = Family,
              values_from = reads) %>%
  column_to_rownames(var = "sample") %>%
  dplyr::select(-Muridae, Hominidae)

HEV_dat1 <- HEV_dat %>%
  dplyr::select(where(~ !is.numeric(.) || sum(.) >1 ))

zeros <- HEV_dat[which(rowSums(HEV_dat) == 0),] #12

zeros$sample <- rownames(zeros)

HEV_meta <- HEV_meta %>%
  anti_join(zeros, by = "sample") %>%
  distinct(sample, Method, Island, Habitat, Microhabitat, Year, avg_mm)

HEV_dat <- HEV_dat[, which(colSums(HEV_dat) != 0)] #70 to 31 Orders
HEV_dat <- HEV_dat[which(rowSums(HEV_dat) != 0),] #70 to 61 samples

#I’m using the metaMDS function from the vegan package
nmds1 <- metaMDS(HEV_dat, distance = "jaccard", binary=TRUE, k=2, trymax = 1000)
#then look at a stress plot. This is to evaluate whether the NMDS is actually a good representation of the structure of your data.
#I don’t know a ton about evaluating these, but it should be roughly a linear stair step
stressplot(nmds1)

#making a data frame of the points from your NMDS
nmds1_df<-data.frame(MDS1=nmds1$points[,1], MDS2=nmds1$points[,2])
#combining with metadata for plotting
nmds1_df<-cbind(HEV_meta, nmds1_df)

#using ggplot to plot it
#only points
ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=avg_mm))+
  geom_point(size = 3)+
  coord_fixed()+
  theme_bw()


ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=Habitat))+
  geom_point()+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()

ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=Microhabitat))+
  geom_point()+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()

nmds1_df$Year <- as.factor(nmds1_df$Year)
ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=Year))+
  geom_point()+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()

quantile(nmds1_df$avg_mm)
#    0%    25%    50%    75%   100% 
#2.500 12.025 16.800 20.925 27.500 

nmds1_df$quantiles <- ifelse(nmds1_df$avg_mm > 0 & nmds1_df$avg_mm <= 2.5, 0,
                             ifelse(nmds1_df$avg_mm > 2.5 & nmds1_df$avg_mm <= 12.025, 1,
                                    ifelse(nmds1_df$avg_mm > 12.025 & nmds1_df$avg_mm <= 16.800, 2,
                                           ifelse(nmds1_df$avg_mm > 16.800 & nmds1_df$avg_mm <= 20.925, 3, 4
                                           ))))


nmds1_df$quantiles <- factor(nmds1_df$quantiles, levels = c("0", "1", "2", "3", "4"))

ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=quantiles))+
  geom_point()+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()

#########################
#Euborellia####
###########################
EUB <- dna_fam %>%
  filter(sample_str == "EUB")

EUB$sample <- str_sub(EUB$sample, start = 1, end =-2) #remove "a" at end

EUB <- EUB %>%
  left_join(meta, by= c("sample" = "Extraction.ID"))

EUB_meta <- EUB %>%
  dplyr::select(sample, Method, Island, Habitat, Microhabitat, Year, avg_mm)

EUB_dat <- EUB %>%
  dplyr::select(sample, reads, Family) %>%
  group_by(sample, Family) %>%
  summarise(reads = sum(reads)) %>%
  group_by(sample, Family) %>%
  summarise(reads = ifelse(reads > 0, 1, 0)) %>%
  pivot_wider(names_from = Family,
              values_from = reads) %>%
  column_to_rownames(var = "sample") %>%
  dplyr::select(-Muridae, Hominidae)

EUB_dat1 <- EUB_dat %>%
  dplyr::select(where(~ !is.numeric(.) || sum(.) >1 ))

zeros <- EUB_dat1[which(rowSums(EUB_dat1) == 0),] #3

zeros$sample <- rownames(zeros)

EUB_meta <- EUB_meta %>%
  anti_join(zeros, by = "sample") %>%
  distinct(sample, Method, Island, Habitat, Microhabitat, Year, avg_mm)

EUB_dat <- EUB_dat[, which(colSums(EUB_dat) != 0)] #67 yo 8...
EUB_dat <- EUB_dat[which(rowSums(EUB_dat) != 0),] #22 to 19

#I’m using the metaMDS function from the vegan package
nmds1 <- metaMDS(EUB_dat, distance = "jaccard", binary = TRUE, k=2, trymax = 1000)
#then look at a stress plot. This is to evaluate whether the NMDS is actually a good representation of the structure of your data.
#I don’t know a ton about evaluating these, but it should be roughly a linear stair step
stressplot(nmds1)



#making a data frame of the points from your NMDS
nmds1_df<-data.frame(MDS1=nmds1$points[,1], MDS2=nmds1$points[,2])
#combining with metadata for plotting
nmds1_df<-cbind(EUB_meta, nmds1_df)

#using ggplot to plot it
#only points
ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=avg_mm))+
  geom_point(size = 3)+
  coord_fixed()+
  theme_bw()


ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=Habitat))+
  geom_point()+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()

ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=Microhabitat))+
  geom_point()+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()


quantile(nmds1_df$avg_mm)
#    0%    25%    50%    75%   100% 
#2.980  7.325  9.900 11.400 15.800

nmds1_df$quantiles <- ifelse(nmds1_df$avg_mm > 0 & nmds1_df$avg_mm <= 2.980, 0,
                             ifelse(nmds1_df$avg_mm > 2.980 & nmds1_df$avg_mm <= 7.325, 1,
                                    ifelse(nmds1_df$avg_mm >7.325 & nmds1_df$avg_mm <= 9.900, 2,
                                           ifelse(nmds1_df$avg_mm > 9.900 & nmds1_df$avg_mm <= 11.400, 3, 4
                                           ))))


nmds1_df$quantiles <- factor(nmds1_df$quantiles, levels = c("0", "1", "2", "3", "4"))

ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=quantiles))+
  geom_point()+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()

#########################
#Phisis####
###########################
PHH <- dna_fam %>%
  filter(sample_str == "PHH")

PHH$sample <- str_sub(PHH$sample, start = 1, end =-2) #remove "a" at end

PHH <- PHH %>%
  left_join(meta, by= c("sample" = "Extraction.ID"))

PHH_meta <- PHH %>%
  dplyr::select(sample, Method, Island, Habitat, Microhabitat, Year, avg_mm)

PHH_dat <- PHH %>%
  dplyr::select(sample, reads, Family) %>%
  group_by(sample, Family) %>%
  summarise(reads = sum(reads)) %>%
  group_by(sample, Family) %>%
  summarise(reads = ifelse(reads > 0, 1, 0)) %>%
  pivot_wider(names_from = Family,
              values_from = reads) %>%
  column_to_rownames(var = "sample") %>%
  dplyr::select(-Muridae, Hominidae) #54 of 67

PHH_dat1 <- PHH_dat %>%
  dplyr::select(where(~ !is.numeric(.) || sum(.) >1 )) #54 of 11

zeros <- PHH_dat[which(rowSums(PHH_dat) == 0),] #4

zeros$sample <- rownames(zeros)

PHH_meta <- PHH_meta %>%
  anti_join(zeros, by = "sample") %>%
  distinct(sample, Method, Island, Habitat, Microhabitat, Year, avg_mm)

PHH_dat <- PHH_dat[, which(colSums(PHH_dat) != 0)] #20
PHH_dat <- PHH_dat[which(rowSums(PHH_dat) != 0),] #50

#I’m using the metaMDS function from the vegan package
nmds1 <- metaMDS(PHH_dat, distance = "jaccard", binary = TRUE, k=2, trymax = 1000)
#then look at a stress plot. This is to evaluate whether the NMDS is actually a good representation of the structure of your data.
#I don’t know a ton about evaluating these, but it should be roughly a linear stair step
stressplot(nmds1)



#making a data frame of the points from your NMDS
nmds1_df<-data.frame(MDS1=nmds1$points[,1], MDS2=nmds1$points[,2])
#combining with metadata for plotting
nmds1_df<-cbind(PHH_meta, nmds1_df)

#using ggplot to plot it
#only points
ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=avg_mm))+
  geom_point(size = 3)+
  coord_fixed()+
  theme_bw()


ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=Habitat))+
  geom_point()+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()

ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=Microhabitat))+
  geom_point()+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()


quantile(nmds1_df$avg_mm)
#    0%    25%    50%    75%   100% 
#6.35 11.20 13.10 14.30 17.70 

nmds1_df$quantiles <- ifelse(nmds1_df$avg_mm > 0 & nmds1_df$avg_mm <= 6.35, 0,
                             ifelse(nmds1_df$avg_mm > 6.35 & nmds1_df$avg_mm <= 11.20, 1,
                                    ifelse(nmds1_df$avg_mm >11.20 & nmds1_df$avg_mm <= 13.10, 2,
                                           ifelse(nmds1_df$avg_mm > 13.10 & nmds1_df$avg_mm <= 14.30, 3, 4
                                           ))))


nmds1_df$quantiles <- factor(nmds1_df$quantiles, levels = c("0", "1", "2", "3", "4"))

ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=quantiles))+
  geom_point()+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()

#########################
#Neoscona####
###########################
NEO <- dna_fam %>%
  filter(sample_str == "NEO")

NEO$sample <- str_sub(NEO$sample, start = 1, end =-2) #remove "a" at end

NEO <- NEO %>%
  left_join(meta, by= c("sample" = "Extraction.ID"))

NEO_meta <- NEO %>%
  dplyr::select(sample, Method, Island, Habitat, Microhabitat, Year, avg_mm)%>%
  distinct(sample, Method, Island, Habitat, Microhabitat, Year, avg_mm)

NEO_dat <- NEO %>%
  dplyr::select(sample, reads, Family) %>%
  group_by(sample, Family) %>%
  summarise(reads = sum(reads)) %>%
  group_by(sample, Family) %>%
  summarise(reads = ifelse(reads > 0, 1, 0)) %>%
  pivot_wider(names_from = Family,
              values_from = reads) %>%
  column_to_rownames(var = "sample") %>%
  dplyr::select(-Muridae, Hominidae) #26 of 66

NEO_dat1 <- NEO_dat %>%
  dplyr::select(where(~ !is.numeric(.) || sum(.) >1 )) #26 of 15

zeros <- NEO_dat[which(rowSums(NEO_dat) == 0),] #0

NEO_dat <- NEO_dat[, which(colSums(NEO_dat) != 0)] #66 to 30
NEO_dat <- NEO_dat[which(rowSums(NEO_dat) != 0),] #none 

#I’m using the metaMDS function from the vegan package
nmds1 <- metaMDS(NEO_dat, distance = "jaccard", binary = TRUE, k=2, trymax = 1000)
#then look at a stress plot. This is to evaluate whether the NMDS is actually a good representation of the structure of your data.
#I don’t know a ton about evaluating these, but it should be roughly a linear stair step
stressplot(nmds1)



#making a data frame of the points from your NMDS
nmds1_df<-data.frame(MDS1=nmds1$points[,1], MDS2=nmds1$points[,2])
#combining with metadata for plotting
nmds1_df<-cbind(NEO_meta, nmds1_df)

#using ggplot to plot it
#only points
ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=avg_mm))+
  geom_point(size = 3)+
  coord_fixed()+
  theme_bw()


ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=Habitat))+
  geom_point()+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()

ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=Microhabitat))+
  geom_point()+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()


quantile(nmds1_df$avg_mm)
#    0%    25%    50%    75%   100% 
#1.700000 4.033333 5.700000 7.400000 9.300000 

nmds1_df$quantiles <- ifelse(nmds1_df$avg_mm > 0 & nmds1_df$avg_mm <= 1.7, 0,
                             ifelse(nmds1_df$avg_mm > 1.7 & nmds1_df$avg_mm <= 4.03, 1,
                                    ifelse(nmds1_df$avg_mm >4.03 & nmds1_df$avg_mm <= 5.7, 2,
                                           ifelse(nmds1_df$avg_mm > 5.7 & nmds1_df$avg_mm <= 7.4, 3, 4
                                           ))))


nmds1_df$quantiles <- factor(nmds1_df$quantiles, levels = c("0", "1", "2", "3", "4"))

ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=quantiles))+
  geom_point()+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()

#########################
#SCYtodes####
###########################
SCY <- dna_fam %>%
  filter(sample_str == "SCY")

SCY$sample <- str_sub(SCY$sample, start = 1, end =-2) #remove "a" at end

SCY <- SCY %>%
  left_join(meta, by= c("sample" = "Extraction.ID"))

SCY_meta <- SCY %>%
  dplyr::select(sample, Method, Island, Habitat, Microhabitat, Year, avg_mm)%>%
  distinct(sample, Method, Island, Habitat, Microhabitat, Year, avg_mm)

SCY_dat <- SCY %>%
  dplyr::select(sample, reads, Family) %>%
  group_by(sample, Family) %>%
  summarise(reads = sum(reads)) %>%
  group_by(sample, Family) %>%
  summarise(reads = ifelse(reads > 0, 1, 0)) %>%
  pivot_wider(names_from = Family,
              values_from = reads) %>%
  column_to_rownames(var = "sample") %>%
  dplyr::select(-Muridae, Hominidae) #26 of 66

SCY_dat1 <- SCY_dat %>%
  dplyr::select(where(~ !is.numeric(.) || sum(.) >1 )) #26 of 15

zeros <- SCY_dat[which(rowSums(SCY_dat) == 0),] #0

SCY_dat <- SCY_dat[, which(colSums(SCY_dat) != 0)] #67 to 15
SCY_dat <- SCY_dat[which(rowSums(SCY_dat) != 0),] #none 

#I’m using the metaMDS function from the vegan package
nmds1 <- metaMDS(SCY_dat, distance = "jaccard", binary = TRUE, k=2, trymax = 1000)
#then look at a stress plot. This is to evaluate whether the NMDS is actually a good representation of the structure of your data.
#I don’t know a ton about evaluating these, but it should be roughly a linear stair step
stressplot(nmds1)



#making a data frame of the points from your NMDS
nmds1_df<-data.frame(MDS1=nmds1$points[,1], MDS2=nmds1$points[,2])
#combining with metadata for plotting
nmds1_df<-cbind(SCY_meta, nmds1_df)

#using ggplot to plot it
#only points
ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=avg_mm))+
  geom_point(size = 3)+
  coord_fixed()+
  theme_bw()


ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=Habitat))+
  geom_point()+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()

ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=Microhabitat))+
  geom_point()+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()


quantile(nmds1_df$avg_mm)
#    0%    25%    50%    75%   100% 
#2.366667  3.333333  5.480000  9.400000 11.100000 

nmds1_df$quantiles <- ifelse(nmds1_df$avg_mm > 0 & nmds1_df$avg_mm <= 2.366667 , 0,
                             ifelse(nmds1_df$avg_mm > 2.366667  & nmds1_df$avg_mm <= 3.333333, 1,
                                    ifelse(nmds1_df$avg_mm >3.333333 & nmds1_df$avg_mm <= 5.480000, 2,
                                           ifelse(nmds1_df$avg_mm > 5.480000 & nmds1_df$avg_mm <= 9.400000 , 3, 4
                                           ))))


nmds1_df$quantiles <- factor(nmds1_df$quantiles, levels = c("0", "1", "2", "3", "4"))

ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=quantiles))+
  geom_point()+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()

#########################
#CEN####
###########################
CEN <- dna_fam %>%
  filter(sample_str == "CEN")

CEN$sample <- str_sub(CEN$sample, start = 1, end =-2) #remove "a" at end

CEN <- CEN %>%
  left_join(meta, by= c("sample" = "Extraction.ID"))

CEN_meta <- CEN %>%
  dplyr::select(sample, Method, Island, Habitat, Microhabitat, Year, avg_mm)%>%
  distinct(sample, Method, Island, Habitat, Microhabitat, Year, avg_mm)

CEN_dat <- CEN %>%
  dplyr::select(sample, reads, Family) %>%
  group_by(sample, Family) %>%
  summarise(reads = sum(reads)) %>%
  group_by(sample, Family) %>%
  summarise(reads = ifelse(reads > 0, 1, 0)) %>%
  pivot_wider(names_from = Family,
              values_from = reads) %>%
  column_to_rownames(var = "sample") %>%
  dplyr::select(-Muridae, Hominidae) #26 of 66

CEN_dat1 <- CEN_dat %>%
  dplyr::select(where(~ !is.numeric(.) || sum(.) >1 )) #26 of 15

zeros <- CEN_dat[which(rowSums(CEN_dat) == 0),] #0

CEN_dat <- CEN_dat[, which(colSums(CEN_dat) != 0)] #65 to 7
CEN_dat <- CEN_dat[which(rowSums(CEN_dat) != 0),] #none 

#I’m using the metaMDS function from the vegan package
nmds1 <- metaMDS(CEN_dat, distance = "jaccard", binary = TRUE, k=2, trymax = 1000)
#then look at a stress plot. This is to evaluate whether the NMDS is actually a good representation of the structure of your data.
#I don’t know a ton about evaluating these, but it should be roughly a linear stair step
stressplot(nmds1)



#making a data frame of the points from your NMDS
nmds1_df<-data.frame(MDS1=nmds1$points[,1], MDS2=nmds1$points[,2])
#combining with metadata for plotting
nmds1_df<-cbind(CEN_meta, nmds1_df)

#using ggplot to plot it
#only points
ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=avg_mm))+
  geom_point(size = 3)+
  coord_fixed()+
  theme_bw()


ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=Habitat))+
  geom_point()+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()

ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=Microhabitat))+
  geom_point()+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()


quantile(nmds1_df$avg_mm)
#    0%    25%    50%    75%   100% 
#15.000 19.825 25.350 36.500 40.800

nmds1_df$quantiles <- ifelse(nmds1_df$avg_mm > 0 & nmds1_df$avg_mm <= 15.000 , 0,
                             ifelse(nmds1_df$avg_mm > 15.000  & nmds1_df$avg_mm <= 19.825, 1,
                                    ifelse(nmds1_df$avg_mm >19.825 & nmds1_df$avg_mm <= 25.350, 2,
                                           ifelse(nmds1_df$avg_mm > 25.350 & nmds1_df$avg_mm <= 36.500 , 3, 4
                                           ))))


nmds1_df$quantiles <- factor(nmds1_df$quantiles, levels = c("0", "1", "2", "3", "4"))

ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=quantiles))+
  geom_point()+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()

#########################
#SME####
###########################
SME <- dna_fam %>%
  filter(sample_str == "SME")

SME$sample <- str_sub(SME$sample, start = 1, end =-2) #remove "a" at end

SME <- SME %>%
  left_join(meta, by= c("sample" = "Extraction.ID"))

SME_meta <- SME %>%
  dplyr::select(sample, Method, Island, Habitat, Microhabitat, Year, avg_mm)%>%
  distinct(sample, Method, Island, Habitat, Microhabitat, Year, avg_mm)

SME_dat <- SME %>%
  dplyr::select(sample, reads, Family) %>%
  group_by(sample, Family) %>%
  summarise(reads = sum(reads)) %>%
  group_by(sample, Family) %>%
  summarise(reads = ifelse(reads > 0, 1, 0)) %>%
  pivot_wider(names_from = Family,
              values_from = reads) %>%
  column_to_rownames(var = "sample") %>%
  dplyr::select(-Muridae, Hominidae) #26 of 66

SME_dat1 <- SME_dat %>%
  dplyr::select(where(~ !is.numeric(.) || sum(.) >1 )) #26 of 15

zeros <- SME_dat[which(rowSums(SME_dat) == 0),] #0

SME_dat <- SME_dat[, which(colSums(SME_dat) != 0)] #65 to 9
SME_dat <- SME_dat[which(rowSums(SME_dat) != 0),] #none 

#I’m using the metaMDS function from the vegan package
nmds1 <- metaMDS(SME_dat, distance = "jaccard", binary = TRUE, k=2, trymax = 1000)
#then look at a stress plot. This is to evaluate whether the NMDS is actually a good representation of the structure of your data.
#I don’t know a ton about evaluating these, but it should be roughly a linear stair step
stressplot(nmds1)



#making a data frame of the points from your NMDS
nmds1_df<-data.frame(MDS1=nmds1$points[,1], MDS2=nmds1$points[,2])
#combining with metadata for plotting
nmds1_df<-cbind(SME_meta, nmds1_df)

#using ggplot to plot it
#only points
ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=avg_mm))+
  geom_point(size = 3)+
  coord_fixed()+
  theme_bw()


ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=Habitat))+
  geom_point()+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()

ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=Microhabitat))+
  geom_point()+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()


quantile(nmds1_df$avg_mm)
#    0%    25%    50%    75%   100% 
#15.000 19.825 25.350 36.500 40.800

nmds1_df$quantiles <- ifelse(nmds1_df$avg_mm > 0 & nmds1_df$avg_mm <= 15.000 , 0,
                             ifelse(nmds1_df$avg_mm > 15.000  & nmds1_df$avg_mm <= 19.825, 1,
                                    ifelse(nmds1_df$avg_mm >19.825 & nmds1_df$avg_mm <= 25.350, 2,
                                           ifelse(nmds1_df$avg_mm > 25.350 & nmds1_df$avg_mm <= 36.500 , 3, 4
                                           ))))


nmds1_df$quantiles <- factor(nmds1_df$quantiles, levels = c("0", "1", "2", "3", "4"))

ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=quantiles))+
  geom_point()+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()

#########################
#LRs####
###########################
LRS <- dna_fam %>%
  filter(sample_str == "LRS")

LRS$sample <- str_sub(LRS$sample, start = 1, end =-2) #remove "a" at end

LRS <- LRS %>%
  left_join(meta, by= c("sample" = "Extraction.ID"))

LRS_meta <- LRS %>%
  dplyr::select(sample, Method, Island, Habitat, Microhabitat, Year, avg_mm)%>%
  distinct(sample, Method, Island, Habitat, Microhabitat, Year, avg_mm)

LRS_dat <- LRS %>%
  dplyr::select(sample, reads, Family) %>%
  group_by(sample, Family) %>%
  summarise(reads = sum(reads)) %>%
  group_by(sample, Family) %>%
  summarise(reads = ifelse(reads > 0, 1, 0)) %>%
  pivot_wider(names_from = Family,
              values_from = reads) %>%
  column_to_rownames(var = "sample") %>%
  dplyr::select(-Muridae, Hominidae) #26 of 66

LRS_dat1 <- LRS_dat %>%
  dplyr::select(where(~ !is.numeric(.) || sum(.) >1 )) #26 of 15

zeros <- LRS_dat[which(rowSums(LRS_dat) == 0),] #0

LRS_dat <- LRS_dat[, which(colSums(LRS_dat) != 0)] #67 to 8
LRS_dat <- LRS_dat[which(rowSums(LRS_dat) != 0),] #none 

#I’m using the metaMDS function from the vegan package
nmds1 <- metaMDS(LRS_dat, distance = "jaccard", binary = TRUE, k=2, trymax = 1000)
#then look at a stress plot. This is to evaluate whether the NMDS is actually a good representation of the structure of your data.
#I don’t know a ton about evaluating these, but it should be roughly a linear stair step
stressplot(nmds1)



#making a data frame of the points from your NMDS
nmds1_df<-data.frame(MDS1=nmds1$points[,1], MDS2=nmds1$points[,2])
#combining with metadata for plotting
nmds1_df<-cbind(LRS_meta, nmds1_df)

#using ggplot to plot it
#only points
ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=avg_mm))+
  geom_point(size = 3)+
  coord_fixed()+
  theme_bw()


ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=Habitat))+
  geom_point()+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()

ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=Microhabitat))+
  geom_point()+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()


quantile(nmds1_df$avg_mm)
#    0%    25%    50%    75%   100% 
#1.233333 1.677083 1.872500 1.927500 1.950000

nmds1_df$quantiles <- ifelse(nmds1_df$avg_mm > 0 & nmds1_df$avg_mm <= 1.233333, 0,
                             ifelse(nmds1_df$avg_mm > 1.233333 & nmds1_df$avg_mm <= 1.677083, 1,
                                    ifelse(nmds1_df$avg_mm >1.677083 & nmds1_df$avg_mm <= 1.872500, 2,
                                           ifelse(nmds1_df$avg_mm > 1.872500 & nmds1_df$avg_mm <= 1.927500 , 3, 4
                                           ))))


nmds1_df$quantiles <- factor(nmds1_df$quantiles, levels = c("0", "1", "2", "3", "4"))

ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=quantiles))+
  geom_point()+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()


#########################
#PAN####
###########################
PAN <- dna_fam %>%
  filter(sample_str == "PAN")

PAN$sample <- str_sub(PAN$sample, start = 1, end =-2) #remove "a" at end

PAN <- PAN %>%
  left_join(meta, by= c("sample" = "Extraction.ID"))

PAN_meta <- PAN %>%
  dplyr::select(sample, Method, Island, Habitat, Microhabitat, Year, avg_mm)%>%
  distinct(sample, Method, Island, Habitat, Microhabitat, Year, avg_mm)

PAN_dat <- PAN %>%
  dplyr::select(sample, reads, Family) %>%
  group_by(sample, Family) %>%
  summarise(reads = sum(reads)) %>%
  group_by(sample, Family) %>%
  summarise(reads = ifelse(reads > 0, 1, 0)) %>%
  pivot_wider(names_from = Family,
              values_from = reads) %>%
  column_to_rownames(var = "sample") %>%
  dplyr::select(-Muridae, Hominidae) #26 of 66

PAN_dat1 <- PAN_dat %>%
  dplyr::select(where(~ !is.numeric(.) || sum(.) >1 )) #26 of 15

zeros <- PAN_dat[which(rowSums(PAN_dat) == 0),] #0
zeros$sample <- rownames(zeros)

PAN_meta <- PAN_meta %>%
  anti_join(zeros, by = "sample")

PAN_dat <- PAN_dat[, which(colSums(PAN_dat) != 0)] #65 to 9
PAN_dat <- PAN_dat[which(rowSums(PAN_dat) != 0),] #none 

#I’m using the metaMDS function from the vegan package
nmds1 <- metaMDS(PAN_dat, distance = "jaccard", binary = TRUE, k=2, trymax = 1000)
#then look at a stress plot. This is to evaluate whether the NMDS is actually a good representation of the structure of your data.
#I don’t know a ton about evaluating these, but it should be roughly a linear stair step
stressplot(nmds1)



#making a data frame of the points from your NMDS
nmds1_df<-data.frame(MDS1=nmds1$points[,1], MDS2=nmds1$points[,2])
#combining with metadata for plotting
nmds1_df<-cbind(PAN_meta, nmds1_df)

#using ggplot to plot it
#only points
ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=avg_mm))+
  geom_point(size = 3)+
  coord_fixed()+
  theme_bw()


ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=Habitat))+
  geom_point()+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()

ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=Microhabitat))+
  geom_point()+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()


quantile(nmds1_df$avg_mm)
#    0%    25%    50%    75%   100% 
#34.6 38.0 40.2 41.0 44.1 

nmds1_df$quantiles <- ifelse(nmds1_df$avg_mm > 0 & nmds1_df$avg_mm <= 34.6 , 0,
                             ifelse(nmds1_df$avg_mm > 34.6  & nmds1_df$avg_mm <= 38.0, 1,
                                    ifelse(nmds1_df$avg_mm >38.0 & nmds1_df$avg_mm <= 40.2, 2,
                                           ifelse(nmds1_df$avg_mm > 40.2 & nmds1_df$avg_mm <= 41.0  , 3, 4
                                           ))))


nmds1_df$quantiles <- factor(nmds1_df$quantiles, levels = c("0", "1", "2", "3", "4"))

ggplot(nmds1_df, aes(x=MDS1,y=MDS2, color=quantiles))+
  geom_point()+
  stat_ellipse(level = 0.90)+
  coord_fixed()+
  theme_bw()



