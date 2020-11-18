##########################
# 1. Is there stage structure? -----
# Ana Miller-ter Kuile
# October 20, 2020
###########################

# this script does k-means clustering to
#determine if there are clusters in three
#predator species for which we have the largest
#sample sizes. maybe extend to others too? not sure
#then looks at whether these stages are predicted by 
#body size. 


# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 
                  "glmmTMB", "emmeans",
                  "MuMIn", "DHARMa",
                  "effects", "ggeffects",
                  "vegan", "recluster",
                  "phytools", "stats",
                  "cluster", "data.tree")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Import data -------------------------------------------------------------

data <- read.csv(here("data", "outputs", "8_final_dataset",
                      "pred_prey_sizes_tp_DNAinteractions.csv"))

stage <- data %>%
  dplyr::select(-X, -X.x, -X.y, -reads) %>%
  mutate(pred_mass_mg = exp(pred_log_mass_mg)) %>%
  filter(sample_str %in% c("HEV", "NEO", "PHH"))

# HEV ---------------------------------------------------------------------
sppData <- split(stage, stage$sample_str)

process <- function(x, df = sppData){
#process <- function(df){
    
  mat <- df[[x]] %>% 
  #mat <- df %>% 
  dplyr::select(sample, Family) %>% #select sample and family
  group_by(sample, Family) %>% #group by them
  mutate(presence = 1) %>% #give presence value
  mutate(Family = factor(Family)) %>% #factor family
  mutate(presence = replace_na(presence, 0)) %>% #thus, allowing zeros to be set for NA
    ungroup() %>%
  pivot_wider(names_from = Family, #pivot to matrix where columns are diet
              values_from = presence, #and rows are samples
              values_fill = 0) %>% #fill any missing values with 0
  column_to_rownames(var = "sample") #set the sample to rowname for matrix format

  mat_cons <- mat[colSums(mat) > 1] #remove diet that only occurs once
  mat_cons <- mat_cons[rowSums(mat_cons) > 0,] #remove any samples that sets to 0
  
  #jaccard distance for these individuals
  dist <- recluster.dist(mat_cons, dist="jaccard")
  
  #using the most standard clustering algorithm, 
  #UPGMA = Unweighted Pair-Group Method Using Arithmetic Averages
  clust<-hclust(dist, "average")
  #plot(clust) #visualize that cluster
  cut <- cutree(clust, h = 0.2) #cut this cluster at 0.2, meaning total similarity
  dendro <- plot(clust, labels = as.character(cut)) #dendro plot
  cut_big <- cutree(clust, h = 0.75) #less conservative, but at least one species shared
  dendro_big <- plot(clust, labels = as.character(cut_big)) #dendro big plot
  
  cut_df <- as.data.frame(cutree(clust, h=0.2)) #makes conservative clusters a DF
  big_cut_df <- as.data.frame(cutree(clust, h=0.75)) #makes bigger clusters a DF
  
  cut_df_cs <- cut_df %>%
    rownames_to_column(var = "sample") %>% #make same a column again
    rename("cluster_cons" = "cutree(clust, h = 0.2)") %>% #rename the cluster column
    mutate(cluster_cons = as.factor(cluster_cons)) %>% #factor variable
    group_by(cluster_cons) %>% #group by it
    summarise(cluster_cons_n = n()) %>% #how many in cluster
    filter(cluster_cons_n > 1) %>% #remove clusters of one sample
    ungroup()
  
  cut_df <- cut_df %>%
    rownames_to_column(var = "sample") %>% #make sammple a column again
    rename("cluster_cons" = "cutree(clust, h = 0.2)") %>% #rename the cluster column
    mutate(cluster_cons = as.factor(cluster_cons)) %>% #factor variable
    left_join(cut_df_cs, by = "cluster_cons")
    
  big_cut_df_cs <- big_cut_df %>%
    rownames_to_column(var = "sample") %>% #make same a column again
    rename("cluster_big" = "cutree(clust, h = 0.75)") %>% #rename the cluster column
    mutate(cluster_big = as.factor(cluster_big)) %>% #factor variable
    group_by(cluster_big) %>% #group by it
    summarise(cluster_big_n = n()) %>% #how many samples in cluster
    filter(cluster_big_n > 1) %>% #remove clusters of one sample
    ungroup()
  
  big_cut_df <- big_cut_df %>%
    rownames_to_column(var = "sample") %>%
    rename("cluster_big" = "cutree(clust, h = 0.75)") %>%
    mutate(cluster_big = as.factor(cluster_big)) %>%
    left_join(big_cut_df_cs, by = "cluster_big")
  
  clustered_DF <- stage %>%
    left_join(cut_df, by = "sample") %>% #join back up with the full DF
    left_join(big_cut_df, by = "sample") %>% #join back to full df
    filter(!is.na(cluster_big))

  return(clustered_DF)
}

## Running on each data frame in the split list
out <- lapply(1:length(sppData), FUN = process)

stages_split <- bind_rows(out)


# Visualize ---------------------------------------------------------------

stages_split %>%
  filter(cluster_cons_n > 1) %>%
  ungroup() %>%
  distinct(sample, cluster_cons, pred_mass_mg, sample_str) %>%
  ggplot(aes(x = cluster_cons, y = pred_mass_mg)) +
  geom_boxplot() +
  facet_wrap(~sample_str, scale = "free") +
  theme_classic() +
  labs(x= "Cluster on 100% similarity", y = "Predator mass (mg)") +
  theme(axis.text.x = element_blank(),
        axis.title.y= element_text(size = 30),
        axis.title.x= element_text(size = 30, vjust = -1),
        axis.text.y = element_text(size = 25),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        strip.text = element_text(size = 25))

stages_split %>%
  filter(cluster_big_n > 1) %>%
  ungroup() %>%
  distinct(sample, cluster_big, pred_mass_mg, sample_str) %>%
  ggplot(aes(x = cluster_big, y = pred_mass_mg)) +
  geom_boxplot() +
  facet_wrap(~sample_str, scale = "free") +
  theme_classic() +
  labs(x= "Cluster on 25% similarity", y = "Predator mass (mg)") +
  theme(axis.text.x = element_blank(),
        axis.title.y= element_text(size = 30),
        axis.title.x= element_text(size = 30, vjust = -1),
        axis.text.y = element_text(size = 25),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        strip.text = element_text(size = 25))

stages_split %>%
  distinct(sample, pred_mass_mg, sample_str) %>%
  ggplot(aes(x = pred_mass_mg, fill = sample_str)) +
  geom_histogram(position = "identity") + 
  theme_classic() +
  scale_x_log10() +
  scale_y_continuous(breaks = c(2,4,6,8)) +
  facet_grid(sample_str~.) +
  theme(legend.position = "none",
        axis.text = element_text(size = 25),
        axis.title = element_text(size =30),
        strip.text = element_text(size = 15)) +
  labs(x = "Predator mass (mg)", y = "Count")






