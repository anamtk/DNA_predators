##########################
# 1. Is there stage structure? -----
# Ana Miller-ter Kuile
# October 20, 2020
###########################

# this script does hierarchical clustering to
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
                  "cluster", "data.tree",
                  "ggbeeswarm", "ggdendro")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Import data -------------------------------------------------------------

data <- read.csv(here("data", 
                      "outputs",  
                      "8_final_dataset", 
                      "pred_prey_sizes_DNAinteractions.csv"))

stage <- data %>%
  dplyr::select(-X, -X.x, -X.y, -reads) %>%
  mutate(pred_mass_mg = exp(pred_log_mass_mg)) %>%
  filter(sample_str %in% c("HEV", "NEO", "PHH"))

stage %>%
  distinct(sample, sample_str) %>%
  group_by(sample_str) %>%
  summarise(n = n())

# 100% clusters ---------------------------------------------------------------------
cons_clusters <- function(x, df = sppData){
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
 
  cut_df <- as.data.frame(cutree(clust, h=0.2)) #makes conservative clusters a DF

  cut_df_cs <- cut_df %>%
    rownames_to_column(var = "sample") %>% #make same a column again
    rename("cluster_cons" = "cutree(clust, h = 0.2)") %>% #rename the cluster column
    mutate(cluster_cons = as.factor(cluster_cons)) %>% #factor variable
    group_by(cluster_cons) %>% #group by it
    summarise(cluster_cons_n = n()) %>% #how many in cluster
    filter(cluster_cons_n > 2) %>% #remove clusters of 2 or fewer
    ungroup()
  
  cut_df <- cut_df %>%
    rownames_to_column(var = "sample") %>% #make sammple a column again
    rename("cluster_cons" = "cutree(clust, h = 0.2)") %>% #rename the cluster column
    mutate(cluster_cons = as.factor(cluster_cons)) %>% #factor variable
    left_join(cut_df_cs, by = "cluster_cons")
    
  clustered_DF <- stage %>%
    left_join(cut_df, by = "sample") %>% #join back up with the full DF
    ungroup() %>%
    distinct(sample, 
             cluster_cons,  
             pred_mass_mg, 
             sample_str,
             cluster_cons_n) %>%
    filter(!is.na(cluster_cons_n)) #remove any with 1 or fewere indivs

  return(clustered_DF)
}

# 50% clusters ------------------------------------------------------------

big_clusters <- function(x, df = sppData){
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
  cut_big <- cutree(clust, h = 0.5) #less conservative, but at least one species shared
  dendro_big <- plot(clust, labels = as.character(cut_big)) #dendro big plot
  
  big_cut_df <- as.data.frame(cutree(clust, h=0.5)) #makes bigger clusters a DF

  big_cut_df_cs <- big_cut_df %>%
    rownames_to_column(var = "sample") %>% #make same a column again
    rename("cluster_big" = "cutree(clust, h = 0.5)") %>% #rename the cluster column
    mutate(cluster_big = as.factor(cluster_big)) %>% #factor variable
    group_by(cluster_big) %>% #group by it
    summarise(cluster_big_n = n()) %>% #how many samples in cluster
    filter(cluster_big_n > 2) %>% #remove clusters of 2 or fewer
    ungroup()
  
  big_cut_df <- big_cut_df %>%
    rownames_to_column(var = "sample") %>%
    rename("cluster_big" = "cutree(clust, h = 0.5)") %>%
    mutate(cluster_big = as.factor(cluster_big)) %>%
    left_join(big_cut_df_cs, by = "cluster_big")
  
  clustered_DF <- stage %>%
    left_join(big_cut_df, by = "sample") %>% #join back to full df
    filter(!is.na(cluster_big)) %>%
    ungroup() %>%
    distinct(sample,  
             cluster_big, 
             pred_mass_mg, 
             sample_str,
             cluster_big_n) %>%
    filter(!is.na(cluster_big_n)) #remove any with one or fewere indivs
  
  return(clustered_DF)
}

# Run functions -----------------------------------------------------------
#split data by species
sppData <- split(stage, stage$sample_str)
## Running on each data frame in the split list
cluster_cons <- lapply(1:length(sppData), FUN = cons_clusters)

## Running on each data frame in the split list
cluster_big <- lapply(1:length(sppData), FUN = big_clusters)

# GLM per species -------------------------------------------------------
#DF per species for statistics
HEV_cons <- cluster_cons[[1]]
NEO_cons <- cluster_cons[[2]]
PHH_cons <- cluster_cons[[3]]

HEV_big <- cluster_big[[1]]
NEO_big <- cluster_big[[2]]
PHH_big <- cluster_big[[3]]

hist(HEV_big$pred_mass_mg)
hist(HEV_cons$pred_mass_mg)

m <- lm(pred_mass_mg ~ cluster_cons,
               data = HEV_cons)

summary(m)
fit <- simulateResiduals(m, plot=T)
em <- emmeans(m, "cluster_cons")
pairs(em)
tukey <- as.data.frame(pairs(em))

tukey <- tukey %>%
  mutate(sig = ifelse(p.value < 0.05, "sig", "non-sig"))

HEV_cons %>%
  mutate(cluster_cons = fct_reorder(cluster_cons, pred_mass_mg, .fun='mean')) %>%
  ggplot(aes(x = cluster_cons, y = pred_mass_mg)) +
  geom_boxplot() +
  theme_bw() +
  scale_y_log10() +
  labs(x= "Cluster on % similarity", y = "Predator mass (mg)") +
  theme(axis.title.y= element_text(size = 30),
        axis.title.x= element_text(size = 30, vjust = -1),
        axis.text.y = element_text(size = 25))

#pairwise comparisons
tukey %>%
  arrange(estimate, p.value) %>%
  mutate(contrast = factor(contrast, levels = contrast)) %>%
  ggplot(aes(x = contrast, y = estimate, color = sig)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_point(size = 2) +
  geom_linerange(aes(ymin = estimate - SE, 
                     ymax = estimate + SE), 
                 size = 1) +
  scale_color_manual(values = c("sig" = "#229AAA",
                                "non-sig" = "#878787")) +
  theme_bw() +
  coord_flip() +
  labs(x = "Pairwise contrast", y = "Difference") +
  theme(legend.position = "none",
        axis.title = element_text(size = 25))

# ONe DF ---------------------------------------------------------------


#get one DF as opposed to list
cons_stages <- bind_rows(cluster_cons)
#how many clustered:
cons_stages %>%
  group_by(sample_str) %>%
  summarise(n = n())
#HEV: 29/53 
#NEO: 5/24
#PHH: 25/42
#get one DF as opposed to list
big_stages <- bind_rows(cluster_big)
#how many clustered:
big_stages %>%
  group_by(sample_str) %>%
  summarise(n = n())
#HEV: 49/53 
#NEO: 21/24
#PHH: 40/42

# Visualize ---------------------------------------------------------------

cons_stages %>%
  mutate(cluster_cons = fct_reorder(cluster_cons, pred_mass_mg, .fun='mean')) %>%
  ggplot(aes(x = cluster_cons, y = pred_mass_mg)) +
  geom_boxplot() +
  facet_wrap(~sample_str, scale = "free") +
  theme_classic() +
  labs(x= "Cluster on 100% similarity", y = "Predator mass (mg)") +
  theme(axis.title.y= element_text(size = 30),
        axis.title.x= element_text(size = 30, vjust = -1),
        axis.text.y = element_text(size = 25),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        strip.text = element_text(size = 25))

big_stages %>%
  mutate(cluster_big = fct_reorder(cluster_big, pred_mass_mg, .fun='mean')) %>%
  ggplot(aes(x = cluster_big, y = pred_mass_mg)) +
  geom_boxplot() +
  facet_wrap(~sample_str, scale = "free") +
  theme_classic() +
  labs(x= "Cluster on 50% similarity", y = "Predator mass (mg)") +
  theme(axis.text.x = element_blank(),
        axis.title.y= element_text(size = 30),
        axis.title.x= element_text(size = 30, vjust = -1),
        axis.text.y = element_text(size = 25),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        strip.text = element_text(size = 25))


prey_richness <- stage %>%
  distinct(sample, sample_str, Family) %>%
  group_by(sample, sample_str) %>%
  summarise(richness = n()) 



ggplot(aes(x = sample_str, y = richness)) +
  geom_boxplot() +
  geom_beeswarm() +
  theme_bw() 

plot(clust, labels = as.character(cut))
ggdendrogram(clust) +
  theme_void() +
  labs(x = "Predator individual", y = "Jaccard similarity") +
  theme(axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 20),
        axis.title.y = element_text(angle = 90),
        axis.text.x = element_text(size = 15, angle = 90))



