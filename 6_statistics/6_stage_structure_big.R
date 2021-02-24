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
                  "ggbeeswarm", "ggdendro",
                  "patchwork")

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

prey_richness <- stage %>%
  distinct(sample, sample_str, Family) %>%
  group_by(sample, sample_str) %>%
  summarise(richness = n()) 

# Dendrogram Function -----------------------------------------------------

dendro <- function(x, df = sppData){
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
  
  return(clust)
  
}

# 50% clusters function ---------------------------------------------------------------------
big_clusters <- function(x, hclust = cluster_objs){
  #process <- function(df){
  cluster <- hclust[[x]]
  
  cut_df <- as.data.frame(cutree(cluster, h=0.5)) #makes conservative clusters a DF
  
  cut_df_cs <- cut_df %>%
    rownames_to_column(var = "sample") %>% #make same a column again
    rename("cluster_big" = "cutree(cluster, h = 0.5)") %>% #rename the cluster column
    mutate(cluster_big = as.factor(cluster_big)) %>% #factor variable
    group_by(cluster_big) %>% #group by it
    summarise(cluster_big_n = n()) %>% #how many in cluster
    filter(cluster_big_n > 1) %>% #remove clusters of 2 or fewer
    ungroup()
  
  cut_df <- cut_df %>%
    rownames_to_column(var = "sample") %>% #make sammple a column again
    rename("cluster_big" = "cutree(cluster, h = 0.5)") %>% #rename the cluster column
    mutate(cluster_big = as.factor(cluster_big)) %>% #factor variable
    left_join(cut_df_cs, by = "cluster_big")
  
  clustered_DF <- stage %>%
    left_join(cut_df, by = "sample") %>% #join back up with the full DF
    ungroup() %>%
    distinct(sample, 
             cluster_big,  
             pred_mass_mg, 
             sample_str,
             cluster_big_n) %>%
    filter(!is.na(cluster_big_n)) #remove any with 1 or fewere indivs
  
  return(clustered_DF)
}

# Run functions -----------------------------------------------------------
#split data by species
sppData <- split(stage, stage$sample_str)
## Running on each data frame in the split list
cluster_objs <- lapply(1:length(sppData), FUN = dendro)

cluster_dfs <- lapply(1:length(cluster_objs), FUN = big_clusters)

# GLM per species -------------------------------------------------------
#DF per species for statistics
HEV_50 <- cluster_dfs[[1]]
NEO_50 <- cluster_dfs[[2]]
PHH_50 <- cluster_dfs[[3]]

pred_cols <- c("HEV" = "#EEB00C", "NEO" = "#158D8E", "PHH" = "#114C54")
# HEV ---------------------------------------------------------------------

hist(HEV_50$pred_mass_mg)

m_h <- glm(pred_mass_mg ~ cluster_big,
           data = HEV_50,
           na.action= "na.fail")
dredge(m_h)

m_h2 <- glm(pred_mass_mg ~ 1,
           data = HEV_100,
           na.action= "na.fail")

fit <- simulateResiduals(m_h2, plot=T)

#NEO: "#EEB00C" 
#PHH: "#158D8E"
#HEV: "#114C54"

clusters_HEV <- HEV_50 %>%
  mutate(cluster_big = fct_reorder(cluster_big, pred_mass_mg, .fun='mean')) %>%
  ggplot(aes(x = cluster_big, y = pred_mass_mg)) +
  geom_boxplot(color = "#114C54", size = 1) +
  geom_point(color = "#114C54") +
  theme_bw() +
  scale_y_log10() +
  labs(x= "50% similarity clusters", y = "Predator mass (mg)") +
  theme(axis.title.y= element_text(size = 20),
        axis.title.x= element_text(size = 20),
        axis.text = element_text(size = 15))

# NEO ---------------------------------------------------------------------
hist(NEO_50$pred_mass_mg)

m_n <- lm(pred_mass_mg ~ cluster_big,
          data = NEO_50,
          na.action = "na.fail")
dredge(m_n)
m_n2 <- lm(pred_mass_mg ~ 1,
          data = NEO_50,
          na.action = "na.fail")
fit <- simulateResiduals(m_n2, plot=T)

#NEO: "#EEB00C" 
#PHH: "#158D8E"
#HEV: "#114C54"

clusters_NEO <- NEO_50 %>%
  mutate(cluster_big = fct_reorder(cluster_big, pred_mass_mg, .fun='mean')) %>%
  ggplot(aes(x = cluster_big, y = pred_mass_mg)) +
  geom_boxplot(color = "#EEB00C", size = 1) +
  geom_point(color = "#EEB00C") +
  theme_bw() +
  scale_y_log10() +
  labs(x= "50% similarity clusters", y = "Predator mass (mg)") +
  theme(axis.title.y= element_text(size = 20),
        axis.title.x= element_text(size = 20),
        axis.text = element_text(size = 15))

# PHH ---------------------------------------------------------------------

hist(PHH_50$pred_mass_mg)

m_p <- lm(pred_mass_mg ~ cluster_big,
          data = PHH_50,
          na.action = "na.fail")
dredge(m_p)

m_p2 <- lm(pred_mass_mg ~ 1,
          data = PHH_50,
          na.action = "na.fail")

fit <- simulateResiduals(m_p2, plot=T)

#NEO: "#EEB00C" 
#PHH: "#158D8E"
#HEV: "#114C54"

clusters_PHH <- PHH_50 %>%
  mutate(cluster_big = fct_reorder(cluster_big, pred_mass_mg, .fun='mean')) %>%
  ggplot(aes(x = cluster_big, y = pred_mass_mg)) +
  geom_boxplot(color = "#158D8E", size = 1) +
  geom_point(color = "#158D8E") +
  theme_bw() +
  scale_y_log10() +
  labs(x= "50% similarity clusters", y = "Predator mass (mg)") +
  theme(axis.title.y= element_text(size = 20),
        axis.title.x= element_text(size = 20),
        axis.text = element_text(size = 15))

# One DF ---------------------------------------------------------------

#get one DF as opposed to list
big_stages <- bind_rows(cluster_dfs)
#how many clustered:
big_stages %>%
  group_by(sample_str) %>%
  summarise(n = n())
#HEV: 44/53 
#NEO: 14/24
#PHH: 38/42

# Visualize ---------------------------------------------------------------
pred_cols <- c("HEV" = "#EEB00C", "NEO" = "#158D8E", "PHH" = "#114C54")

rich_HEV <- big_stages %>%
  filter(sample_str == "HEV") %>%
  left_join(prey_richness, by = c("sample", "sample_str")) %>%
  ggplot(aes(x = sample_str, y= richness)) +
  geom_boxplot(size = 1, color = "#114C54") + 
  geom_beeswarm(size = 2, color = "#114C54") +
  labs(y = "Prey richness") +
  theme_bw() +
  scale_y_continuous(breaks = c(0, 1, 2, 3)) +
  ylim(0,3) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 20))

rich_NEO <- big_stages %>%
  filter(sample_str == "NEO") %>%
  left_join(prey_richness, by = c("sample", "sample_str")) %>%
  ggplot(aes(x = sample_str, y= richness)) +
  geom_boxplot(size = 1, color = "#EEB00C") + 
  geom_beeswarm(size = 2, color = "#EEB00C") +
  labs(y = "Prey richness") +
  theme_bw() +
  scale_y_continuous(breaks = c(0, 1, 2, 3)) +
  ylim(0,3) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 20))

rich_PHH <- big_stages %>%
  filter(sample_str == "PHH") %>%
  left_join(prey_richness, by = c("sample", "sample_str")) %>%
  ggplot(aes(x = sample_str, y= richness)) +
  geom_boxplot(size = 1, color = "#158D8E") + 
  geom_beeswarm(size = 2, color = "#158D8E") +
  labs(y = "Prey richness") +
  theme_bw() +
  scale_y_continuous(breaks = c(0, 1, 2, 3)) +
  ylim(0,3) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 20))


dendro_hev <- ggdendrogram(cluster_objs[[1]]) +
  theme_void() +
  geom_hline(yintercept = 0.52, 
             linetype = "dashed", 
             color = "#114C54", 
             size =1) +
  labs(y = "Jaccard similarity") +
  theme(axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 20),
        axis.title.y = element_text(angle = 90),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

dendro_NEO <- ggdendrogram(cluster_objs[[2]]) +
  theme_void() +
  geom_hline(yintercept = 0.52, 
             linetype = "dashed", 
             color = "#EEB00C", 
             size =1) +
  labs(y = "Jaccard similarity") +
  theme(axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 20),
        axis.title.y = element_text(angle = 90),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

dendro_PHH <- ggdendrogram(cluster_objs[[3]]) +
  theme_void() +
  geom_hline(yintercept = 0.52, 
             linetype = "dashed", 
             color = "#158D8E", 
             size =1) +
  labs(y = "Jaccard similarity") +
  theme(axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 20),
        axis.title.y = element_text(angle = 90),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
