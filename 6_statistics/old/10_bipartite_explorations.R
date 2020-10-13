#################
#bipartite explorations
################

###########################
# Load packages
package.list <- c("here", "tidyverse", "ggplot2", "bipartite")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}
#############################

#############################

pal <- read.csv(here("data", "outputs",
                     "5_rarefied_taxonomic_sort",
                     "fam_prey_DNA.csv"))
  
pal$sample <- str_sub(pal$sample, end=-2)

meta <- read.csv(here("data",
                      "Sample_metadata.csv"))

#############################

pal1 <- pal %>%
  filter(reads > 0) %>%
  group_by(sample, Family) %>%
  summarise(reads = sum(reads)) %>%
  mutate(presence = ifelse(reads > 0, 1, 0)) %>%
  dplyr::select(-reads) %>%
  pivot_wider(names_from = "sample",
              values_from = "presence") %>%
  filter(!Family %in% c("Hominidae", "Muridae", "Salmonidae", "Sulidae"))

pal1[is.na(pal1)] <- 0

pal1 <- as.data.frame(pal1)

rownames(pal1) <- pal1$Family

pal1 <- pal1[,-1]

plotweb(pal1)
visweb(pal1)

plotPAC(PAC(pal1), outby=0.9)
mod <- computeModules(pal1)
plotModuleWeb(mod)

###############################
#HEV####
##########################

meta <- meta %>%
  group_by(Extraction.ID) %>%
  summarise(length_mm = mean(Length_mm))

hev_sizes <- pal %>%
  filter(reads > 0) %>%
  filter(pred_Family == "Sparassidae") %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  group_by(sample, length_mm) %>%
  summarise(reads) %>%
  dplyr::select(-reads) %>%
  distinct(sample, length_mm) %>%
  arrange(length_mm)

samples <- hev_sizes$sample

hev <- pal %>%
  filter(reads > 0) %>%
  filter(pred_Family == "Sparassidae") %>%
  group_by(sample, Family) %>%
  summarise(reads = sum(reads)) %>%
  mutate(presence = ifelse(reads > 0, 1, 0)) %>%
  dplyr::select(-reads) %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  arrange(length_mm) %>%
  mutate(length_mm = factor(length_mm)) %>%
  ungroup() %>%
  dplyr::select(-sample) %>%
  group_by(length_mm, Family) %>%
  summarise(presence = sum(presence)) %>%
  pivot_wider(names_from = "length_mm",
              values_from = "presence") %>%
  filter(!Family %in% c("Hominidae", "Muridae", "Salmonidae", "Sulidae"))

hev[is.na(hev)] <- 0

hev <- as.data.frame(hev)

rownames(hev) <- hev$Family

hev <- hev[,-1]

plotweb(hev)
str(hev)

visweb(hev)

mod <- computeModules(hev)
plotModuleWeb(mod)

hev1 <- pal %>%
  filter(pred_Family == "Sparassidae") %>%
  group_by(sample, Family) %>%
  summarise(reads = sum(reads)) %>%
  mutate(presence = ifelse(reads > 0, 1, 0)) %>%
  dplyr::select(-reads) %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  arrange(length_mm) %>%
  filter(!Family %in% c("Hominidae", "Muridae", "Salmonidae", "Sulidae"))

hev_fam <- hev1 %>%
  group_by(Family) %>%
  summarise(sum = sum(presence)) %>%
  filter(sum > 0)

hev_vec <- hev_fam$Family

hev1 <- hev1 %>%
  filter(Family %in% hev_vec)

hev1$length_mm <- as.factor(hev1$length_mm)
hev1$presence <- as.factor(hev1$presence)

ggplot(hev1, aes(x = length_mm, y = Family, fill = presence)) +
  geom_tile(color = "black") + 
  coord_equal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("white", "black"))

###############################
#NEO####
##########################

meta <- meta %>%
  group_by(Extraction.ID) %>%
  summarise(length_mm = mean(Length_mm))

neo_sizes <- pal %>%
  filter(reads > 0) %>%
  filter(pred_Species == "Neoscona theisi") %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  group_by(sample, length_mm) %>%
  summarise(reads) %>%
  dplyr::select(-reads) %>%
  distinct(sample, length_mm) %>%
  arrange(length_mm)

neo <- pal %>%
  filter(reads > 0) %>%
  filter(pred_Species == "Neoscona theisi") %>%
  group_by(sample, Family) %>%
  summarise(reads = sum(reads)) %>%
  mutate(presence = ifelse(reads > 0, 1, 0)) %>%
  dplyr::select(-reads) %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  arrange(length_mm) %>%
  mutate(length_mm = factor(length_mm)) %>%
  ungroup() %>%
  dplyr::select(-sample) %>%
  group_by(length_mm, Family) %>%
  summarise(presence = sum(presence)) %>%
  pivot_wider(names_from = "length_mm",
              values_from = "presence") %>%
  filter(!Family %in% c("Hominidae", "Muridae", "Salmonidae", "Sulidae", "Sparassidae"))

neo[is.na(neo)] <- 0

neo <- as.data.frame(neo)

rownames(neo) <- neo$Family

neo <- neo[,-1]

plotweb(neo)
str(neo)

visweb(neo)

mod <- computeModules(neo)
plotModuleWeb(mod)

neo1 <- pal %>%
  filter(pred_Species == "Neoscona theisi") %>%
  group_by(sample, Family) %>%
  summarise(reads = sum(reads)) %>%
  mutate(presence = ifelse(reads > 0, 1, 0)) %>%
  dplyr::select(-reads) %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  arrange(length_mm) %>%
  filter(!Family %in% c("Hominidae", "Muridae", "Salmonidae", "Sulidae", "Sparassidae"))

neo_fam <- neo1 %>%
  group_by(Family) %>%
  summarise(sum = sum(presence)) %>%
  filter(sum > 0)

neo_vec <- neo_fam$Family

neo1 <- neo1 %>%
  filter(Family %in% neo_vec)

neo1$length_mm <- as.factor(neo1$length_mm)
neo1$presence <- as.factor(neo1$presence)

ggplot(neo1, aes(x = length_mm, y = Family, fill = presence)) +
  geom_tile(color = "black") + 
  coord_equal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("white", "black"))

###############################
#PHH####
##########################

phh_sizes <- pal %>%
  filter(reads > 0) %>%
  filter(pred_Species == "Phisis holdhausi") %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  group_by(sample, length_mm) %>%
  summarise(reads) %>%
  dplyr::select(-reads) %>%
  distinct(sample, length_mm) %>%
  arrange(length_mm)

phh <- pal %>%
  filter(reads > 0) %>%
  filter(pred_Species == "Phisis holdhausi") %>%
  group_by(sample, Family) %>%
  summarise(reads = sum(reads)) %>%
  mutate(presence = ifelse(reads > 0, 1, 0)) %>%
  dplyr::select(-reads) %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  arrange(length_mm) %>%
  mutate(length_mm = factor(length_mm)) %>%
  ungroup() %>%
  dplyr::select(-sample) %>%
  group_by(length_mm, Family) %>%
  summarise(presence = sum(presence)) %>%
  pivot_wider(names_from = "length_mm",
              values_from = "presence") %>%
  filter(!Family %in% c("Hominidae", "Muridae", "Salmonidae", "Sulidae", "Sparassidae"))

phh[is.na(phh)] <- 0

phh <- as.data.frame(phh)

rownames(phh) <- phh$Family

phh <- phh[,-1]

plotweb(phh)


visweb(phh)

mod <- computeModules(phh)
plotModuleWeb(mod)

phh1 <- pal %>%
  filter(pred_Species == "Phisis holdhausi") %>%
  group_by(sample, Family) %>%
  summarise(reads = sum(reads)) %>%
  mutate(presence = ifelse(reads > 0, 1, 0)) %>%
  dplyr::select(-reads) %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  arrange(length_mm) %>%
  filter(!Family %in% c("Hominidae", "Muridae", "Salmonidae", "Sulidae", "Sparassidae"))

phh_fam <- phh1 %>%
  group_by(Family) %>%
  summarise(sum = sum(presence)) %>%
  filter(sum > 0)

phh_vec <- phh_fam$Family

phh1 <- phh1 %>%
  filter(Family %in% phh_vec)

phh1$length_mm <- as.factor(phh1$length_mm)
phh1$presence <- as.factor(phh1$presence)

ggplot(phh1, aes(x = length_mm, y = Family, fill = presence)) +
  geom_tile(color = "black") + 
  coord_equal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("white", "black"))

###############################
#SCY####
##########################

scy_sizes <- pal %>%
  filter(reads > 0) %>%
  filter(pred_Species == "Scytodes longipes") %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  group_by(sample, length_mm) %>%
  summarise(reads) %>%
  dplyr::select(-reads) %>%
  distinct(sample, length_mm) %>%
  arrange(length_mm)

scy <- pal %>%
  filter(reads > 0) %>%
  filter(pred_Species == "Scytodes longipes") %>%
  group_by(sample, Family) %>%
  summarise(reads = sum(reads)) %>%
  mutate(presence = ifelse(reads > 0, 1, 0)) %>%
  dplyr::select(-reads) %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  arrange(length_mm) %>%
  mutate(length_mm = factor(length_mm)) %>%
  ungroup() %>%
  dplyr::select(-sample) %>%
  group_by(length_mm, Family) %>%
  summarise(presence = sum(presence)) %>%
  pivot_wider(names_from = "length_mm",
              values_from = "presence") %>%
  filter(!Family %in% c("Hominidae", "Muridae", "Salmonidae", "Sulidae", "Sparassidae"))

scy[is.na(scy)] <- 0

scy <- as.data.frame(scy)

rownames(scy) <- scy$Family

scy <- scy[,-1]

plotweb(scy)


visweb(scy)

mod <- computeModules(scy)
plotModuleWeb(mod)

scy1 <- pal %>%
  filter(pred_Species == "Scytodes longipes") %>%
  group_by(sample, Family) %>%
  summarise(reads = sum(reads)) %>%
  mutate(presence = ifelse(reads > 0, 1, 0)) %>%
  dplyr::select(-reads) %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  arrange(length_mm) %>%
  filter(!Family %in% c("Hominidae", "Muridae", "Salmonidae", "Sulidae", "Sparassidae"))

scy_fam <- scy1 %>%
  group_by(Family) %>%
  summarise(sum = sum(presence)) %>%
  filter(sum > 0)

scy_vec <- scy_fam$Family

scy1 <- scy1 %>%
  filter(Family %in% scy_vec)

scy1$length_mm <- as.factor(scy1$length_mm)
scy1$presence <- as.factor(scy1$presence)

ggplot(scy1, aes(x = length_mm, y = Family, fill = presence)) +
  geom_tile(color = "black") + 
  coord_equal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("white", "black"))


###############################
#PAN####
##########################

pan_sizes <- pal %>%
  filter(reads > 0) %>%
  filter(pred_ID == "Pantala flavescens") %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  group_by(sample, length_mm) %>%
  summarise(reads) %>%
  dplyr::select(-reads) %>%
  distinct(sample, length_mm) %>%
  arrange(length_mm)

pan <- pal %>%
  filter(reads > 0) %>%
  filter(pred_ID == "Pantala flavescens") %>%
  group_by(sample, Family) %>%
  summarise(reads = sum(reads)) %>%
  mutate(presence = ifelse(reads > 0, 1, 0)) %>%
  dplyr::select(-reads) %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  arrange(length_mm) %>%
  mutate(length_mm = factor(length_mm)) %>%
  ungroup() %>%
  dplyr::select(-sample) %>%
  group_by(length_mm, Family) %>%
  summarise(presence = sum(presence)) %>%
  pivot_wider(names_from = "length_mm",
              values_from = "presence") %>%
  filter(!Family %in% c("Hominidae", "Muridae", "Salmonidae", "Sulidae", "Sparassidae"))

pan[is.na(pan)] <- 0

pan <- as.data.frame(pan)

rownames(pan) <- pan$Family

pan <- pan[,-1]

plotweb(pan)

visweb(pan)

mod <- computeModules(pan)
plotModuleWeb(mod)

pan1 <- pal %>%
  filter(pred_ID == "Pantala flavescens") %>%
  group_by(sample, Family) %>%
  summarise(reads = sum(reads)) %>%
  mutate(presence = ifelse(reads > 0, 1, 0)) %>%
  dplyr::select(-reads) %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  arrange(length_mm) %>%
  filter(!Family %in% c("Hominidae", "Muridae", "Salmonidae", "Sulidae", "Sparassidae"))

pan_fam <- pan1 %>%
  group_by(Family) %>%
  summarise(sum = sum(presence)) %>%
  filter(sum > 0)

pan_vec <- pan_fam$Family

pan1 <- pan1 %>%
  filter(Family %in% pan_vec)

pan1$length_mm <- as.factor(pan1$length_mm)
pan1$presence <- as.factor(pan1$presence)

ggplot(pan1, aes(x = length_mm, y = Family, fill = presence)) +
  geom_tile(color = "black") + 
  coord_equal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("white", "black"))


###############################
#SME####
##########################

sme_sizes <- pal %>%
  filter(reads > 0) %>%
  filter(pred_ID == "Smeringopus pallidus") %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  group_by(sample, length_mm) %>%
  summarise(reads) %>%
  dplyr::select(-reads) %>%
  distinct(sample, length_mm) %>%
  arrange(length_mm)

sme <- pal %>%
  filter(reads > 0) %>%
  filter(pred_ID == "Smeringopus pallidus") %>%
  group_by(sample, Family) %>%
  summarise(reads = sum(reads)) %>%
  mutate(presence = ifelse(reads > 0, 1, 0)) %>%
  dplyr::select(-reads) %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  arrange(length_mm) %>%
  mutate(length_mm = factor(length_mm)) %>%
  ungroup() %>%
  dplyr::select(-sample) %>%
  group_by(length_mm, Family) %>%
  summarise(presence = sum(presence)) %>%
  pivot_wider(names_from = "length_mm",
              values_from = "presence") %>%
  filter(!Family %in% c("Hominidae", "Muridae", "Salmonidae", "Sulidae", "Sparassidae"))

sme[is.na(sme)] <- 0

sme <- as.data.frame(sme)

rownames(sme) <- sme$Family

sme <- sme[,-1]

plotweb(sme)

visweb(sme)

#Not part of analysis because of low data ####
###############################
#LRS####
##########################

lrs_sizes <- pal %>%
  filter(reads > 0) %>%
  filter(pred_ID == "Keijia mneon") %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  group_by(sample, length_mm) %>%
  summarise(reads) %>%
  dplyr::select(-reads) %>%
  distinct(sample, length_mm) %>%
  arrange(length_mm)

lrs <- pal %>%
  filter(reads > 0) %>%
  filter(pred_ID == "Keijia mneon") %>%
  group_by(sample, Family) %>%
  summarise(reads = sum(reads)) %>%
  mutate(presence = ifelse(reads > 0, 1, 0)) %>%
  dplyr::select(-reads) %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  arrange(length_mm) %>%
  mutate(length_mm = factor(length_mm)) %>%
  ungroup() %>%
  dplyr::select(-sample) %>%
  group_by(length_mm, Family) %>%
  summarise(presence = sum(presence)) %>%
  pivot_wider(names_from = "length_mm",
              values_from = "presence") %>%
  filter(!Family %in% c("Hominidae", "Muridae", "Salmonidae", "Sulidae", "Sparassidae"))

lrs[is.na(lrs)] <- 0

lrs <- as.data.frame(lrs)

rownames(lrs) <- lrs$Family

lrs <- lrs[,-1]

plotweb(lrs)

visweb(lrs)

###############################
#EUB####
##########################

eub_sizes <- pal %>%
  filter(reads > 0) %>%
  filter(pred_ID == "Euborellia annulipes") %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  group_by(sample, length_mm) %>%
  summarise(reads) %>%
  dplyr::select(-reads) %>%
  distinct(sample, length_mm) %>%
  arrange(length_mm)

eub <- pal %>%
  filter(reads > 0) %>%
  filter(pred_ID == "Euborellia annulipes") %>%
  group_by(sample, Family) %>%
  summarise(reads = sum(reads)) %>%
  mutate(presence = ifelse(reads > 0, 1, 0)) %>%
  dplyr::select(-reads) %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  arrange(length_mm) %>%
  mutate(length_mm = factor(length_mm)) %>%
  ungroup() %>%
  dplyr::select(-sample) %>%
  group_by(length_mm, Family) %>%
  summarise(presence = sum(presence)) %>%
  pivot_wider(names_from = "length_mm",
              values_from = "presence") %>%
  filter(!Family %in% c("Hominidae", "Muridae", "Salmonidae", "Sulidae", "Sparassidae"))

eub[is.na(eub)] <- 0

eub <- as.data.frame(eub)

rownames(eub) <- eub$Family

eub <- eub[,-1]

plotweb(eub)

visweb(eub)

###############################
#CEN####
##########################

cen_sizes <- pal %>%
  filter(reads > 0) %>%
  filter(pred_ID == "Geophilomorpha") %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  group_by(sample, length_mm) %>%
  summarise(reads) %>%
  dplyr::select(-reads) %>%
  distinct(sample, length_mm) %>%
  arrange(length_mm)

cen <- pal %>%
  filter(reads > 0) %>%
  filter(pred_ID == "Geophilomorpha") %>%
  group_by(sample, Family) %>%
  summarise(reads = sum(reads)) %>%
  mutate(presence = ifelse(reads > 0, 1, 0)) %>%
  dplyr::select(-reads) %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  arrange(length_mm) %>%
  mutate(length_mm = factor(length_mm)) %>%
  ungroup() %>%
  dplyr::select(-sample) %>%
  group_by(length_mm, Family) %>%
  summarise(presence = sum(presence)) %>%
  pivot_wider(names_from = "length_mm",
              values_from = "presence") %>%
  filter(!Family %in% c("Hominidae", "Muridae", "Salmonidae", "Sulidae", "Sparassidae"))

cen[is.na(cen)] <- 0

cen <- as.data.frame(cen)

rownames(cen) <- cen$Family

cen <- cen[,-1]

plotweb(cen)

visweb(cen)

mod <- computeModules(cen)
plotModuleWeb(mod)

cen1 <- pal %>%
  filter(pred_ID == "Geophilomorpha") %>%
  group_by(sample, Family) %>%
  summarise(reads = sum(reads)) %>%
  mutate(presence = ifelse(reads > 0, 1, 0)) %>%
  dplyr::select(-reads) %>%
  left_join(meta, by = c("sample" = "Extraction.ID")) %>%
  arrange(length_mm) %>%
  filter(!Family %in% c("Hominidae", "Muridae", "Salmonidae", "Sulidae", 
                        "Sparassidae", "Tettigoniidae", "Pholcidae"))

cen_fam <- cen1 %>%
  group_by(Family) %>%
  summarise(sum = sum(presence)) %>%
  filter(sum > 0)

cen_vec <- cen_fam$Family

cen1 <- cen1 %>%
  filter(Family %in% cen_vec)

cen1$length_mm <- as.factor(cen1$length_mm)
cen1$presence <- as.factor(cen1$presence)

ggplot(cen1, aes(x = length_mm, y = Family, fill = presence)) +
  geom_tile(color = "black") + 
  coord_equal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("white", "black"))


