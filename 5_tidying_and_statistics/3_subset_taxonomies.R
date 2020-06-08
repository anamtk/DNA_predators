###########################
#Subset taxonomies
#Ana Miller-ter Kuile
#June 8, 2020
###########################

#Need to subset prey taxonomies per predator from the rarefied dataset. 
#This will include giving the predators names that match the ASV ID, 
#so that I can filter out the ASVs where these match from any diet 
#analysis

###########################
#Load Packages ####
###########################
library(here)
library(tidyverse)
library(fuzzyjoin) #fuzzy inner join for string detection join

###########################
#Load Data ####
###########################
#Prey DNA taxonomies
taxa <- read.csv(here("data", "outputs", "1_taxonomic_assignment", "ASV_taxonomies.csv"))

preds <- read.csv(here("data", "Predator_IDs.csv"))
#load community data
comm <- read.csv(here("data", "outputs", "4_rarefied", "community_rare.csv"))

comm <- comm %>%
  rename("ASV" = "X.1") %>%
  dplyr::select(-X) 

###########################
#Manipulate Data to Long and Attach Predator IDs ####
###########################
#make long by sample and read abundance
comm_long <- comm %>%
  gather(sample, reads, HEV01a:HEV99d)

#this  binds to the predator ID DF, detecting the string from the sample_str
#column in the predator DF in each sample name. 
comm_long <- comm_long %>%
  fuzzy_inner_join(preds, by = c("sample" = "sample_str"), match_fun = str_detect)

###########################
#Attach prey IDs for each ASV ####
###########################
taxa_comm <- comm_long %>%
  left_join(taxa, by = "ASV")

###########################
#subset by species to sort ####
###########################

#For each predator, may need to revisit some of the high-read things that
#DON'T match to family or lower, as these are most likely predator DNA
#My thought for correcting for this is to create a master predator DNA
#DF, and then fit a distribution to these data, then use  this (ala Jerde)
#to predict if ASV reads above a certain amount should be subset as predator
#if they match to order or Class of the predator species in question

###########################
##HEV ####
###########################

#subset HEV predators
HEV <- taxa_comm %>%
  filter(sample_str == "HEV")

#ID matched to species?
#HEV %>% 
#  filter(pred_ID == unique_ID)

#ID matched to genus?
#HEV %>% 
 # filter(pred_Genus == unique_ID)

#ID matched to family?
heva <- HEV %>% 
  filter(pred_Family == unique_ID)

#ID matched on bold?
hevb <- HEV %>% 
  filter(pred_ID == ID_bold)

###########################
##NEO ####
###########################
#NEO
NEO <- taxa_comm %>%
  filter(sample_str == "NEO")

#ID matched to species?
neoa <- NEO %>% 
  filter(pred_ID == unique_ID)

#ID matched to genus?
neob <- NEO %>% 
  filter(pred_Genus == unique_ID)

#ID matched to family?
#NEO %>% 
#  filter(pred_Family == unique_ID)

#ID matched to bold?
neoc <- NEO %>% 
  filter(pred_ID == ID_bold)

###########################
##SCY ####
###########################
#SCY
SCY <- taxa_comm %>%
  filter(sample_str == "SCY")

#ID matched to species?
#SCY %>% 
#  filter(pred_ID == unique_ID)

#ID matched to genus?
#SCY %>% 
#  filter(pred_Genus == unique_ID)

#ID matched to family?
#SCY %>% 
#filter(pred_Family == unique_ID)

#ID matched to bold?
scya <- SCY %>% 
  filter(pred_ID == ID_bold)

###########################
##CEN ####
###########################
#CEN 
#subset centipede
CEN <- taxa_comm %>%
  filter(sample_str == "CEN")

#ID matched to order?
cena <- CEN %>% 
  filter(pred_ID == unique_ID)

#ID matched to bold?
#CEN %>% 
#  filter(pred_ID == ID_bold)

###########################
##PHH ####
###########################

#PHH subset 
PHH <- taxa_comm %>%
  filter(sample_str == "PHH")

#ID matched to species?
#PHH %>% 
#  filter(pred_ID == unique_ID)

#ID matched to genus?
##PHH %>% 
 #filter(pred_Genus == unique_ID)

#ID matched to family?
#PHH %>% 
# filter(pred_Family == unique_ID)

#ID matched to order? - not sure i trust this, want to be able to distinguish pred from non
phha <- PHH %>%
  filter(pred_Order == unique_ID)

#ID matched to bold order?
phhb <- PHH %>% 
  filter(pred_Order == ID_bold)

###########################
##SME ####
###########################

#SME subset 
SME <- taxa_comm %>%
  filter(sample_str == "SME")

#ID matched to species?
smea <- SME %>% 
  filter(pred_ID == unique_ID)

#ID matched to genus?
#SME %>% 
# filter(pred_Genus == unique_ID)

#ID matched to family?
#SME %>% 
# filter(pred_Family == unique_ID)

#ID matched to bold?
smeb <- SME %>% 
  filter(pred_ID == ID_bold)

###########################
##EUB ####
###########################
#EUB subset 
EUB <- taxa_comm %>%
  filter(sample_str == "EUB")

#ID matched to species?
#EUB %>% 
#  filter(pred_ID == unique_ID)

#ID matched to genus?
##EUB %>% 
# filter(pred_Genus == unique_ID)

#ID matched to family?
#EUB %>% 
# filter(pred_Family == unique_ID)

#Order?
#EUB %>%
#  filter(pred_Order == unique_ID)

#ID matched to bold?
#EUB %>% 
#  filter(pred_ID == ID_bold)

#Order matched bold ID?
euba <- EUB %>% 
  filter(pred_Order == ID_bold)

###########################
##LRS ####
###########################

#LRS subset 
LRS <- taxa_comm %>%
  filter(sample_str == "LRS")

#ID matched to species?
lrsa <- LRS %>% 
  filter(pred_ID == unique_ID)

#ID matched to genus?
##LRS %>% 
 #filter(pred_Genus == unique_ID)

#ID matched to family?
#LRS %>% 
# filter(pred_Family == unique_ID)

#ID matched to bold?
#LRS %>% 
#  filter(pred_ID == ID_bold)

###########################
##PAN ####
###########################
#PAN subset 
PAN <- taxa_comm %>%
  filter(sample_str == "PAN")

#ID matched to species?
pana <- PAN %>% 
  filter(pred_ID == unique_ID)

#ID matched to genus?
#PAN %>% 
# filter(pred_Genus == unique_ID)

#ID matched to family?
##PAN %>% 
#filter(pred_Family == unique_ID)

#Order?
#PAN %>%
#  filter(pred_Order == unique_ID)

#ID matched to bold?
panb <- PAN %>% 
  filter(pred_ID == ID_bold)

###########################
#Bind them all together and only take the unique ones ####
###########################
known_pred <- heva %>%
  bind_rows(hevb) %>%
  bind_rows(neoa) %>%
  bind_rows(neob) %>%
  bind_rows(neoc) %>%
  bind_rows(scya) %>%
  bind_rows(cena) %>%
  bind_rows(phha) %>%
  bind_rows(phhb) %>%
  bind_rows(smea) %>%
  bind_rows(smeb) %>%
  bind_rows(euba) %>%
  bind_rows(lrsa) %>%
  bind_rows(pana) %>%
  bind_rows(panb) %>%
  group_by(ASV, sample, reads, pred_ID, unique_ID) %>%
  distinct()

###########################
#Use abundances to inform unknowns####
###########################

#Based on this, it looks like there are a lot of zeros, but also a lot
#of high-read ASVs. IF we subset just thoes high read ASVs, we'll get the
#ASVs from samples that maybe didn't get close enough taxonomies to subset
#to be able to subset predator DNA that has high read and a higher-level
#taxonomic assignment that *could* match predator (or prey) DNA

pred_dist <- known_pred %>%
  filter(reads > 0)

hist(pred_dist$reads)

pred_reads <- as.data.frame(quantile(pred_dist$reads, c(.86, .87, .88, .89, .9, .91, .92, .93, .94, .95, 
                            .96, .97, .98, .99, 1)))

pred_reads <- pred_reads %>%
  rename("reads" = "quantile(pred_dist$reads, c(0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1))") %>%
  rownames_to_column(var = "quantile") %>%
  mutate(number = 1:n()) %>%
  mutate(number = number+85)

ggplot(pred_reads, aes(x = number, y = reads)) +
  geom_point() +theme_bw() +
  geom_hline(yintercept = 6158, linetype = "dashed") +
  labs(y = "Sequence reads", x = "Quantile")

#According to this, I'd be pretty confident that anything above 6158 sequence
#reads AND that matched a higher-order taxonomy of the predator species 
#in question, that it would likely ALSO be predator DNA. 

###########################
#Revisit predator-subset DF to remove these high abundances ####
###########################

#none in this DF above our cutoff
HEV_prey <- HEV %>%
  anti_join(heva, by = c("sample", "ASV")) %>%
  anti_join(hevb, by = c("sample", "ASV"))

#here, too, nothing above our cutoff
PAN_prey <- PAN %>%
  anti_join(pana, by = c("sample", "ASV")) %>%
  anti_join(panb, by = c("sample", "ASV"))

#some still in here 
CEN_prey <- CEN %>%
  anti_join(cena, by = c("sample", "ASV")) %>%
  filter(Class != "Chilopoda") %>% #some chilopoda were not removed above
  filter(reads < 6158) #some "arthropod" DNA had 10,000+ reads, probs predator

#save that additional predator DNA for later
cenb <- CEN %>%
  anti_join(cena, by = c("sample", "ASV")) %>%
  filter(Class == "Chilopoda") %>% 
  filter(reads > 6158) 

#nothing in here
EUB_prey <- EUB %>%
  anti_join(euba, by = c("sample", "ASV"))

#still some high reads in here
LRS_prey <- LRS %>%
  anti_join(lrsa, by = c("sample", "ASV")) %>%
  filter(reads < 6158)

#save the additional predaot DNA for later
lrsb <- LRS %>%
  anti_join(lrsa, by = c("sample", "ASV")) %>%
  filter(reads > 6158)

#still some in here too
NEO_prey <- NEO %>%
  anti_join(neoa, by = c("sample", "ASV")) %>%
  anti_join(neob, by = c("sample", "ASV")) %>%
  anti_join(neoc, by = c("sample", "ASV")) %>%
  filter(reads < 6158)

#save the additional predator DNA for later
neod <- NEO %>%
  anti_join(neoa, by = c("sample", "ASV")) %>%
  anti_join(neob, by = c("sample", "ASV")) %>%
  anti_join(neoc, by = c("sample", "ASV")) %>%
  filter(reads > 6158)

#none in here
PHH_prey <- PHH %>%
  anti_join(phha, by = c("sample", "ASV")) %>%
  anti_join(phhb, by = c("sample", "ASV")) 

#still some high ones in here
SCY_prey <- SCY %>%
  anti_join(scya, by = c("sample", "ASV")) %>%
  filter(reads < 6158)

#save predator DNA for later 
scyb <- SCY %>%
  anti_join(scya, by = c("sample", "ASV")) %>%
  filter(reads > 6158)

#none in here
SME_prey <- SME %>%
  anti_join(smea, by = c("sample", "ASV")) %>%
  anti_join(smeb, by = c("sample", "ASV"))

###########################
#Bind them all back together ####
###########################
#predator - add additional 
predator_DNA <- known_pred %>%
  bind_rows(cenb) %>%
  bind_rows(lrsb) %>%
  bind_rows(neod) %>%
  bind_rows(scyb)

other_DNA <- HEV_prey %>%
  bind_rows(PAN_prey) %>%
  bind_rows(CEN_prey) %>%
  bind_rows(EUB_prey) %>%
  bind_rows(LRS_prey) %>%
  bind_rows(NEO_prey) %>%
  bind_rows(PHH_prey) %>%
  bind_rows(SCY_prey) %>%
  bind_rows(SME_prey) #319352

prey_DNA <- other_DNA %>%
  filter(unique_ID != "NA") %>% #remove unassigned, total 152861
  group_by(sample, unique_ID, pred_ID, Phylum, Class, Order, Family, 
           Genus, Species, ID_level) %>% #select only variables i wnat going forward
  summarise(reads = sum(reads)) %>% #summarise reads by unique ID
  rename("Phylum_prey" = "Phylum",
         "Class_prey" = "Class",
         "Order_prey" = "Order",
         "Family_prey" = "Family",
         "Genus_prey" = "Genus",
         "Species_prey" = "Species") #rename taxonomy columns so they clearly correspond to prey (not predator)

write.csv(predator_DNA, here("data", "outputs", "5_rarefied_taxonomic_sort", "predator_DNA.csv"))

write.csv(prey_DNA, here("data", "outputs", "5_rarefied_taxonomic_sort", "prey_DNA.csv"))

write.csv(other_DNA, here("data", "outputs", "5_rarefied_taxonomic_sort", "non-predator_DNA.csv"))


