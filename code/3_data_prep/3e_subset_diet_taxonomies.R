###########################
#Subset taxonomies
#Ana Miller-ter Kuile
#June 8, 2020
###########################

#this code subsets prey taxonomies per predator from the rarefied dataset. 


###########################
#Load Packages ####
###########################
package.list <- c("here", "tidyverse", "fuzzyjoin")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

###########################
#Load Data ####
###########################
#Prey DNA taxonomies
taxa <- read.csv(here("data", 
                      "outputs", 
                      "3b_taxonomic_assignment", 
                      "c_final_dataset",
                      "ASV_taxonomies_summed_wIndiv.csv"))
taxa <- taxa %>%
  dplyr::select(-X)

preds <- read.csv(here("data",
                       "raw_data",
                       "Predator_IDs.csv"))
#load community data
comm <- read.csv(here("data", 
                      "outputs", 
                      "3c_rarefied", 
                      "community_rare.csv"))

comm <- comm %>%
  rename("ASV" = "X")

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
  left_join(taxa, by = "ASV") %>%
  filter(Domain == "Eukaryota") %>% #remove NA taxonomies from this community
  mutate(run = substr(sample, nchar(sample)-1+1, nchar(sample))) #indicates the 
#run it was run on

taxa_comm %>%
  filter(reads > 0) %>%
  group_by(ASV) %>%
  tally()
comm_long %>%
  filter(reads > 0) %>%
  group_by(ASV) %>%
  tally() %>%
  tally()

comm_long %>%
  filter(reads > 0) %>%
  group_by(ASV) %>%
  summarise(reads = sum(reads)) %>%
  left_join(taxa, by = "ASV") %>%
  filter(!is.na(ID_level)) %>%
  tally()
###########################
#subset by species to sort ####
###########################

#a: HEVa, NEO, SCY
#b: CEN, PHH, SME, HEVb
#c: EUB, ISO, LRS, PAN, HEVc
#d: HEVd

###########################
##HEV ####
###########################

#subset HEV predators
HEV <- taxa_comm %>%
  filter(sample_str == "HEV")

#subset predator DNA from this:
HEV_pred <- HEV %>%
  ungroup() %>%
  filter(Order == "Araneae") %>% #only order for spiders
  filter(Family %in% c("Sparassidae", "")) #either matched to the family or no family (being conservative)

#a: HEVa, NEO, SCY
HEVa_other_pred <- HEV %>%
  filter(run == "a") %>%
  ungroup() %>%
  filter(Order == "Araneae") %>%
  filter(Family %in% c("Araneidae", "Scytodidae"))


###########################
##NEO ####
###########################
#NEO
NEO <- taxa_comm %>%
  filter(sample_str == "NEO")

NEO_pred <- NEO %>%
  ungroup() %>%
  filter(Order == "Araneae") %>% #only order for spiders
  filter(Family %in% c("Araneidae", "", "Theridiidae")) #either matched to the family or no family (being conservative)

#a: HEVa, NEO, SCY
NEO_other_pred <- NEO %>%
  ungroup() %>%
  filter(Order == "Araneae") %>% #only order for spiders
  filter(Family %in% c("Scytodidae", "Sparassidae")) #either matched to the family or no family (being conservative)

###########################
##SCY ####
###########################
#SCY
SCY <- taxa_comm %>%
  filter(sample_str == "SCY")

SCY_pred <- SCY %>%
  ungroup() %>%
  filter(Order == "Araneae") %>% #only order for spiders
  filter(Family %in% c("Scytodidae", "")) #either matched to the family or no family (being conservative)

#a: HEVa, NEO, SCY
SCY_other_pred <- SCY %>%
  ungroup() %>%
  filter(Order == "Araneae") %>% #only order for spiders
  filter(Family %in% c("Sparassidae", "Araneidae")) #either matched to the family or no family (being conservative)

###########################
##CEN ####
###########################
#CEN 
#subset centipede
CEN <- taxa_comm %>%
  filter(sample_str == "CEN")

CEN_pred <- CEN %>%
  ungroup() %>%
  filter(Class == "Chilopoda") #only class of centipedes

#b: CEN, PHH, SME, HEVb
CEN_other_pred <- CEN %>%
  ungroup() %>%
  filter(Order %in% c("Orthoptera", "Araneae")) %>%
  filter(Family %in% c("Sparassidae", "Tettigoniidae", "Pholcidae"))


###########################
##PHH ####
###########################

#PHH subset 
PHH <- taxa_comm %>%
  filter(sample_str == "PHH")

PHH_pred <- PHH %>%
  ungroup() %>%
  filter(Order == "Orthoptera") %>% #only order for katydids
  filter(Family %in% c("Tettigoniidae", "")) #either matched to the family or no family (being conservative)

#b: CEN, PHH, SME, HEVb
PHH_other_pred <- PHH %>%
  ungroup() %>%
  filter(Order %in% c("Geophilomorpha", "Araneae")) %>%
  filter(Family %in% c("", "Sparassidae", "Pholcidae"))

###########################
##SME ####
###########################

#SME subset 
SME <- taxa_comm %>%
  filter(sample_str == "SME")

SME_pred <- SME %>%
  ungroup() %>%
  filter(Order == "Araneae") %>% #only order for spiders
  filter(Family %in% c("Pholcidae", "")) #either matched to the family or no family (being conservative)

#b: CEN, PHH, SME, HEVb
SME_other_pred <- SME %>%
  ungroup() %>%
  filter(Order %in% c("Geophilomorpha", "Araneae")) %>%
  filter(Family %in% c("", "Sparassidae", "Tettigoniidae"))

###########################
##EUB ####
###########################
#EUB subset 
EUB <- taxa_comm %>%
  filter(sample_str == "EUB")

EUB_pred <- EUB %>%
  ungroup() %>%
  filter(Order == "Dermaptera") %>% #only order for earwigs
  filter(Family %in% c("Anisolabididae", "")) #either matched to the family or no family (being conservative)

#c: EUB, ISO, LRS, PAN, HEVc
EUB_other_pred <- EUB %>%
  ungroup() %>%
  filter(Order %in% c("Odonata", "Scorpiones", "Araneae")) %>%
  filter(Family %in% c("Libellulidae", "Buthidae", "Theridiidae", "Sparassidae"))

###########################
##LRS ####
###########################

#LRS subset 
LRS <- taxa_comm %>%
  filter(sample_str == "LRS")

LRS_pred <- LRS %>%
  ungroup() %>%
  filter(Order == "Araneae") %>% #only order for spiders
  filter(Family %in% c("Oonopidae", "")) 
#c: EUB, ISO, LRS, PAN, HEVc
LRS_other_pred <- LRS %>%
  ungroup() %>%
  filter(Order %in% c("Odonata", "Scorpiones", "Dermaptera", "Araneae")) %>%
  filter(Family %in% c("Libellulidae", "Buthidae", "Anisolabididae", "Sparassidae"))

###########################
##PAN ####
###########################
#PAN subset 
PAN <- taxa_comm %>%
  filter(sample_str == "PAN")

PAN_pred <- PAN %>%
  ungroup() %>%
  filter(Order == "Odonata") %>% #only order for dragonflies
  filter(Family %in% c("Libellulidae", "")) #either matched to the family or no family (being conservative)

#c: EUB, ISO, LRS, PAN, HEVc
PAN_other_pred <- PAN %>%
  ungroup() %>%
  filter(Order %in% c("Scorpiones", "Dermaptera", "Araneae")) %>%
  filter(Family %in% c("Theridiidae", "Buthidae", "Anisolabididae", "Sparassidae"))

###########################
#Bind them all together and only take the unique ones ####
###########################
known_pred <- CEN_pred%>%
  bind_rows(EUB_pred) %>%
  bind_rows(HEV_pred) %>%
  bind_rows(LRS_pred) %>%
  bind_rows(NEO_pred) %>%
  bind_rows(PAN_pred) %>%
  bind_rows(PHH_pred) %>%
  bind_rows(SCY_pred) %>%
  bind_rows(SME_pred)

other_pred <- CEN_other_pred %>%
  bind_rows(EUB_other_pred) %>%
  bind_rows(HEVa_other_pred) %>%
  bind_rows(LRS_other_pred) %>%
  bind_rows(NEO_other_pred) %>%
  bind_rows(PAN_other_pred) %>%
  bind_rows(PHH_other_pred) %>%
  bind_rows(SCY_other_pred) %>%
  bind_rows(SME_other_pred) %>%
  filter(reads > 0)
  

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
  geom_hline(yintercept = 12098, linetype = "dashed") +
  geom_hline(yintercept = 3554, linetype = "dashed") +
  labs(y = "Sequence reads", x = "Quantile")

#According to this, I'd be pretty confident that anything above 12098 sequence
#reads AND that matched a higher-order taxonomy of the predator species 
#in question, that it would likely ALSO be predator DNA. OR
#anything *below 3554 that matched to the predator higher taxonomy is 
#likely prey.
#from looking at this, though - it doesn't seem to be that relevant
#as we will be deleting those higher level taxonomies from our analysis anyway

###########################
#Bind them all back together ####
###########################
other_DNA <- taxa_comm %>%
  anti_join(known_pred, by = c("sample", "ASV"))

other_DNA_conservative <- taxa_comm %>%
  anti_join(known_pred, by = c("sample", "ASV")) %>%
  anti_join(other_pred, by = c("sample", "ASV"))

other_DNA_conservative %>%
  filter(reads >0) %>%
  group_by(ASV) %>%
  tally() %>%
  tally()

prey_fam <- other_DNA %>%
  filter(Family != "")

prey_fam_conservative <- other_DNA_conservative %>%
  filter(Family != "")

prey_fam_conservative %>%
  filter(reads >0) %>%
  group_by(ASV) %>%
  tally() %>%
  tally()

###########################
#Export ####
###########################

write.csv(known_pred, here("data", 
                           "outputs", 
                           "3d_rarefied_taxonomic_sort", 
                           "predator_DNA.csv"))

write.csv(other_DNA, here("data", 
                          "outputs",
                          "3d_rarefied_taxonomic_sort", 
                          "all_prey_DNA.csv"))

write.csv(prey_fam_conservative, here("data", 
                                      "outputs", 
                                      "3d_rarefied_taxonomic_sort", 
                                      "fam_prey_DNA_conservative.csv"))

