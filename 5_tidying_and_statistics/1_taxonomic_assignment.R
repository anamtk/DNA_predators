#####################
#Taxonomic Assignments####
#Ana Miller -ter Kuile
#May 29, 2020
#####################

#This script compiles taxonomic assignments for all Palmyra predators
#from both NCBI and BOLD database searches. It will create one
#taxonomic ID column based on combinations of the best taxonomic
#assignment between both databases
#the output should include that taxonomy, and then data on 
#species, genus, family, and order from either GenBank or BOLD

#####################
#Load Packages ####
#####################
library(here) #easy file paths
library(tidyverse) #tidy data

#####################
#Import and tidy NCBI Data####
#####################

#NCBI IDs
IDs_ncbi <- read.csv(here("2_taxonomic_assignment", "taxonomies", "ncbi.csv"))

#set all columns that are factors as characters so we can use 
#the if-else statement below
IDs_ncbi <- IDs_ncbi %>%
  mutate_if(is.factor, as.character)

#set an ID column based on the lowest taxonomic assignment
IDs_ncbi$ID_ncbi <- ifelse(IDs_ncbi$Species != "", IDs_ncbi$Species,
                           ifelse(IDs_ncbi$Genus != "", IDs_ncbi$Genus,
                                  ifelse(IDs_ncbi$Family != "", IDs_ncbi$Family,
                                         ifelse(IDs_ncbi$Order != "", IDs_ncbi$Order,
                                                ifelse(IDs_ncbi$Class != "", IDs_ncbi$Class,
                                                       ifelse(IDs_ncbi$Phylum != "", IDs_ncbi$Phylum, IDs_ncbi$Domain))))))

#set a level at which the ID has been made
IDs_ncbi$ID_level_ncbi <- ifelse(IDs_ncbi$Species != "", "Species",
                            ifelse(IDs_ncbi$Genus != "", "Genus",
                                   ifelse(IDs_ncbi$Family != "", "Family",
                                          ifelse(IDs_ncbi$Order != "", "Order",
                                                 ifelse(IDs_ncbi$Class != "", "Class",
                                                        ifelse(IDs_ncbi$Phylum != "", "Phylum", "Domain"))))))

#remove the weird thing that MEGAN puts in front of the ID in that column
#so that it is just the ID from that column
IDs_ncbi$ID_ncbi <- sapply(str_split(IDs_ncbi$ID_ncbi, "_"), function(x){return(x[[3]])})

#before we bind to bold, I want to ensure I keep track of which taxonomies
#come from where, so I'm going to rename the taxonomy columns to indicate
#which DB they came from
IDs_ncbi <- IDs_ncbi %>%
  rename(Domain_ncbi = Domain,
         Phylum_ncbi = Phylum,
         Class_ncbi = Class,
         Order_ncbi = Order,
         Family_ncbi = Family,
         Genus_ncbi = Genus,
         Species_ncbi = Species)

#####################
#Import and tidy BOLD Data####
#####################

#BOLD IDs from the IDEngine in batches
bold1 <- read.csv(here("2_taxonomic_assignment", "taxonomies", "BOLD_0.csv"))
bold2 <- read.csv(here("2_taxonomic_assignment", "taxonomies", "BOLD_1.csv"))
bold3 <- read.csv(here("2_taxonomic_assignment", "taxonomies", "BOLD_2.csv"))
bold4 <- read.csv(here("2_taxonomic_assignment", "taxonomies", "BOLD_3.csv"))
bold5 <- read.csv(here("2_taxonomic_assignment", "taxonomies", "BOLD_4.csv"))
bold6 <- read.csv(here("2_taxonomic_assignment", "taxonomies", "BOLD_5.csv"))
bold7 <- read.csv(here("2_taxonomic_assignment", "taxonomies", "BOLD_6.csv"))
bold8 <- read.csv(here("2_taxonomic_assignment", "taxonomies", "BOLD_7.csv"))
bold9 <- read.csv(here("2_taxonomic_assignment", "taxonomies", "BOLD_8.csv"))
bold10 <- read.csv(here("2_taxonomic_assignment", "taxonomies", "BOLD_9.csv"))
bold11 <- read.csv(here("2_taxonomic_assignment", "taxonomies", "BOLD_10.csv"))
bold12 <- read.csv(here("2_taxonomic_assignment", "taxonomies", "BOLD_11.csv"))
bold13 <- read.csv(here("2_taxonomic_assignment", "taxonomies", "BOLD_12.csv"))
bold14 <- read.csv(here("2_taxonomic_assignment", "taxonomies", "BOLD_13.csv"))
bold15 <- read.csv(here("2_taxonomic_assignment", "taxonomies", "BOLD_14.csv"))
bold16 <- read.csv(here("2_taxonomic_assignment", "taxonomies", "BOLD_15.csv"))
bold17 <- read.csv(here("2_taxonomic_assignment", "taxonomies", "BOLD_16.csv"))
bold18 <- read.csv(here("2_taxonomic_assignment", "taxonomies", "BOLD_17.csv"))

#bind these all together
all_bold <- bold1 %>%
  bind_rows(bold2) %>%
  bind_rows(bold3) %>%
  bind_rows(bold4) %>%
  bind_rows(bold5) %>%
  bind_rows(bold6) %>%
  bind_rows(bold7) %>%
  bind_rows(bold8) %>%
  bind_rows(bold9) %>%
  bind_rows(bold10) %>%
  bind_rows(bold11) %>%
  bind_rows(bold12) %>%
  bind_rows(bold13) %>%
  bind_rows(bold14) %>%
  bind_rows(bold15) %>%
  bind_rows(bold16) %>%
  bind_rows(bold17) %>%
  bind_rows(bold18) %>%
  filter(!Best.ID == "No match") %>%
  rename("ID_bold" = "Best.ID", "ASV" = "Query.ID") 

#I exported these and assigned taxonomic levels to everything via internet searches
#export so I can add taxonomic levels to this database
write.csv(all_bold, here("2_taxonomic_assignment", "taxonomies", "bold.csv"))

#I created a new CSV that includes all the taxonomic levels I compiled from internet
#searches
IDs_bold <- read.csv(here("2_taxonomic_assignment", "taxonomies", "bold_wID.csv"))

#now I want to assign the level at which those IDs were made, as in the above DF
#from NCBI
IDs_bold <- IDs_bold %>%
  mutate_if(is.factor, as.character)

#set a level at which the ID has been made
IDs_bold$ID_level_bold <- ifelse(IDs_bold$Species != "", "Species",
                                 ifelse(IDs_bold$Genus != "", "Genus",
                                        ifelse(IDs_bold$Family != "", "Family",
                                               ifelse(IDs_bold$Order != "", "Order",
                                                      ifelse(IDs_bold$Class != "", "Class",
                                                             ifelse(IDs_bold$Phylum != "", "Phylum", "Domain"))))))


#Again,  I want to keep track of where the taxonomies are coming from, so I
#am going to rename columns
IDs_bold <- IDs_bold %>%
  dplyr::select(-Search.DB, -Top.., -Low.., -X, -X.1) %>% #remove weird rows from BOLD
  rename(Domain_bold = Domain,
         Phylum_bold = Phylum, #rename columns
         Class_bold = Class,
         Order_bold = Order,
         Family_bold = Family,
         Genus_bold = Genus,
         Species_bold = Species)

#####################
#Combine datasets of IDs ####
#####################
  
#now we can join these into one DF with both bold and NCBI identifications and
#do some initial cleaning
IDs <- IDs_ncbi %>%
  full_join(IDs_bold, by = "ASV") %>% #join them both
  filter(!Type %in% c("unclear", "non-diet")) %>% #remove things that got multiple taxonomic
#matches from BOLD
  replace(., is.na(.), "") #replace NA values with nothing for ease of sorting later
#####################
#Make a combined ID from best knowledge of both databases####
#####################

#This dataset now includes four types of IDs:
#1. Those which matched to NCBI and not BOLD
#2. Those which matched to BOLD but not NCBI
#3. Those which matched to both and are the same
#4. Those which matched to both and are different

#I am going to break the dataframe into these four component parts (for now)
#just because that makes it easier for me to see what the composition of
#these four types are

#NCBI only
IDs_ncbi_only <- IDs %>%
  filter(ID_ncbi != "" & ID_bold == "")

#get stats on this:
IDs_ncbi_only %>%
  tally()

IDs_ncbi_only %>%
  group_by(ID_level_ncbi) %>%
  tally()

#BOLD only
IDs_bold_only <- IDs %>%
  filter(ID_bold != "" & ID_ncbi == "")

#get stats on this:
IDs_bold_only %>%
  tally()

IDs_bold_only %>%
  group_by(ID_level_bold) %>%
  tally()

#Both Databases
IDs_both <- IDs %>%
  filter(ID_bold != "", ID_ncbi != "")

IDs_both %>%
  tally()

IDs_both %>%
  group_by(ID_level_bold, ID_level_ncbi) %>%
  tally()

IDs_same_level <- IDs %>%
  filter(ID_level_ncbi == ID_level_bold) %>%
  dplyr::select(ID_ncbi, ID_bold)

View(IDs_species)
