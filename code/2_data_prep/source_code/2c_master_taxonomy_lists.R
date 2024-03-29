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

#Following the initial outputs, I also BLASTed individual ASVs and took
#corresponding family-level assignemnts for those that had consistent family-level
# matches for the first 10+ genbank matches

#####################
#Load Packages ####
#####################
package.list <- c("here", "tidyverse")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

#####################
#Import and tidy NCBI Data####
#####################

#NCBI IDs
IDs_ncbi <- read.csv(here("data",
                          "raw_data",
                          "3_taxonomies",
                          "NCBI",
                          "ncbi.csv"))

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

IDs_ncbi <- IDs_ncbi %>%
  mutate(Domain = str_sub(Domain, 4, -1)) %>%
  mutate(Phylum = str_sub(Phylum, 4, -1)) %>%
  mutate(Class = str_sub(Class, 4, -1)) %>%
  mutate(Order = str_sub(Order, 4, -1)) %>%
  mutate(Family = str_sub(Family, 4, -1)) %>%
  mutate(Genus = str_sub(Genus, 4, -1)) %>%
  mutate(Species = str_sub(Species, 4, -1))

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
bold1 <- read.csv(here("data",
                       "raw_data",
                       "3_taxonomies",
                       "BOLD",
                       "BOLD_raw",
                       "BOLD_0.csv"))
bold2 <- read.csv(here("data",
                       "raw_data",
                       "3_taxonomies",
                       "BOLD",
                       "BOLD_raw", 
                       "BOLD_1.csv"))
bold3 <- read.csv(here("data",
                       "raw_data",
                       "3_taxonomies",
                       "BOLD",
                       "BOLD_raw",
                       "BOLD_2.csv"))
bold4 <- read.csv(here("data",
                       "raw_data",
                       "3_taxonomies",
                       "BOLD",
                       "BOLD_raw",
                       "BOLD_3.csv"))
bold5 <- read.csv(here("data",
                       "raw_data",
                       "3_taxonomies",
                       "BOLD",
                       "BOLD_raw",
                       "BOLD_4.csv"))
bold6 <- read.csv(here("data",
                       "raw_data",
                       "3_taxonomies",
                       "BOLD",
                       "BOLD_raw",
                       "BOLD_5.csv"))
bold7 <- read.csv(here("data",
                       "raw_data",
                       "3_taxonomies",
                       "BOLD",
                       "BOLD_raw",
                       "BOLD_6.csv"))
bold8 <- read.csv(here("data",
                       "raw_data",
                       "3_taxonomies",
                       "BOLD",
                       "BOLD_raw",
                       "BOLD_7.csv"))
bold9 <- read.csv(here("data",
                       "raw_data",
                       "3_taxonomies",
                       "BOLD",
                       "BOLD_raw",
                       "BOLD_8.csv"))
bold10 <- read.csv(here("data",
                        "raw_data",
                        "3_taxonomies",
                        "BOLD",
                        "BOLD_raw",
                        "BOLD_9.csv"))
bold11 <- read.csv(here("data",
                        "raw_data",
                        "3_taxonomies",
                        "BOLD",
                        "BOLD_raw",
                        "BOLD_10.csv"))
bold12 <- read.csv(here("data",
                        "raw_data",
                        "3_taxonomies",
                        "BOLD",
                        "BOLD_raw",
                        "BOLD_11.csv"))
bold13 <- read.csv(here("data",
                        "raw_data",
                        "3_taxonomies",
                        "BOLD",
                        "BOLD_raw",
                        "BOLD_12.csv"))
bold14 <- read.csv(here("data",
                        "raw_data",
                        "3_taxonomies",
                        "BOLD",
                        "BOLD_raw",
                        "BOLD_13.csv"))
bold15 <- read.csv(here("data",
                        "raw_data",
                        "3_taxonomies",
                        "BOLD",
                        "BOLD_raw",
                        "BOLD_14.csv"))
bold16 <- read.csv(here("data",
                        "raw_data",
                        "3_taxonomies",
                        "BOLD",
                        "BOLD_raw",
                        "BOLD_15.csv"))
bold17 <- read.csv(here("data",
                        "raw_data",
                        "3_taxonomies",
                        "BOLD",
                        "BOLD_raw",
                        "BOLD_16.csv"))
bold18 <- read.csv(here("data",
                        "raw_data",
                        "3_taxonomies",
                        "BOLD",
                        "BOLD_raw",
                        "BOLD_17.csv"))

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
write.csv(all_bold, here("2c_taxonomic_assignment", 
                         "a_export_for_manual_entry",
                         "bold.csv"))

#I created a new CSV that includes all the taxonomic levels I compiled from internet
#searches
IDs_bold <- read.csv(here("data",
                          "raw_data",
                          "3_taxonomies",
                          "BOLD",
                          "bold_wID.csv"))

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
#Find sources and matches/mismatches of taxonomies####
#####################

#This dataset now includes four types of IDs:
#1. Those which matched to NCBI and not BOLD
#2. Those which matched to BOLD but not NCBI
#3. Those which matched to both and are the same
#4. Those which matched to both and are different

#I am going to break the dataframe into these four component parts (for now)
#just because that makes it easier for me to see what the composition of
#these four types are

#While I'm at it, I'm also going to give each set of data a master "unique_ID" 
#that will be it's identifier later. I think that since most of the data came
#from NCBI, I'm going to defer to the NCBI ID here (and will use the database
#coverage supplement to support this decision)

#NCBI only
IDs_ncbi_only <- IDs %>%
  filter(ID_ncbi != "" & ID_bold == "") %>%
  mutate(unique_ID = ID_ncbi)

#get stats on this:
IDs_ncbi_only %>%
  tally()

IDs_ncbi_only %>%
  group_by(ID_level_ncbi) %>%
  tally()

#BOLD only
IDs_bold_only <- IDs %>%
  filter(ID_bold != "" & ID_ncbi == "") %>%
  mutate(unique_ID = ID_bold)

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

IDs_same_level <- IDs_both %>%
  filter(ID_level_ncbi == ID_level_bold) 

#some of these are the same, but with different 
#names or spellings (e.g. different genus, but
#same species, or with an "sp" after a family)
IDs_same_level %>%
  dplyr::select(ASV, ID_ncbi, ID_bold) %>% #select columns for easy viewing
  filter(ID_ncbi != ID_bold) #look at mismatches to find real vs. not real

#based on this, only two ASVs have been assigned to different things
#in the two databases:
#ASV_1123
#ASV_1171

same_unmatched <- c("ASV_1123", "ASV_1171") #mismatched taxnomies

IDs_same_level <- IDs_same_level %>%
  filter(!ASV %in% same_unmatched) %>% #remove mismatches
  mutate(unique_ID = ID_ncbi) #add Unique_ID based on ID_ncbi, since it gave more data

#now all I have left are things not matched to the same level:
IDs_diff_level <- IDs_both %>%
  filter(ID_level_ncbi != ID_level_bold) 

#View(IDs_diff_level) to see unmatched
IDs_diff_level %>%
  dplyr::select(ASV, ID_ncbi, ID_bold)

#These are the differences 
diff_unmatched <- c("ASV_123", "ASV_186", "ASV_872", "ASV_1230", "ASV_1027")

IDs_diff_level <- IDs_diff_level %>%
  filter(!ASV %in% diff_unmatched) %>% #remove mismatched taxonomies
  mutate(unique_ID = ID_ncbi) #create a unique ID based on NCBI

#####################
#Combine all consistent taxonomies and create a master taxonomic ID column####
#####################

#combine final taxonomic assignments
all_IDs <- IDs_ncbi_only %>%
  bind_rows(IDs_bold_only) %>%
  bind_rows(IDs_same_level) %>%
  bind_rows(IDs_diff_level) 

#Now let's make a master level of every taxonomic assignment based on NCBI first
all_IDs$Domain <- ifelse(all_IDs$Domain_ncbi == "", all_IDs$Domain_bold, 
                         all_IDs$Domain_ncbi)
all_IDs$Phylum <- ifelse(all_IDs$Phylum_ncbi == "", all_IDs$Phylum_bold, 
                         all_IDs$Phylum_ncbi)
all_IDs$Class <- ifelse(all_IDs$Class_ncbi == "", all_IDs$Class_bold, 
                         all_IDs$Class_ncbi)
all_IDs$Order <- ifelse(all_IDs$Order_ncbi == "", all_IDs$Order_bold, 
                         all_IDs$Order_ncbi)
all_IDs$Family <- all_IDs$Family_ncbi #since nothing was below order level for BOLD
all_IDs$Genus <- all_IDs$Genus_ncbi #nothing below order for BOLD only
all_IDs$Species <- all_IDs$Species_ncbi #nothing below order for BOLD only

#going to keep the BOLD IDs in, since they are more specific, so they can be used
#to sort out predator DNA when we sort DNA, but still deferring to the NCBI
#ID for the unique_ID assigned to each ASV
all_IDs <- all_IDs %>%
  dplyr::select(ASV, unique_ID, Domain, Phylum, Class, Order, Family, Genus, Species,
                Class_bold, Order_bold, Family_bold, Genus_bold, Species_bold, ID_bold) %>%
  mutate(Domain = ifelse(Domain == "", "Eukaryota", Domain))

all_IDs$ID_level <- ifelse(all_IDs$Species != "", "Species",
                                 ifelse(all_IDs$Genus != "", "Genus",
                                        ifelse(all_IDs$Family != "", "Family",
                                               ifelse(all_IDs$Order != "", "Order",
                                                      ifelse(all_IDs$Class != "", "Class",
                                                             ifelse(all_IDs$Phylum != "", "Phylum", "Domain"))))))

#Summaries of these data now
all_IDs %>%
  tally()

all_IDs$ID_level <- factor(all_IDs$ID_level, levels = c("Domain", "Phylum",
                                                           "Class", "Order", "Family",
                                                           "Genus", "Species"))

all_IDs %>%
  group_by(ID_level) %>%
  tally()

# Export taxonomic assignments --------------------------------------------

write.csv(all_IDs, here("data", 
                        "outputs", 
                        "2c_taxonomic_assignment",
                        "a_export_for_manual_entry",
                        "ASV_taxonomies.csv"))

#####################
#Post-individual BLAST subsetting ####
#####################
#I blasted individual unmatched sequences on BLAST and assigned to a family level

taxa <- read.csv(here("data", 
                      "outputs",
                      "2c_taxonomic_assignment",
                      "b_import_after_manual_entry",
                      "ASV_taxonomies_wIndiv.csv"))

taxa <- taxa %>%
  filter(unique_ID !="non-diet")

#108 additional taxonomies to the family level from the sorting done 
#at individual level
taxa %>%
  filter(BLAST_unique_ID != "") %>%
  tally()

taxa_sort <- taxa %>%
  dplyr::select(ASV, Domain, Phylum, Class, Order, Family, Genus, Species, 
                Class_bold, Order_bold, Family_bold, Genus_bold, Species_bold) %>%
  mutate(Class = ifelse(Class != "", Class, Class_bold),
         Order = ifelse(Order != "", Order, Order_bold),
         Family = ifelse(Family != "", Family, Family_bold),
         Genus = ifelse(Genus != "", Genus, Genus_bold),
         Species = ifelse(Species != "", Species, Species_bold)) %>%
  dplyr::select(ASV, Domain, Phylum, Class, Order, Family, Genus, Species)

taxa_sort$ID_level <- ifelse(taxa_sort$Family != "", "Family",
                             ifelse(taxa_sort$Order != "", "Order",
                                    ifelse(taxa_sort$Class != "", "Class",
                                           ifelse(taxa_sort$Phylum != "", "Phylum", "Domain"))))

taxa_sort %>%
  tally()

taxa_sort %>%
  group_by(ID_level) %>%
  summarise(proportion = n()/781)

#####################
#Export ####
#####################

write.csv(taxa_sort, here("data", 
                          "outputs",
                          "2c_taxonomic_assignment",
                          "c_final_dataset",
                          "ASV_taxonomies_summed_wIndiv.csv"))

####################
#See how many ASVs are assigned target taxonomies####
#####################

all_ASVS <- read.csv(here("data", 
                          "raw_data",
                          "1_denoised_data", 
                          "dada2",  
                          "ASVs_counts_all.tsv"), sep = "\t")

all_ASVS <- all_ASVS %>%
  dplyr::select(X) %>%
  rename("ASV" = "X")

all <- all_ASVS %>%
  tally()

assigned <- all_IDs %>%
  tally()

#Total assigned to potential prey:
assigned/all #45% assigned to potential prey



