library(bold)
library(here)
here()
library(tidyverse)

#read tab delimited file with headers
plants <- read.delim(here("Database_coverage", "plant_species-Oct19.txt"), header = FALSE, stringsAsFactors = FALSE, sep = "\t", quote = "\"")
animals <- read.delim(here("Database_coverage", "pal-species_nov19.txt"), header = FALSE, stringsAsFactors = FALSE, sep = "\t", quote = "\"")
families <- read.delim(here("Database_coverage", "pal-families.txt"), header = FALSE, stringsAsFactors = FALSE, sep = "\t", quote = "\"")

#check number of rows
nrow(plants)
nrow(animals)
nrow(families)

#get all genes for plants in list and put them in a tab delimited file and a dataframe called results. 
#All plant request too long. Would only take about 160 at a time. 
#The file appends, so you can run the list in blocks.
animals1 <- animals[1:160,]
animals2 <- animals[161:268,]

for (animalName in animals1){
  results <- bold_seqspec(taxon = animalName)
  write.table(results, file='results-BOLD-alltaxa.tsv', sep="\t", append = TRUE , row.names = F, col.names = TRUE)
}

for (animalName in animals2){
  results <- bold_seqspec(taxon = animalName)
  write.table(results, file='results-BOLD-alltaxa.tsv', sep="\t", append = TRUE , row.names = F, col.names = TRUE)
}

#check results
Bold_species <- read.delim("results-BOLD-alltaxa.tsv", header = FALSE, stringsAsFactors = FALSE, sep = "\t", quote = "\"")

#check number of rows
nrow(Bold_species)

colnames(Bold_species) <- as.character(Bold_species[1, ])

Bold_species <- Bold_species[-1,]

bold <- Bold_species %>%
  dplyr::select(institution_storing, phylum_name, class_name, order_name, 
         family_name, genus_name, species_name) %>%
  filter(!institution_storing == "Mined from GenBank, NCBI") %>%
  group_by(order_name, family_name, genus_name, species_name) %>%
  summarise(number_of_records = n()) %>%
  mutate(ID_bold = ifelse(species_name == "", genus_name, species_name))


coiNCBI <- read.delim(here("Database_coverage", "palmyra-coi.txt"), header = FALSE, stringsAsFactors = FALSE, sep = "\t", quote = "\"")

co1NCBI <- read.delim(here("Database_coverage", "palmyra-co1.txt"), header = FALSE, stringsAsFactors = FALSE, sep = "\t", quote = "\"")


colnames(animals) <- "species"

animals$genus <- sapply(str_split(animals$species, " "), function(x){return(x[[1]])})

genus <- bold %>%
  full_join(animals, by = c("genus_name" = "genus"))



