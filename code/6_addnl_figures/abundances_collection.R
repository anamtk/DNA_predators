# Relative abundances in sample collection
# Ana Miller-ter Kuile
# September 1, 2021

#this code examines the relative abudnance (somewhat akin to abundance) of our predator 
# species studied and other species (predatory and not) in teh insect collection


# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", "readxl")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Import data -------------------------------------------------------------

data <- read.csv(here('data', 
                      "raw_data", 
                      "4_body_size_data",
                      "palmyra_nodes.csv"))

sheets <- here("data",
               "raw_data",
               "2_sample_data",
               "Palmyra_2017_Arthropod_Collection.xlsx")

tab_names <- excel_sheets(path = sheets)

list_all <- lapply(tab_names, function(x) read_excel(path = sheets, sheet = x))

all <-  bind_rows(list_all)

collection <- all %>%
  mutate(ID = case_when(ID %in% c("Soil centipede", "Centipede", "Soil Centipede") ~ "Mecistocephalus sp",
                        ID == "LRS-2" ~ "Opopaea sp",
                        TRUE ~ ID)) %>%
  group_by(Order, ID) %>%
  summarise(total = n()) %>%
  filter(total > 2) %>%
  mutate(studied = case_when(ID %in% c("Neoscona theisi", "Phisis holdhausi",
                                       "Isometrus maculatus", "Euborellia annulipes",
                                       "Scytodes sp.", "Opopaea sp",
                                       "Heteropoda venatoria", "Smeringopus pallidus",
                                       "Mecistocephalus sp", "Pantala flavescens") ~ "yes",
                             TRUE ~ "no")) %>%
  mutate(abund = case_when(total < 10 ~ "lower abundance",
                               TRUE ~ "higher abundance")) %>%
  mutate(category = case_when(abund == "higher abundance" & ID %in% c("Neoscona theisi", "Phisis holdhausi",
                                                                       "Isometrus maculatus", "Euborellia annulipes",
                                                                       "Scytodes sp.", "Opopaea sp",
                                                                       "Heteropoda venatoria", "Smeringopus pallidus",
                                                                       "Mecistocephalus sp", "Pantala flavescens",
                                                                       "Juvenile Salticidae", "Ara07",
                                                                       "Chrysosoma sp.", "Ara06", "LRS-3") ~ "predatory",
                               abund == "lower abundance" ~ "less than 10 individuals",
                               TRUE ~ "not predatory, unknown, or no species ID")) %>%
  filter(!is.na(ID)) %>%
  mutate(ID = fct_reorder(ID, desc(total)))
  


# Graph abundance in collection -------------------------------------------



ggplot(collection, (aes(x = reorder(ID, -total), y = total, fill = category, color = studied))) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_color_manual(values = c("white", "black")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "top") +
  labs(x = "Morphospecies ID", y = "Total samples")


# Summary stats of collection ---------------------------------------------

collection %>%
  filter(studied == "yes") %>%
  ungroup() %>%
  summarise(mean = mean(total),
            sd = sd(total),
            n = n(),
            se = sd/sqrt(n),
            min = min(total))
  
  