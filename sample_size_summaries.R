library(tidyverse)
library(ggplot2)
library(here)

data <- read.csv(here("Sample_metadata.csv"))
str(data)
data$Date.Collected <- as.character(data$Date.Collected)
data$Extraction.ID <- as.character(data$Extraction.ID)
data$Isotope_ID <- as.character(data$Isotope_ID)
str(data)

data %>%
  group_by(Year, Island) %>%
  summarize(sum = n())

species_num <- data %>%
  group_by(Date.Collected, Island, Habitat, Method) %>%
  summarize(species_num = n_distinct(ID)) 

indiv_num <- data %>%
  group_by(Date.Collected, Island, Habitat, Method) %>%
  summarize(indiv_num = n())

communities <- species_num %>%
  left_join(indiv_num, by = c("Date.Collected", "Island", "Habitat", "Method"))

fogging_comm <- communities %>%
  filter(Method == "Fogging") %>%
  unite("Tree_ID", Island, Habitat, Date.Collected, remove = FALSE)

soil <- data %>%
  filter(ID %in% c("Centipede", "Euborellia annulipes"))

microhabitats <- data %>%
  filter(ID %in% c("Euborellia annulipes", "Heteropoda venatoria")) %>%
  unite("Micro_Specific", Habitat, Microhabitat, remove = FALSE)

ggplot(microhabitats, aes(x = Length_mm, color = Microhabitat)) +
  geom_freqpoly(size = 2) +
  facet_wrap(~ID) +theme_bw() +
  scale_color_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
                                "#fb9a99","#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a",
                                "#ffff99")) 

soil_method <- data %>%
  filter(Method == "Soil")
soil$Method[which(soil$Method == "Hand")] <- "Soil"

#number of areas each species was found
data %>%
  filter(Source == "FIELD") %>%
  filter(Sterilized == "NS" | is.na(Sterilized)) %>%
  group_by(ID) %>%
  summarize(communities = n_distinct(Microhabitat))

data <- data %>%
  filter(Source == "FIELD") %>%
  filter(Sterilized == "NS" | is.na(Sterilized))

hab_generalists <- soil %>%
  group_by(ID, Habitat, Method) %>%
  summarize(individuals = n())

species_indivs <- data %>%
  group_by(ID) %>%
  summarize(individuals = n())

data1 <- data %>%
  left_join(species_indivs, by = "ID")

hist(data$Length_mm)

ggplot(data1, aes(x = Length_mm, fill = ID)) +
  geom_histogram(binwidth = 0.5) +
  scale_fill_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
    "#fb9a99","#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a",
    "#ffff99")) + theme_bw() +
  facet_wrap(~ID) +
  geom_text(aes(x = 39, y = 15, label = individuals), color = "black") +
  annotate("text", x = 33, y = 15, label = "N =") +
  labs(x = "Predator body length (mm)", y = "Individual count", 
       title = "Predator body size distributions and sample sizes (N)")


ggplot(data, aes(x = Length_mm, after_stat(density), color = ID)) +
  geom_freqpoly(size = 2) +
  scale_color_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
                               "#fb9a99","#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a",
                               "#ffff99")) + theme_bw()

ggplot(data, aes(x = Length_mm)) +
  geom_histogram()  + theme_bw()

fogging <- data %>%
  filter(Method == "Fogging") %>%
  group_by(Island, Habitat, Date.Collected) %>%
  mutate(group_id = group_indices()) %>%
  unite("Environment", Island, Habitat, remove = FALSE) %>%
  unite("Tree", group_id, Environment, remove = FALSE) %>%
  unite("Tree_ID", Island, Habitat, Date.Collected, remove = FALSE) %>%
  left_join(fogging_comm, by = "Tree_ID")

ggplot(fogging, aes(x = Length_mm, color = ID)) +
  geom_freqpoly(size = 2, alpha = 0.8) +
  scale_color_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
                                "#fb9a99","#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a",
                                "#ffff99")) +
  facet_wrap(~Tree) + theme_bw() +
  geom_text(aes(x = 17, y = 10, label = species_num), color = "black") +
  geom_text(aes(x = 17, y = 12, label = indiv_num), color = "black") +
  annotate("text", x = 15, y = 10, label = "SR =  ") +
  annotate("text", x = 15, y = 12, label = " N =  ") +
  labs(x = "Predator Body Length (mm)", y = "Number of Individuals", 
       title = "Predator community density and body size distributions of 6 tree canopies
       (Including total sample size (N) and species richness (SR))")

ggplot(fogging, aes(x = Length_mm)) +
  geom_density() +
  facet_wrap(~Tree) + theme_bw() +
  geom_text(aes(x = 17, y = 0.2, label = species_num), color = "black") +
  geom_text(aes(x = 17, y = 0.25, label = indiv_num), color = "black") +
  annotate("text", x = 15, y = 0.2, label = "SR = ") +
  annotate("text", x = 15, y = 0.25, label = "N =")
