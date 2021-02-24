

data <- read.csv(here("data", "outputs", "5_rarefied_taxonomic_sort", "fam_prey_DNA_conservative.csv"))

prey <- data %>%
  distinct(Domain, Phylum, Class, Order, Family)

prey_fams <- prey %>%
  dplyr::select(Family)

sizes <- read.csv(here("data", "size_data", "Palmyra_BS_Sep2020.csv"))
str(sizes)

unique(sizes$Family)

sizes <- sizes %>%
  dplyr::select(Morphospecies, 
                Phylum, 
                Class, 
                Order.1, 
                Family,
                Body_Length_Mean_mm,
                Body_Mass_Mean_mg,
                Body_Mass_SD_mg,
                Body_Mass_N)

missing <- sizes %>%
  filter(is.na(Body_Mass_Mean_mg)) %>%
  distinct(Morphospecies, 
           Phylum, 
           Class, 
           Order.1, 
           Family,
           Body_Length_Mean_mm,
           Body_Mass_Mean_mg,
           Body_Mass_SD_mg,
           Body_Mass_N) %>%
  filter(Phylum != "")
str(missing)

missing_sp <- missing %>%
  semi_join(prey_fams, by = "Family")

write.csv(missing_sp, here("data", "outputs", "palmyra_bs_for_chelsea_feb162021.csv"))





