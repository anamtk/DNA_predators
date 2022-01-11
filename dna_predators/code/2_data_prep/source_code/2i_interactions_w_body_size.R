##########################
# Prey Average Body Masses -----
# Ana Miller-ter Kuile
# October 8, 2020
###########################

# this script summarises the prey mass by family
#and by order (and class) for the prey reference list, and then
#combines with the interaction DF, and then figures out what
#is missing

###########################
# Load packages-----
package.list <- c("here", "tidyverse")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}
#############################

#############################
#Import Data -------
#############################
#prey size references
prey_ref <- read.csv(here("data", 
                          "outputs", 
                          "2g_master_body_size_lists", 
                          "prey_mass_length.csv"))

#individual predator sizes
preds <- read.csv(here("data", 
                       "outputs", 
                       "2h_pred_mass_length", 
                       "DNA_pred_mass_length.csv"))

#interactions for each predator
ints <- read.csv(here("data", 
                      "outputs",
                      "2e_rarefied_taxonomic_sort",
                      "fam_prey_DNA_conservative.csv"))

#remove run ID
ints$sample <- str_sub(ints$sample, end=-2)

#get rid of 0 reads, and then combine by sample, class, order, family
#multiple ASVs with the same species assignemnt
ints <- ints %>%
  filter(reads > 0) %>%
  group_by(sample, sample_str, Class, Order, Family) %>%
  summarise(reads = sum(reads))
  
#############################
#Average and min reference prey size by family, order, and class
#############################

#for some orders (Entognatha, eg) I only have order-level match
#but for many I'll want to go with the Family-level summary
#I'll make each and then combine both with the final DF

ref_fam <- prey_ref %>%
  group_by(Order, Family) %>%
  summarise(mean_prey_length_mm = mean(Length_mm, na.rm =T),
            mean_prey_mass_mg = mean(Mass_mg, na.rm = T),
            min_prey_length_mm = min(Length_mm, na.rm =T),
            min_prey_mass_mg = min(Mass_mg, na.rm = T))

prey_ref %>%
  group_by(Source) %>%
  summarise(total = n())

prey_ref %>%
  summarise(total = n())

family_n <- prey_ref %>%
  group_by(Order, Family) %>%
  summarise(total = n())

ref_ord <- prey_ref %>%
  group_by(Order) %>%
  summarise(mean_prey_length_mm = mean(Length_mm, na.rm =T),
            mean_prey_mass_mg = mean(Mass_mg, na.rm = T),
            min_prey_length_mm = min(Length_mm, na.rm =T),
            min_prey_mass_mg = min(Mass_mg, na.rm = T))

ref_cls <- prey_ref %>%
  group_by(Class) %>%
  summarise(mean_prey_length_mm = mean(Length_mm, na.rm =T),
            mean_prey_mass_mg = mean(Mass_mg, na.rm = T),
            min_prey_length_mm = min(Length_mm, na.rm =T),
            min_prey_mass_mg = min(Mass_mg, na.rm = T))

#############################
#Combine prey sizes with interaction dataframe
#############################

#family first, and then match order or class to those not assigned
ints_fam <- ints %>%
  left_join(ref_fam, by = c("Order", "Family"))

#order level sizes
ints_ord <- ints_fam %>%
  filter(is.na(mean_prey_length_mm)) %>%
  dplyr::select(-mean_prey_length_mm, -mean_prey_mass_mg, 
                -min_prey_length_mm, - min_prey_mass_mg) %>%
  left_join(ref_ord, by = c("Order")) 

#who is still missing? (they are ones that have different classifications)
ints_ord %>%
  filter(is.na(mean_prey_length_mm)) %>%
  distinct(Class, Order, Family)

#give these manually, which is just the mites and collembola. small but IMPORTANT
ints_cls <- ints_ord %>%
  filter(is.na(mean_prey_length_mm)) %>%
  mutate(mean_prey_mass_mg = ifelse(Order == "Sarcoptiformes", 0.00383400, mean_prey_mass_mg)) %>%
  mutate(min_prey_mass_mg = ifelse(Order == "Sarcoptiformes", 0.000633000, min_prey_mass_mg)) %>%
  mutate(mean_prey_mass_mg = ifelse(Class == "Collembola", 0.03473734, mean_prey_mass_mg)) %>%
  mutate(min_prey_mass_mg = ifelse(Class == "Collembola", 0.003383756, min_prey_mass_mg)) 

#############################
#Merge prey size assignments into on DF ---------
#############################

ints_f <- ints_fam %>%
  filter(!is.na(mean_prey_mass_mg))

ints_o <- ints_ord %>%
  filter(!is.na(mean_prey_mass_mg))

ints_c <- ints_cls %>%
  filter(!is.na(mean_prey_mass_mg))

#a total of 335 interactions given size assignements
ints_size <- ints_f %>%
  bind_rows(ints_o) %>%
  bind_rows(ints_c)

#############################
#Merge prey size assignments with preds ---------
#############################

p_p_sizes <- preds %>%
  rename("sample" = "Extraction.ID") %>%
  group_by(sample, hunting_mode, venom, webs) %>%
  summarise(pred_length_mm = mean(Length_mm, na.rm=T),
            pred_log_length_mm = mean(log_length, na.rm=T),
            pred_log_mass_mg = mean(log_mass, na.rm=T)) %>%
  left_join(ints_size, by = c("sample")) %>%
  filter(!is.na(mean_prey_mass_mg)) %>%
  mutate(mean_prey_log_mass_mg = log10(mean_prey_mass_mg),
         min_prey_log_mass_mg = log10(min_prey_mass_mg))

p_p_sizes <- p_p_sizes %>%
  filter(sample != "EUB36")

#############################
#Export for analyses ---------
#############################

write.csv(p_p_sizes, here("data", 
                          "outputs",  
                          "2i_final_dataset", 
                          "pred_prey_sizes_DNAinteractions.csv"))









