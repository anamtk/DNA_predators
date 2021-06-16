##########################
# Predator mass-length relationships-----
# Ana Miller-ter Kuile
# October 7, 2020
###########################

# this script builds species-specific mass-length
#relationships for each predator species in my 
#samples so I can convert lengths to masses for analyses

###########################
# Load packages-----
package.list <- c("here", "tidyverse", "glmmTMB", 
                  "MuMIn", "emmeans", "DHARMa", "stringr", 
                  "ggeffects", "fuzzyjoin")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}
#############################

#############################
# Load and combine datasets --------
#############################

pred_size <- read.csv(here("data", 
                           "outputs", 
                           "2g_master_body_size_lists", 
                           "pred_mass_length.csv"))

pred_size <- pred_size %>%
  mutate(Family = ifelse(Order == "Orthoptera", 
                         "Tettigonidae", Family)) %>%
  mutate(Family = ifelse(Family == "Tettigonidae", 
                         "Tettigoniidae", Family))

pred_id <- read.csv(here("data", 
                         "raw_data", 
                         "2_sample_data",
                         "Predator_IDs.csv"))

pred_id <- pred_id %>%
  dplyr::select(pred_Family, sample_str, hunting_mode, venom, webs)

pred_size <- pred_size %>%
  left_join(pred_id, by = c("Family" = "pred_Family"))

meta <- read.csv(here("data", 
                      "raw_data",
                      "2_sample_data",
                      "Sample_metadata.csv"))

#############################
# Model of relationship --------
#############################
#The model should have random intercepts by species and random slopes such that the relationship
#is allowed to vary by species
ggplot(pred_size, aes(x = Length_mm, y = Mass_mg, color = Source)) +
  geom_point() +
  facet_wrap(~sample_str, scale = "free")
#mutate to log10- log10 relationship
pred_size <- pred_size %>%
  mutate(log_mass = log10(Mass_mg),
         log_length = log10(Length_mm))

ggplot(pred_size, aes(x = log_length, y = log_mass, color = sample_str)) +
  geom_point() +
  facet_wrap(~sample_str, scale = "free") 

m1 <- glmmTMB(log_mass ~ log_length + (log_length|sample_str),
              data = pred_size)

#very strong relationship!
summary(m1)
r.squaredGLMM(m1)

#predict for pretty graph
me <- ggpredict(m1, terms = c("log_length", "sample_str"), type = "random")

#graph the predicted values
plot(me, add.data = TRUE) +
  facet_wrap(~group)

plot(me, add.data = TRUE) 

#############################
# Predict for my predators using model --------
#############################
#take metadata on predators and attach sample_str
#then convert log length and remove the scorpions

predators <- meta %>% 
  fuzzy_inner_join(pred_id, by = c("Extraction.ID" = "sample_str"), 
                   match_fun = str_detect) %>%
  dplyr::select(Extraction.ID, Length_mm, sample_str, 
                hunting_mode, venom, webs) %>%
  mutate(log_length = log10(Length_mm)) %>%
  filter(!sample_str == "ISO")
max(predators$log_length, na.rm = T)
max(predators$Length_mm, na.rm= T)
#now predict the masses of each of those predators based on model 1
log_mass <- predict(m1, newdata = predators)

#add that to the dataframe
predators <- predators %>%
  mutate(log_mass = log_mass)

#now they are values matched to the lines of the length-mass relationship
ggplot() +
  geom_point(data = pred_size, aes(x = log_length, y = log_mass), 
             color = "grey", alpha = 0.6) +
  geom_point(data = predators, aes(x = log_length, y = log_mass),
             color = "black") +
  geom_line(data = predators, aes(x = log_length, y = log_mass,
                                  color = sample_str),
            size = 1, alpha = 0.6) +
  theme_bw() +
  facet_wrap(~sample_str, scales="free")

#############################
# Export for later --------
#############################
  
write.csv(predators, here("data", 
                          "outputs", 
                          "2h_pred_mass_length", 
                          "DNA_pred_mass_length.csv"))

#############################
# Vis by species of body size (maybe needed later?) --------
#############################
ggplot(pred_size, aes(x = Length_mm, y = Mass_mg, color = sample_str)) +
  geom_smooth(aes(color = sample_str), method = "lm", se =F) +
  geom_point(size = 3) +
  theme_bw() +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~sample_str, scale = "free")

ggplot(pred_size, aes(x = Length_mm, y = Mass_mg, color = sample_str)) +
  geom_smooth(aes(color = sample_str), method = "lm", se =F) +
  geom_point(size = 3) +
  theme_bw() +
  scale_x_log10() +
  scale_y_log10() 

pred_size %>%
  filter(!is.na(Mass_mg)) %>%
  group_by(sample_str) %>%
  summarise(total = n())

pred_size %>%
  filter(!is.na(Mass_mg)) %>%
  summarise(total = n())
  


