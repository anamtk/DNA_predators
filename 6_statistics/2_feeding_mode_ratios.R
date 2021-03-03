##########################
# 2. Are prey:predator ratios explained by feeding mode? -----
# Ana Miller-ter Kuile
# October 29, 2020
###########################

# this script analyzes whether predator feeding
#mode (active vs. non-active) influences the
#size of prey they eat beyond just a 
#species-specific size mode. OR, maybe 
#there is somethign that is varyign at the 
#species level that matters instead, regardless
#of hunting mode. 

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 
                  "glmmTMB", "emmeans",
                  "MuMIn", "DHARMa",
                  "effects", "ggeffects")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Import data -------------------------------------------------------------

data <- read.csv(here("data", "outputs", "8_final_dataset",
                      "pred_prey_sizes_tp_DNAinteractions.csv"))

size <- data %>%
  dplyr::select(-X, -X.x, -X.y, -reads) %>%
  mutate(pred_mass_mg = exp(pred_log_mass_mg)) %>%
  mutate(mean_ratio = mean_prey_mass_mg/pred_mass_mg,
         min_ratio = min_prey_mass_mg/pred_mass_mg,
         log_mean_ratio = log(mean_ratio),
         log_min_ratio = (log(min_ratio)))

# Size relationship visualizations -----------------------------------------------
#mean ratio
ggplot(size, aes(x = sample_str, y = mean_ratio, fill = Feeding_mode)) +
  geom_hline(yintercept  = 1, linetype = "dashed") +
  geom_boxplot() +
  scale_y_log10() +
  theme_bw()

#min ratio
ggplot(size, aes(x = sample_str, y = min_ratio, fill= Feeding_mode)) +
  geom_hline(yintercept  = 1, linetype = "dashed") +
  geom_boxplot() +
  scale_y_log10() +
  theme_bw()

# Mean size ratio model -----------------------------------------------
#does feeding mode or species influence influence the ratio of prey size to predator size?
m5 <- glmmTMB(log_mean_ratio ~ sample_str + Feeding_mode + (1|sample),
              data = size, 
              REML = FALSE)

#The full model asks if a combination of species
#identity and feeding mode predict size relationships
dredge(m5)

#best model says that just species matters for this relationship
#regardless of feeding mode
m6 <- glmmTMB(log_mean_ratio ~ sample_str + (1|sample),
              data = size)

fit <- simulateResiduals(m6, plot = T)

me <- ggpredict(m6, terms = c("sample_str"), type = "random")
plot(me, add.data = TRUE) 

summary(m6)
r.squaredGLMM(m6)

# Min size ratio model -----------------------------------------------
#does feeding mode or species influence influence the ratio of prey size to predator size?
m7 <- glmmTMB(log_min_ratio ~ Feeding_mode + sample_str + (1|sample),
              data = size, 
              REML = FALSE)

#The full model asks if a combination of species
#identity and feeding mode predict size relationships
dredge(m7)

#best model says that just species matters for this relationship
#regardless of feeding mode
m8 <- glmmTMB(log_min_ratio ~ sample_str + (1|sample),
              data = size)

fit <- simulateResiduals(m8, plot = T)

me <- ggpredict(m8, terms = c("sample_str"), type = "random")
plot(me, add.data = TRUE)

summary(m8)
r.squaredGLMM(m8)
