##########################
# 1. What are predator-prey size relationships -----
# Ana Miller-ter Kuile
# October 8, 2020
###########################

# this script analyzes predator-prey body size relationships,
#asking the question, first:
#1. do larger predator individuals across species eat larger prey?
#(both with mean and min prey family sizes)
#and then asking 
#2. For predator species which we have large sample sizes,
#do larger predator individuals within species eat larger prey?
#(both with mean and min prey family sizes)

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
  mutate(pred_mass_mg = exp(pred_log_mass_mg)) 


# Size relationship visualizations -----------------------------------------------

size %>%
  distinct(sample, sample_str) %>%
  group_by(sample_str) %>%
  summarise(sample_size = n())

size %>%
  distinct(sample_str, Family) %>%
  group_by(sample_str) %>%
  summarise(sample_size = n())

#Do larger predators eat larger prey?
#mean of prey species
ggplot(size, aes(x = pred_mass_mg, y = mean_prey_mass_mg, color = sample_str)) +
  geom_abline(slope = 1, linetype = "dashed") +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se =F) +
  #scale_x_log10() +
  #scale_y_log10() +
  theme_bw() +
  facet_wrap(~sample_str, scale = "free")

#min of prey species
ggplot(size, aes(x = pred_mass_mg, y = min_prey_mass_mg, color = sample_str)) +
  geom_abline(slope = 1, linetype = "dashed") +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se =F) +
  #scale_x_log10() +
  #scale_y_log10() +
  theme_bw() +
  facet_wrap(~sample_str, scale = "free")

# All species Body size model mean ---------------------------------------------------------
hist(size$mean_prey_log_mass_mg)
m1 <- glmmTMB(mean_prey_log_mass_mg ~ pred_mass_mg + (1|sample_str) + (1|sample),
              data = size)
summary(m1)

dredge(m1)

plot(allEffects(m1))

fit <- simulateResiduals(m1, plot = T)

me <- ggpredict(m1, terms = c("pred_mass_mg", "sample_str"), type = "random")
plot(me) 

r.squaredGLMM(m1)

# All species Body size model min ---------------------------------------------------------

m2 <- glmmTMB(min_prey_log_mass_mg ~ pred_mass_mg + (1|sample_str) + (1|sample),
              data = size)
summary(m2)
dredge(m2)
plot(allEffects(m2))
fit <- simulateResiduals(m2, plot = T)

me <- ggpredict(m2, terms = c("pred_mass_mg", "sample_str"), type = "random")
plot(me) 

r.squaredGLMM(m2)


# subset species Body size model mean ---------------------------------------------------------

species <- size %>%
  filter(sample_str %in% c("HEV", "PHH", "NEO"))

hist(species$mean_prey_log_mass_mg)

m3 <- glmmTMB(mean_prey_log_mass_mg ~ pred_mass_mg*sample_str + (1|sample),
              data = species)
summary(m3)

dredge(m3)

plot(allEffects(m3))

fit <- simulateResiduals(m3, plot = T)

me <- ggpredict(m3, terms = c("pred_mass_mg", "sample_str"), type = "random")
plot(me) 

r.squaredGLMM(m3)

# subset species Body size model min ---------------------------------------------------------

hist(species$min_prey_log_mass_mg)

species %>%
  group_by(sample_str) %>%
  summarise(range = max(pred_mass_mg) - min(pred_mass_mg))

m4 <- glmmTMB(min_prey_log_mass_mg ~ pred_mass_mg*sample_str + (1|sample),
              data = species)
summary(m4)

dredge(m4)

plot(allEffects(m4))

fit <- simulateResiduals(m4, plot = T)

r.squaredGLMM(m4)

