##########################
# 1. What are predator-prey size relationships -----
# Ana Miller-ter Kuile
# October 8, 2020
###########################

# this script analyzes predator-prey body size relationships,
#asking the question:
#1. Does species identity, body size, or 
#both determine prey size across 
#all predator sizes?


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


# Size relationship quick visualizations -----------------------------------------------

#Does predator identity or size determine prey size?
#mean of prey species
ggplot(size, aes(x = pred_mass_mg, y = mean_prey_mass_mg, color = sample_str)) +
  geom_abline(slope = 1, linetype = "dashed") +
  geom_point(size = 3) +
  #scale_x_log10() +
  #scale_y_log10() +
  theme_bw() +
  facet_wrap(~sample_str, scale = "free")

ggplot(size, aes(x = pred_mass_mg, y = mean_prey_mass_mg, color = sample_str)) +
  geom_abline(slope = 1, linetype = "dashed") +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se =F) +
  #scale_x_log10() +
  #scale_y_log10() +
  theme_bw()# +
# facet_wrap(~sample_str, scale = "free")

#min of prey species
ggplot(size, aes(x = pred_mass_mg, y = min_prey_mass_mg, color = sample_str)) +
  geom_abline(slope = 1, linetype = "dashed") +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se =F) +
  #scale_x_log10() +
  #scale_y_log10() +
  theme_bw() 

ggplot(size, aes(x = pred_mass_mg, y = min_prey_mass_mg, color = sample_str)) +
  geom_abline(slope = 1, linetype = "dashed") +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se =F) +
  #scale_x_log10() +
  #scale_y_log10() +
  theme_bw() +
  facet_wrap(~sample_str, scale = "free")

# Body size model mean ---------------------------------------------------------

#is the prey size determined by some combination of predator identity
#and predator size?
m1 <- glmmTMB(mean_prey_log_mass_mg ~ pred_log_mass_mg*sample_str + (1|sample),
              data = size,
              REML = FALSE)

#interaction model: both slope and intercept vary by species
#mass + species model: only intercept varies by species
#mass model: only mass matters
#species model: only species matters

dredge(m1)

#best is mass + species model
m2 <- glmmTMB(mean_prey_log_mass_mg ~ pred_log_mass_mg + sample_str + (1|sample),
              data = size)

fit <- simulateResiduals(m2, plot = T)

me <- ggpredict(m2, terms = c("pred_log_mass_mg", "sample_str"), type = "random")
plot(me, add.data = TRUE) +
  geom_abline(slope = 1, linetype = "dashed") +
  facet_wrap(~group) 

plot(me, add.data = TRUE, ci = FALSE) +
  geom_abline(slope = 1, linetype = "dashed") 

#in this summary,
#intercept Estimate gives that base intercept for the model
#the intercept for each species is the sum of that and the 
#species intercept
#the power law relationship is given by the estimate of 
#pred_log_mass_mg, which is sublinear with a relationship of
#y = a + x^0.41259
summary(m2)
r.squaredGLMM(m2)

# Body size model min ---------------------------------------------------------

m3 <- glmmTMB(min_prey_log_mass_mg ~ pred_log_mass_mg*sample_str + (1|sample),
              data = size,
              REML = FALSE)

dredge(m3)

m4 <- glmmTMB(min_prey_log_mass_mg ~ pred_log_mass_mg + sample_str + (1|sample),
              data = size)

fit <- simulateResiduals(m4, plot = T)

me2 <- ggpredict(m4, terms = c("pred_log_mass_mg", "sample_str"))
plot(me2, add.data = TRUE) +
  geom_abline(slope = 1, linetype = "dashed") +
  facet_wrap(~group) 

plot(me2, add.data = TRUE, ci=FALSE) +
  geom_abline(slope = 1, linetype = "dashed")

#in this summary,
#intercept Estimate gives that base intercept for the model
#the intercept for each species is the sum of that and the 
#species intercept
#the power law relationship is given by the estimate of 
#pred_log_mass_mg, which is sublinear with a relationship of
#y = a + x^0.26466
summary(m4)
r.squaredGLMM(m4)





