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
                  "effects", "ggeffects",
                  "calecopal", "patchwork",
                  "emmeans")

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


# How many preds? -----------------------------------------------
size %>%
  distinct(sample) %>%
  summarise(total = n())

# Body size model selection ---------------------------------------------------------

#is the prey size determined by some combination of predator identity
#and predator size?
m1 <- glmmTMB(mean_prey_log_mass_mg ~ pred_log_mass_mg*sample_str + (1|sample),
              data = size,
              REML = FALSE)

#interaction model: both slope and intercept vary by species
#mass + species model: slope invariant, only intercept varies by species
#mass model: only mass matters
#species model: only species matters

dredge(m1)

#best is mass + species model
m2 <- glmmTMB(mean_prey_log_mass_mg ~ pred_log_mass_mg + sample_str + (1|sample),
              data = size)

fit <- simulateResiduals(m2, plot = T)

#in this summary,
#intercept Estimate gives that base intercept for the model
#the intercept for each species is the sum of that and the 
#species intercept
#the power law relationship is given by the estimate of 
#pred_log_mass_mg, which is sublinear with a relationship of
#y = a + x^0.41259


# Model diagnostics -------------------------------------------------------

summary(m2)
r.squaredGLMM(m2)

em <- emmeans(m2, "sample_str")
tukey <- as.data.frame(pairs(em))

tukey <- tukey %>%
  mutate(sig = ifelse(p.value < 0.05, "sig", "non-sig"))

# Figures ------------------------------------------------------------------

pal_kelp <- cal_palette("kelp1", n = 9, type = "continuous")

pred_labels <- c("CEN" = "Geophilomorpha sp.", "EUB" = "E. annulipes", 
                 "HEV" = "H. venatoria", "LRS" = "Oonopidae sp.", 
                 "NEO" = "N. theisi", "PAN" = "P. flavescens",
                 "PHH" = "P. holdhausi", "SCY" = "S. longipes",
                 "SME" = "S. pallidus")

#pred size distribution on its own
(pred_size2 <- size %>%
    distinct(sample, pred_mass_mg) %>%
    ggplot(aes(x = pred_mass_mg)) +
    geom_histogram(bins = 50, alpha = 0.85, color = "black") +
    #scale_y_log10() +
    xlim(-8, 7) +
    #geom_density(fill = "#BE8333", alpha = 0.85) +
    scale_x_log10() +
    theme_bw() +
    labs(x = "Predator mass (mg)", 
         y = "Individuals") +
    theme(axis.text = element_text(size =20),
          axis.title = element_text(size = 25)))

#prey size distribution on its own
(prey_mean_size2 <- size %>%
    distinct(Family, mean_prey_mass_mg) %>%
    ggplot(aes(x = mean_prey_mass_mg)) +
    geom_histogram(bins = 50, alpha = 0.85, color = "black") +
    scale_x_log10() +
    #geom_density(fill = "#BE8333", alpha = 0.85) +
    theme_bw() +
    labs(x = "Mean prey mass (mg)", 
         y = "Prey families") +
    theme(axis.text = element_text(size =20),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_text(size = 25)))

#Does predator identity or size determine prey size?
size_graph_col <- size %>%
  mutate(sample_str = factor(sample_str, levels = c("LRS", "SCY", "NEO", "CEN",
                                                    "SME", "EUB", "PHH", "PAN",
                                                    "HEV"))) %>%
  ggplot(aes(x = pred_mass_mg, y = mean_prey_mass_mg, color = sample_str)) +
  geom_abline(slope = 1, linetype = "dashed", size = 0.75) +
  geom_point(size = 3) +
  geom_abline(slope = 0.41, size = 1) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Predator mass (mg)", 
       y = "Prey mass (mg)",
       color = "Predator species") +
  scale_color_manual(values = pal_kelp,
                     labels = pred_labels) +
  theme_bw() +
  theme(axis.text = element_text(size =20),
        axis.title = element_text(size = 25))
size_graph_col
#without colors:
size_graph_ncol <- size %>%
  mutate(sample_str = factor(sample_str, levels = c("LRS", "SCY", "NEO", "CEN",
                                                    "SME", "EUB", "PHH", "PAN",
                                                    "HEV"))) %>%
  ggplot(aes(x = pred_mass_mg, y = mean_prey_mass_mg)) +
  geom_abline(slope = 1, linetype = "dashed", size = 0.75) +
  geom_point(size = 3) +
  geom_abline(slope = 0.41, size = 1) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Predator mass (mg)", y = "Prey mass (mg)") +
  theme_bw() +
  theme(axis.text = element_text(size =20),
        axis.title = element_text(size = 25))

#sorted by increased average predator size, then showing prey size
species_graph <- size %>%
  mutate(sample_str = fct_reorder(sample_str, pred_mass_mg, .fun='mean')) %>%
  ggplot(aes(x = reorder(sample_str, pred_mass_mg), y = mean_prey_mass_mg, color = sample_str)) +
  geom_boxplot(size = 1) + 
  geom_point() +
  scale_color_manual(values = pal_kelp,
                     labels = pred_labels) +
  theme_bw() +
  scale_y_log10() +
  labs(x = "Predator species", 
       y = "Prey mass (mg)",
       color = "Predator species") +
  theme(axis.text.y = element_text(size =20),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 25))

#pairwise comparisons
pairwise_sp <- tukey %>%
  arrange(estimate, p.value) %>%
  mutate(contrast = factor(contrast, levels = contrast)) %>%
ggplot(aes(x = contrast, y = estimate, color = sig)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_point(size = 2) +
  geom_linerange(aes(ymin = estimate - SE, 
                     ymax = estimate + SE), 
                 size = 1) +
  scale_color_manual(values = c("sig" = "#229AAA",
                                "non-sig" = "#878787")) +
  theme_bw() +
  coord_flip() +
  labs(x = "Pairwise contrast", y = "Difference") +
  theme(legend.position = "none",
        axis.text = element_text(size =15),
        axis.title = element_text(size = 25))

