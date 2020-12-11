##########################
# 1. Do larger predators  have broader range of prey size -----
# Ana Miller-ter Kuile
# October 8, 2020
###########################

# this script analyzes predator-prey body size relationships,
#asking the question:
#Do the range of prey sizes vary by predator size within species?

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

# Range Predator size within an individual -----------------------------------
range_prey <- size %>%
  group_by(sample, sample_str, pred_log_mass_mg, pred_mass_mg) %>%
  summarise(range = max(mean_prey_mass_mg) - min(mean_prey_mass_mg)) %>%
  filter(range > 0) #%>%
  #filter(sample_str %in% c("HEV", "NEO", "PHH"))
hist(range_prey$range)
#some species now have very few individuals, so we should remove
#those here as well.
size %>%
  group_by(sample, sample_str, pred_mass_mg) %>%
  summarise(range = max(mean_prey_mass_mg) - min(mean_prey_mass_mg)) %>%
  filter(range > 0) %>%
  #filter(sample_str %in% c("HEV", "NEO", "PHH")) %>%
  ggplot(aes(x = pred_mass_mg, y = range, color = sample_str)) +
  geom_point(size = 2) +
  theme_bw() +
  facet_wrap(~sample_str, scales = "free")

# All species range prey -----------------------------------------------------
range_prey <- range_prey %>%
  filter(!sample_str %in% c("CEN", "EUB", "LRS")) %>%
  mutate(log_range = log(range))

m <- glmmTMB(log_range ~ pred_log_mass_mg*sample_str,
             data = range_prey,
             na.action = "na.fail")

dredge(m)

m2 <- glmmTMB(log_range ~ pred_log_mass_mg + sample_str,
             data = range_prey,
             na.action = "na.fail")

fit <- simulateResiduals(m2, plot = T)

summary(m2)

# Visualizations ----------------------------------------------------------
x <- 1:100
y <- x^(0.35911)
plot(y ~ x)

pal_kelp <- cal_palette("kelp1", n = 9, type = "continuous")
pal_kelp
#LRS: "#C70000" 
# SCY: "#EA7700" 
#NEO: "#EEB00C" 
#CEN: "#C68A2C" 
#SME: "#89742F" 
#EUB: "#496C3C" 
#PHH: "#158D8E"
#PAN: "#067D8D" 
#HEV: "#114C54"

pred_labels <- c("CEN" = "Geophilomorpha sp.", "EUB" = "E. annulipes", 
                 "HEV" = "H. venatoria", "LRS" = "Oonopidae sp.", 
                 "NEO" = "N. theisi", "PAN" = "P. flavescens",
                 "PHH" = "P. holdhausi", "SCY" = "S. longipes",
                 "SME" = "S. pallidus")

range_graph <- size %>%
  group_by(sample, sample_str, pred_mass_mg) %>%
  summarise(range = max(mean_prey_mass_mg) - min(mean_prey_mass_mg)) %>%
  filter(range > 0) %>%
  filter(!sample_str %in% c("CEN", "EUB", "LRS")) %>%
  #mutate(sample_str = factor(sample_str, 
                             #levels = c("NEO", "PHH", "HEV"))) %>%
  ggplot(aes(x = pred_mass_mg, y = range, color = sample_str)) +
  geom_abline(slope =1, linetype = "dashed") +
  geom_smooth(method = "lm", se = F) +
  geom_point(size = 2) + 
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  labs(x = "Predator mass (mg)", 
       y = "Prey size range (mg)",
       color = "Predator species") +
  #scale_color_manual(values = c("#EEB00C", "#158D8E", "#114C54"),
  #                   labels = pred_labels) +
  theme(axis.text = element_text(size =20),
        axis.title = element_text(size = 25),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 15)) +
  facet_wrap(~sample_str, 
             labeller = labeller(.cols = pred_labels))
range_graph

size %>%
  group_by(sample, sample_str, pred_mass_mg) %>%
  summarise(range = max(mean_prey_mass_mg) - min(mean_prey_mass_mg)) %>%
  filter(range > 0) %>%
  filter(!sample_str %in% c("CEN", "EUB", "LRS")) %>%
  #mutate(sample_str = factor(sample_str, 
  #levels = c("NEO", "PHH", "HEV"))) %>%
  ggplot(aes(x = pred_mass_mg, y = range, color = sample_str)) +
  geom_abline(slope =1, linetype = "dashed") +
  geom_abline(slope = 0.35911) +
  geom_point(size = 2) + 
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  labs(x = "Predator mass (mg)", 
       y = "Prey size range (mg)",
       color = "Predator species") +
  #scale_color_manual(values = c("#EEB00C", "#158D8E", "#114C54"),
  #                   labels = pred_labels) +
  theme(axis.text = element_text(size =20),
        axis.title = element_text(size = 25),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 15)) 

size %>%
  group_by(sample, sample_str, pred_mass_mg) %>%
  summarise(range = max(mean_prey_mass_mg) - min(mean_prey_mass_mg)) %>%
  filter(range > 0) %>%
  filter(!sample_str %in% c("CEN", "EUB", "LRS")) %>%
  #mutate(sample_str = factor(sample_str, 
  #levels = c("NEO", "PHH", "HEV"))) %>%
  ggplot(aes(x = sample_str, y = range, color = sample_str)) +
  geom_boxplot() + 
  theme_bw() +
  scale_y_log10() +
  labs(x = "Predator species", 
       y = "Prey size range (mg)",
       color = "Predator species") +
  #scale_color_manual(values = c("#EEB00C", "#158D8E", "#114C54"),
  #                   labels = pred_labels) +
  theme(axis.text = element_text(size =20),
        axis.title = element_text(size = 25),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 15)) 

