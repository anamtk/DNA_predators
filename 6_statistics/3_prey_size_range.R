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

# SD Predator size within an individual -----------------------------------
SD_prey <- size %>%
  group_by(sample, sample_str, pred_log_mass_mg, pred_mass_mg) %>%
  summarise(SD_mean = sd(mean_prey_mass_mg)) %>%
  filter(!is.na(SD_mean)) %>%
  filter(sample_str %in% c("HEV", "NEO", "PHH"))
hist(SD_prey$SD_mean)
hist(SD_prey$pred_mass_mg)

sample_sizes <- SD_prey %>%
  distinct(sample, sample_str) %>%
  group_by(sample_str) %>%
  summarise(n = n())

SD_prey <- SD_prey %>%
  left_join(sample_sizes, by = "sample_str")

range_prey <- size %>%
  group_by(sample, sample_str, pred_log_mass_mg, pred_mass_mg) %>%
  summarise(range = max(mean_prey_mass_mg) - min(mean_prey_mass_mg)) %>%
  filter(range > 0) %>%
  filter(sample_str %in% c("HEV", "NEO", "PHH"))
hist(range_prey$range)
#some species now have very few individuals, so we should remove
#those here as well.
size %>%
  group_by(sample, sample_str, pred_mass_mg) %>%
  summarise(SD_mean = sd(mean_prey_mass_mg)) %>%
  filter(sample_str %in% c("HEV", "NEO", "PHH")) %>%
  ggplot(aes(x = pred_mass_mg, y = SD_mean, color = sample_str)) +
  geom_point(size = 2) +
  theme_bw() +
  facet_wrap(~sample_str, scales = "free")

# SD size model -----------------------------------------------------------
HEV <- range_prey %>%
  filter(sample_str == "HEV") 

NEO <- range_prey %>%
  filter(sample_str == "NEO")

PHH <- range_prey %>%
  filter(sample_str == "PHH")

m_hev <- glm(range ~ pred_mass_mg,
          data = HEV,
          na.action = "na.fail")

dredge(m_hev)

fit <- simulateResiduals(m_hev, plot = T)

plot(allEffects(m_hev))
summary(m_hev)

m_neo <- glm(range ~ pred_mass_mg,
             data = NEO,
             na.action = "na.fail")

dredge(m_neo)

fit <- simulateResiduals(m_neo, plot = T)

plot(allEffects(m_neo))
summary(m_neo)
hist(log(PHH$range))
m_phh <- glm(log(range) ~ pred_mass_mg,
             data = PHH,
             na.action = "na.fail")

dredge(m_phh)

m_phh2 <- glm(log(range) ~ 1,
             data = PHH,
             na.action = "na.fail")

fit <- simulateResiduals(m_phh2, plot = T)

plot(allEffects(m_phh))
summary(m_phh)

m <- glm(log(range) ~ pred_log_mass_mg*sample_str,
             data = range_prey,
             na.action = "na.fail")

dredge(m)

m_phh2 <- glm(log(range) ~ 1,
              data = PHH,
              na.action = "na.fail")

fit <- simulateResiduals(m, plot = T)

plot(allEffects(m))
summary(m)
r.squaredGLMM(m)

# Visualizations ----------------------------------------------------------
x <- 1:100
y <- x^(-.8511)
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
  filter(sample_str %in% c("NEO", "PHH", "HEV")) %>%
  mutate(sample_str = factor(sample_str, 
                             levels = c("NEO", "PHH", "HEV"))) %>%
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
  scale_color_manual(values = c("#EEB00C", "#158D8E", "#114C54"),
                     labels = pred_labels) +
  theme(axis.text = element_text(size =20),
        axis.title = element_text(size = 25),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 15)) +
  facet_wrap(~sample_str, 
             labeller = labeller(.cols = pred_labels))
