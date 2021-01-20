##########################
# 1. What are PPMR slopes at individual vs. species level -----
# Ana Miller-ter Kuile
# October 8, 2020
###########################

#basically the findings i want to write up 
#(that i think are in there) are “larger predator individuals 
#within/across species do not eat larger prey than smaller 
#individuals and the ratio of prey size to predator size for 
#this dataset follows similar relationships to other 
#published datasets (using a huge dataset I just found on the internet)”


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

pal <- read.csv(here("data", "outputs", "8_final_dataset",
                     "pred_prey_sizes_tp_DNAinteractions.csv"))

size <- pal %>%
  dplyr::select(-X, -X.x, -X.y, -reads) %>%
  mutate(pred_mass_mg = exp(pred_log_mass_mg)) %>%
  dplyr::select(sample_str, pred_mass_mg, mean_prey_mass_mg, min_prey_mass_mg) %>%
  mutate(source = "AMtK",
         foodweb.name = "pal")

brose <- read.csv(here("data", "brose_2019", "data.csv"))
colnames(brose)
#separate only terrestrial aboveground inverts within the 
#sizes of the samples I have
brose_dat <- brose %>%
  filter(con.mass.mean.g. < 1 & con.mass.mean.g. > 0.0002299683) %>%
  filter(interaction.type %in% c("predacious")) %>%
  filter(ecosystem.type %in% c("terrestrial aboveground", "terrestrial belowground")) %>%
  filter(con.metabolic.type == "invertebrate") %>%
  filter(link.citation != "Lafferty et al. (2006)") %>%
  dplyr::select(con.taxonomy, 
                link.citation,
                foodweb.name,
                con.mass.mean.g., 
                res.mass.mean.g., 
                res.mass.min.g.) %>%
  mutate(pred_mass_mg = con.mass.mean.g.*1000,
         mean_prey_mass_mg = res.mass.mean.g. * 1000,
         min_prey_mass_mg = res.mass.min.g. * 1000) %>%
  rename("sample_str" = "con.taxonomy",
         "source" = "link.citation") %>%
  dplyr::select(sample_str, 
                pred_mass_mg, 
                mean_prey_mass_mg, 
                min_prey_mass_mg,
                source,
                foodweb.name)

data <- size %>%
  bind_rows(brose_dat)

brose_dat %>%
  group_by(source) %>%
  distinct(sample_str) %>%
  tally()

# Size ratio visualizations -----------------------------------------------
data <- data %>%
  mutate(source_alpha = ifelse(source == "AMtK", 2, 1))
alphas <- c("AMtK" = 1, "Cattin Blandenier (2004)" = 0.6, "Piechnik et al. (2008), Simberloff and Wilson (1969)" = 0.6)
str(data)
data %>%
  mutate(type = ifelse(source == "AMtK", "DNA", "literature")) %>%
  #mutate(source = factor(source, levels = c("AMtK",
  #                                          "Cattin Blandenier (2004)",
  #                                          "Piechnik et al. (2008), Simberloff and Wilson (1969)"
  #                                          ))) %>%
  ggplot(aes(x = pred_mass_mg, y = mean_prey_mass_mg, color = type)) +
  geom_abline(slope = 1, linetype = "dashed") +
  geom_point(aes(shape = type), size = 3) +
  geom_smooth(method = "lm", se = F) +
  scale_shape_manual(values = c(16, 1)) +
  scale_color_manual(values = c("#000000", "#d9d9d9")) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() 

ggplot(brose_dat, aes(x = pred_mass_mg, y = mean_prey_mass_mg)) +
  geom_abline(slope = 1, linetype = "dashed", size =1) +
  geom_point(size = 3, color = "#d9d9d9", shape = 1) +
  geom_smooth(method = "lm", se = F, color = "black", size =1) +
  geom_point(data = size, aes(x = pred_mass_mg, y = mean_prey_mass_mg), 
             color = "#252525", size = 3) +
  geom_smooth(data = size, aes(x = pred_mass_mg, y = mean_prey_mass_mg), 
              method = "lm", se = F, color= "black", size =1) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  theme(axis.text = element_text(size =25), axis.title = element_text(size =30)) +
  labs(x = "Predator mass (mg)", y = "Prey mass (mg)")

a
data %>%
  mutate(source = factor(source, levels = c("Cattin Blandenier (2004)",
                                            "Piechnik et al. (2008), Simberloff and Wilson (1969)",
                                            "pal"))) %>%
  ggplot(aes(x = pred_mass_mg, y = mean_prey_mass_mg/pred_mass_mg, color = foodweb.name)) +
  geom_abline(slope = 1, linetype = "dashed") +
  geom_point(size = 3, alpha = 0.6) +
  geom_smooth(method = "lm", se = F) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  facet_wrap(~source)

# ratio model -------------------------------------------------------------
ratios <- data %>%
  mutate(ratio= mean_prey_mass_mg/pred_mass_mg,
         type = ifelse(source == "AMtK", "individual", "species"),
         log_ratio = log(ratio),
         log_pred_mass = log(pred_mass_mg))

m1 <- glmmTMB(log_ratio ~ log_pred_mass*type + (1|source) + (1|foodweb.name),
              data = ratios)

summary(m1)
dredge(m1)

m2 <- glmmTMB(log_ratio ~ log_pred_mass + (1|source) + (1|foodweb.name),
              data = ratios)

plot(allEffects(m2))

fit <- simulateResiduals(m2, plot = T)
summary(m2)




