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
                  "emmeans", "gt")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Import data -------------------------------------------------------------

data <- read.csv(here("data", 
                      "outputs",  
                      "3g_final_dataset", 
                      "pred_prey_sizes_DNAinteractions.csv"))

size <- data %>%
  dplyr::select(-X, -reads) %>%
  mutate(pred_mass_mg = 10^(pred_log_mass_mg))

# How many preds? -----------------------------------------------
#how many total preds
size %>%
  distinct(sample) %>%
  summarise(total = n())

#how many interactions per predator individual
size %>%
  group_by(sample) %>%
  summarise(total = n()) %>%
  summarise(mean = mean(total),
            sd = sd(total))

# min and max size per predator species
size %>%
  group_by(sample_str) %>% 
  summarise(min = min(pred_mass_mg),
            max = max(pred_mass_mg))

#table of prey families
size %>%
  distinct(Class, Order, Family) %>%
  arrange(Class, Order, Family) %>% 
  gt() %>%
  tab_header(
    title = "Prey families from DNA diet data")

size %>%
  distinct(Class, Order, Family) %>%
  tally() 

#make a df of this data
prey_fams <- size %>%
  distinct(Class, Order, Family) %>%
  arrange(Class, Order, Family)

#export for supplementary table
write.csv(prey_fams, here("Drafts", 
                          "Figures", 
                          "Supp", 
                          "prey_families.csv"))

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
a <- dredge(m1)

a1 <- a[,c(3:9)]

a1 %>% 
  rename("log10 Predator mass" = "cond(pred_log_mass_mg)",
         "Predator species" = "cond(sample_str)",
         "log10 Predator mass*Predator species" = "cond(pred_log_mass_mg:sample_str)") %>% 
  gt() %>% 
  fmt_number(
    columns = vars("log10 Predator mass", "logLik", "AICc", "delta"),
    decimals = 2) %>% 
  tab_header(
    title = "Model selection of predator-prey size linear model") 

#best is mass + species model
m2 <- glmmTMB(mean_prey_log_mass_mg ~ pred_log_mass_mg + sample_str + (1|sample),
              data = size)

# Model Diagnostics -------------------------------------------------------

fit <- simulateResiduals(m2, plot = T)

#in this summary,
#intercept Estimate gives that base intercept for the model
#the intercept for each species is the sum of that and the 
#species intercept
summary(m2)
#the power law relationship is given by the estimate of 
#pred_log_mass_mg, which is sublinear with a relationship of
#y = a + x^0.34336

# Figures ------------------------------------------------------------------

pal_kelp <- cal_palette("kelp1", n = 9, type = "continuous")

pred_labels <- c("CEN" = "Geophilomorpha sp.", "EUB" = "E. annulipes", 
                 "HEV" = "H. venatoria", "LRS" = "Oonopidae sp.", 
                 "NEO" = "N. theisi", "PAN" = "P. flavescens",
                 "PHH" = "P. holdhausi", "SCY" = "S. longipes",
                 "SME" = "S. pallidus")

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
(size_graph_col <- size %>%
  mutate(sample_str = factor(sample_str, levels = c("LRS", "SCY", "NEO", "CEN",
                                                    "SME", "EUB", "PHH", "PAN",
                                                    "HEV"))) %>%
  ggplot(aes(x = pred_mass_mg, y = mean_prey_mass_mg, color = sample_str)) +
  geom_abline(slope = 1, linetype = "dashed", size = 0.75) +
  geom_point(size = 3) +
  geom_abline(slope = 0.34, size = 1) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Predator mass (mg)", 
       y = "Prey mass (mg)",
       color = "Predator species") +
  scale_color_manual(values = pal_kelp,
                     labels = pred_labels) +
  theme_bw() +
  theme(axis.text = element_text(size =20),
        axis.title = element_text(size = 25),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)))

(size_graph_no_leg <- size %>%
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
          axis.title = element_text(size = 25),
          legend.position = "none",
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15)))
#sorted by increased average predator size, then showing prey size
(species_graph <- size %>%
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
        axis.title = element_text(size = 25),
        legend.position = "none") +
  annotate(geom = "text", x = 4, y = 400, label = "-", size = 8) +
  annotate(geom = "text", x = 8, y = 400, label = "-", size = 8) +
  annotate(geom = "text", x = 6, y = 400, label = "+", size = 8))

(species_graph_nox <- size %>%
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
    theme(axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.title = element_text(size = 25),
          axis.title.y = element_blank(),
          legend.position = "none") +
    annotate(geom = "text", x = 4, y = 400, label = "-", size = 8) +
    annotate(geom = "text", x = 8, y = 400, label = "-", size = 8) +
    annotate(geom = "text", x = 6, y = 400, label = "+", size = 8))

size_graph_col / species_graph +
  plot_layout(guides = 'collect')

size_graph_no_leg + species_graph_nox


(pred_size_hist <- size %>%
  distinct(sample, sample_str, pred_mass_mg) %>%
  mutate(sample_str = fct_reorder(sample_str, pred_mass_mg, .fun='mean')) %>%
  ggplot(aes(x = pred_mass_mg, fill = sample_str)) +
  geom_histogram(color = "black", alpha = 0.85) +
  theme_bw() + 
  scale_x_log10() +
  labs(x = "Predator mass (mg)", y = "Individuals",
       colour = "Predator species") +
  theme(axis.text = element_text(size =20),
        axis.title = element_text(size = 25),
        axis.text.x = element_text(angle = 45, hjust =1)) +
  facet_wrap(~sample_str, labeller = labeller(sample_str = pred_labels)) +
  scale_fill_manual(values = pal_kelp,
                    labels = pred_labels) +
  theme(legend.position = "none",
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 15)))


x <- c(1:100)
y3 <- 0.34*(x)
plot(y3 ~ x)
x4 <- 10^x
y4 <- 10^y3
plot(y4 ~ x4)

