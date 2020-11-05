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
                  "calecopal", "patchwork")

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
size %>%
  summarise(min = min(pred_mass_mg),
            max = max(pred_mass_mg))

size %>%
  distinct(sample) %>%
  summarise(total = n())
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

pal_kelp <- cal_palette("kelp1", n = 9, type = "continuous")

(mean_predprey_plot <- plot(me, ci = FALSE, line.size = 1) +
  geom_abline(slope = 1, linetype = "dashed", size = 0.75) +
  annotate("text", x = 7, y = -3, label = expression(paste("y = ", x^0.41))) +
  annotate("text", x = 7, y = -3.5, label = expression(paste(cR^2, " = 0.35"))) +
  annotate("text", x = 7, y = -4, label = expression(paste(mR^2, " = 0.30"))) +
  scale_color_manual(values = pal_kelp) +
  labs(x = "Predator mass (log(mg))", y = "Prey mass (log(mg))") +
  theme_bw() +
  theme(axis.text = element_text(size = 20), 
        axis.title = element_text(size = 25),
        title = element_blank()))

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

min_predprey_plot <- plot(me2, ci = FALSE, line.size = 1) +
  geom_abline(slope = 1, linetype = "dashed", size = 0.75) +
  annotate("text", x = 7, y = -3, label = expression(paste("y = ", x^0.26))) +
  annotate("text", x = 7, y = -3.5, label = expression(paste(cR^2, " = 0.15"))) +
  annotate("text", x = 7, y = -4, label = expression(paste(mR^2, " = 0.15"))) +
  scale_color_manual(values = pal_kelp) +
  labs(x = "Predator mass (log(mg))", y = "Prey mass (log(mg))") +
  theme_bw() +
  theme(axis.text = element_text(size = 20), 
        axis.title = element_text(size = 25),
        title = element_blank())

#in this summary,
#intercept Estimate gives that base intercept for the model
#the intercept for each species is the sum of that and the 
#species intercept
#the power law relationship is given by the estimate of 
#pred_log_mass_mg, which is sublinear with a relationship of
#y = a + x^0.26466
summary(m4)
r.squaredGLMM(m4)


# Figures ------------------------------------------------------------------
#predator size distribution
(pred_size <- size %>%
  distinct(sample, pred_log_mass_mg) %>%
ggplot(aes(x = pred_log_mass_mg)) +
  geom_histogram(bins = 50, alpha = 0.85, color = "black") +
  #scale_y_log10() +
  #geom_density(fill = "#BE8333", alpha = 0.85) +
  xlim(-2.5, 7.75) +
  theme_bw() +
  labs(y = "Individuals") +
  theme(axis.text = element_text(size =20),
        axis.title = element_text(size = 25),
        axis.title.x = element_blank()) +
  scale_y_reverse())

prey_mean_sz <- size %>%
  #distinct(Family, mean_prey_log_mass_mg) %>%
ggplot(aes(x = mean_prey_log_mass_mg)) +
  geom_histogram(alpha = 0.85, color = "black") +
  #geom_density(fill = "#114C54", alpha = 0.85) +
  #scale_y_log10() +
  scale_x_log10() +
  xlim(-4, 6) +
  theme_bw() +
  labs(y = "Interactions") +
  theme(axis.text = element_text(size =20),
        axis.title = element_text(size = 25),
        axis.title.y = element_blank()) +
  coord_flip()

prey_min_sz <- size %>%
  #distinct(Family, min_prey_log_mass_mg) %>%
ggplot(aes(x = min_prey_mass_mg)) +
  geom_histogram(alpha = 0.85, color = "black") +
  #scale_y_log10() +
  scale_x_log10() +
  #geom_density(fill = "#114C54", alpha = 0.85) +
  theme_bw() +
  labs(x = "Prey mass (mg)", y = "Interactions") +
  theme(axis.text = element_text(size =20),
        axis.title = element_text(size = 25)) +
  coord_flip()

mean_size_graphs <- (mean_predprey_plot + prey_mean_sz + pred_size + plot_spacer() +
  plot_layout(widths = c(3,1), heights = c(3,1)))

min_size_graphs <- (min_predprey_plot + prey_min_sz + pred_size + plot_spacer() +
                       plot_layout(widths = c(3,1), heights = c(3,1)))


#pred size distribution on its own
(pred_size2 <- size %>%
    distinct(sample, pred_log_mass_mg) %>%
    ggplot(aes(x = pred_log_mass_mg)) +
    geom_histogram(bins = 50, alpha = 0.85, color = "black") +
    #scale_y_log10() +
    xlim(-8, 7) +
    #geom_density(fill = "#BE8333", alpha = 0.85) +
    theme_bw() +
    labs(x = "Predator mass (log(mg))", y = "Individuals") +
    theme(axis.text = element_text(size =20),
          axis.title = element_text(size = 25)))

#-1.469814
#6.834101
#-7.36504

#prey size distribution on its own
(prey_mean_size2 <- size %>%
    distinct(Family, mean_prey_log_mass_mg) %>%
    ggplot(aes(x = mean_prey_log_mass_mg)) +
    geom_histogram(bins = 50, alpha = 0.85, color = "black") +
    #scale_y_log10() +
    #geom_density(fill = "#BE8333", alpha = 0.85) +
    theme_bw() +
    xlim(-8, 7) +
    labs(x = "Mean prey mass (log(mg))", y = "Prey families") +
    theme(axis.text = element_text(size =20),
          axis.title = element_text(size = 25)))

(prey_min_size2 <- size %>%
    distinct(Family, min_prey_log_mass_mg) %>%
    ggplot(aes(x = min_prey_log_mass_mg)) +
    geom_histogram(bins = 50, alpha = 0.85, color = "black") +
    #scale_y_log10() +
    xlim(-8, 7) +
    #geom_density(fill = "#BE8333", alpha = 0.85) +
    theme_bw() +
    labs(x = "Min prey mass (log(mg))", y = "Prey families") +
    theme(axis.text = element_text(size =20),
          axis.title = element_text(size = 25)))


