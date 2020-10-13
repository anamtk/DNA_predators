###########################
# Load packages
package.list <- c("here", "tidyverse", "ggplot2", "bipartite", "glmmTMB",
                  "MuMIn", "emmeans", "plotly")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}
#############################

l_nodes <- read.csv(here("Published_webs", "Dryad", "Laigle_etal_2017", "nodes_Laigle.csv"))
l_ints <- read.csv(here("Published_webs", "Dryad", "Laigle_etal_2017", "interactions_Laigle.csv"))

prey <- l_nodes %>%
  dplyr::select(X, mass) %>%
  rename("prey_mass" = "mass",
         "resource" = "X")

pred <- l_nodes %>%
  dplyr::select(X, mass, Web, Poison) %>%
  rename("pred_mass" = "mass",
         "consumer" = "X")

ints <- l_ints %>%
  left_join(prey, by = "resource") %>%
  left_join(pred, by = "consumer") %>%
  filter(interaction > 0) %>%
  mutate(sp_web = ifelse(Web == 1, "web_builder", "non-web_builder")) %>%
  mutate(sp_poison = ifelse(Poison == 1, "poison", "no_poison")) %>%
  filter(!is.na(sp_web)) 

a <- ggplot(ints, aes(x = pred_mass, y = prey_mass, color = sp_web, text = consumer)) +
  geom_point(size = 3) +
  geom_abline(intercept = 0, slope =1) +
  theme_bw() +
  scale_color_manual(values = c("#a6cee3", "#1f78b4")) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~sp_web)

ggplotly(a)

b <- ggplot(ints, aes(x = pred_mass, y = prey_mass, color = sp_poison)) +
  geom_point(size = 3) +
  geom_abline(intercept = 0, slope =1) +
  theme_bw() +
  scale_color_manual(values = c("#a6cee3", "#1f78b4")) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~sp_poison)
b
a
ggplotly(a)

ggplot(ints, aes(x = sp_web, y = prey_mass/pred_mass)) +
  geom_boxplot() +
  geom_hline(yintercept =1, linetype = "dashed") +
  theme_bw() +
  scale_y_log10()





