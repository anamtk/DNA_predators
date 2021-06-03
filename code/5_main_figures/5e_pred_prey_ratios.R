# Ratio figure
# Ana Miller-ter Kuile
# June 2, 2021


# Load data -----------------------------------------------------------

library(here)

source(here("code", 
            "4_analyses",
            "4_source.R"))


# Universal figure design -------------------------------------------------

pal_kelp <- cal_palette("kelp1", n = 9, type = "continuous")

pred_labels <- c("CEN" = "Mecistocephalus sp.", "EUB" = "E. annulipes", 
                 "HEV" = "H. venatoria", "LRS" = "Opopaea sp.", 
                 "NEO" = "N. theisi", "PAN" = "P. flavescens",
                 "PHH" = "P. holdhausi", "SCY" = "S. longipes",
                 "SME" = "S. pallidus")


# Figure ------------------------------------------------------------------


ratios %>%
  mutate(sample_str = fct_reorder(sample_str, pred_mass_mg, .fun='mean')) %>%
  ggplot(aes(x = sample_str, y = ratio, color = sample_str, fill= sample_str)) +
  geom_boxplot(size = 0.75, alpha = 0.6) +
  geom_jitter(width = 0.25, height = 0, shape = 1) +
  theme_bw() +
  scale_color_manual(
    values = c("#C70000", "#EA7700", "#EEB00C", "#89742F", "#114C54",
               "#C68A2C",
               "#496C3C", "#158D8E", "#067D8D")) +
  scale_fill_manual(
    values = c("#C70000", "#EA7700", "#EEB00C", "#89742F", "#114C54",
               "#C68A2C",
               "#496C3C", "#158D8E", "#067D8D")) +
  scale_y_log10(breaks = c(0.01, 1, 100, 10000)) +
  labs(x = "Predator class", y = "Predator:prey mass ratio") +
  theme(axis.text = element_text(size =20),
        axis.title = element_text(size = 25)) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 1) +
  theme(axis.text = element_text(size =20),
        axis.title = element_text(size = 25),
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 15)) +
  facet_grid(.~pred_class, scales = "free_x", space = "free")

