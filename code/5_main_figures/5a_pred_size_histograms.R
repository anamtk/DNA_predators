# Predator size histograms
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


# Figure  ----------------------------------------------------------------

#Predator size histogram Figure
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

