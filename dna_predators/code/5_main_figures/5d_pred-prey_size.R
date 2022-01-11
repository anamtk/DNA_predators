# Main Document Figures
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

# Figures----------------------------------------------------------------


#Does predator identity or size determine prey size?
(size_graph_col <- size %>%
   mutate(sample_str = factor(sample_str, levels = c("LRS", "SCY", "NEO", "CEN",
                                                     "SME", "EUB", "PHH", "PAN",
                                                     "HEV"))) %>%
   ggplot(aes(x = pred_mass_mg, y = mean_prey_mass_mg, color = sample_str)) +
   geom_abline(slope = 1, linetype = "dashed", size = 0.75) +
   geom_point(size = 3) +
   geom_abline(slope = 0.31, size = 1) +
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
    annotate(geom = "text", x = 1, y = 400, label = "+", size = 8) +
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
    annotate(geom = "text", x = 1, y = 400, label = "+", size = 8) +
    annotate(geom = "text", x = 4, y = 400, label = "-", size = 8) +
    annotate(geom = "text", x = 8, y = 400, label = "-", size = 8) +
    annotate(geom = "text", x = 6, y = 400, label = "+", size = 8))

size_graph_col / species_graph +
  plot_layout(guides = 'collect')

size_graph_no_leg + species_graph_nox


