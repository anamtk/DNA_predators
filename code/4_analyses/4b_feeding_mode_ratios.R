##########################
# 2. Are prey:predator ratios explained by feeding mode? -----
# Ana Miller-ter Kuile
# October 29, 2020
###########################

# this script analyzes whether a set of predator traits including 
# hunting mode, web use and venom use, influence prey size

# Import data -------------------------------------------------------------

library(here)

source(here("code", 
            "4_analyses",
            "4_source.R"))


# Models ------------------------------------------------------------------

m_webs <- glmmTMB(log_ratio ~ webs + (1|sample) + (1|sample_str),
                data = ratios)

m_venom <- glmmTMB(log_ratio ~ venom + (1|sample) + (1|sample_str),
               data = ratios)

m_class <- glmmTMB(log_ratio ~ pred_class + (1|sample) + (1|sample_str),
                   data = ratios)

m_null <- glmmTMB(log_ratio ~ 1 + (1|sample) + (1|sample_str),
                  data = ratios)

AICc(m_webs, m_class, m_venom, m_null)


# Model Diagnostics -------------------------------------------------------

simulateResiduals(m_class, plot =T)

summary(m_class)

plot(allEffects(m_class))

pairs(emmeans(m_class, "pred_class"))


# Summaries and tables ----------------------------------------------------

# traits sample sizes table
ratios %>%
  distinct(sample, sample_str, hunting_mode, venom, webs, pred_class) %>%
  group_by(sample_str, hunting_mode, venom, webs, pred_class) %>%
  tally(name = "samples") %>%
  mutate(
    species = case_when(
      sample_str == "CEN" ~ "Mecistocephalus sp.",
      sample_str == "EUB" ~ "E. annulipes",
      sample_str == "HEV" ~ "H. venatoria",
      sample_str == "LRS" ~ "Opopaea sp.",
      sample_str == "NEO" ~ "N. theisi",
      sample_str == "PAN" ~ "P. flavescens",
      sample_str == "PHH" ~ "P. holdhausi", 
      sample_str == "SCY" ~ "S. longipes",
      sample_str == "SME" ~ "S. pallidus"
    )) %>%
  ungroup() %>%
  dplyr::select(species, venom, webs, pred_class, samples) %>%
  rename(Class = pred_class) %>%
  left_join(predator_traits, by = "species") %>%
  gt() %>%
  tab_header(
    title = "Number of samples and interactions per species and traits") 

#model summary table

a <- AICc(m_class, m_null, m_venom, m_webs)

a %>% 
  mutate(delta = AICc - 850.5784) %>% 
  arrange(delta) %>%
  rownames_to_column(var = "model") %>%
  gt() %>% 
  fmt_number(
    columns = vars("df", "AICc", "delta"),
    decimals = 2) %>% 
  tab_header(
    title = "Model selection of predator trait models") 

# prey frequency table
eaten <- size %>%
  group_by(sample_str, Order, Family) %>%
  summarise("Frequency (No. Individuals)" = n())

write.csv(eaten, here("figures", "supp", "tables", "prey_frequency.csv"))
