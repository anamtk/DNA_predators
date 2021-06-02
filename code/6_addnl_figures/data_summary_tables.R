# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", "gt")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Import data -------------------------------------------------------------

data <- read.csv(here("data", 
                      "outputs", 
                      "2i_final_dataset",
                      "pred_prey_sizes_DNAinteractions.csv"))

data2 <- read.csv(here("data", 
                      "outputs", 
                      "2e_rarefied_taxonomic_sort",
                      "fam_prey_DNA_conservative.csv"))

data2$sample <- str_sub(data2$sample, end=-2)

data2 <- data2 %>% 
  dplyr::select(sample_str, pred_ID, sample, run) %>% 
  distinct(sample_str, pred_ID, sample, run)

data <- data %>%
  mutate(pred_mass_mg = exp(pred_log_mass_mg)) %>% 
  dplyr::select(sample, sample_str, pred_mass_mg) %>% 
  distinct(sample, sample_str, pred_mass_mg)

stats <- data %>% 
  group_by(sample_str) %>%
  summarise("min_mg" = min(pred_mass_mg),
            "max_mg" = max(pred_mass_mg),
            "mean_mg" = mean(pred_mass_mg))
  
data <- data %>%
  left_join(data2, by = c("sample", "sample_str"))

data %>%
  group_by(pred_ID, run) %>% 
  summarise(Samples = n()) %>%
  dplyr::select(pred_ID, Samples, run) %>% 
  rename("Species" = "pred_ID",
         "Run" = "run") %>% 
  ungroup() %>%
  mutate(Run = as.factor(Run)) %>%  
  arrange(Run) %>% 
  gt() %>% 
  tab_header(
    title = "Species and sample sizes by sequencing run") %>% 
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(
      columns = vars("Species")))
    
data2 %>%
  dplyr::select(sample_str, pred_ID) %>%
  distinct(sample_str, pred_ID) %>% 
  left_join(stats, by = "sample_str") %>%
  dplyr::select(-sample_str) %>% 
  ungroup() %>%
  arrange(mean_mg) %>%
  rename("Species" = "pred_ID",
         "Min Size (mg)" = "min_mg",
         "Max Size (mg)" = "max_mg",
         "Mean Size (mg)" = "mean_mg") %>%
  gt() %>% 
  fmt_number(
    columns = vars("Min Size (mg)", "Max Size (mg)", "Mean Size (mg)"),
    decimals = 1) %>% 
  tab_header(
    title = "Species size statistics") %>% 
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(
      columns = vars("Species")))

