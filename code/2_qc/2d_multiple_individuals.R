# Samples with multiple individuals
# Ana Miller-ter Kuile
# March 2, 2021

# this script examines whether samples which were run with multiple 
# individuals vary in the number of interactions observed on those 
# samples

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 
                  "glmmTMB", "MuMIn",
                  "DHARMa", "gt")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Load data ---------------------------------------------------------------

data <- read.csv(here("data", 
                      "outputs",
                      "3i_final_dataset",
                      "pred_prey_sizes_DNAinteractions.csv"))

meta <- read.csv(here("data",
                      "raw_data",
                      "Sample_metadata.csv"))

meta <- meta %>%
  dplyr::select(Extraction.ID, No..Individuals) %>%
  rename("sample" = "Extraction.ID") %>%
  distinct(sample, No..Individuals)

samples <- read.csv(here("data",
                      "raw_data",
                      "Sample_metadata.csv"))

# Manipulate datasets --------------------------------------------------------

data <- data %>%
  left_join(meta, by = "sample")

counts <- data %>%
  group_by(sample, sample_str, No..Individuals) %>%
  tally(name = "interactions")

# Data summaries ------------------------------------------------------

counts %>%
  filter(No..Individuals > 1) %>%
  ungroup() %>%
  tally()

counts %>%
  filter(No..Individuals == 1) %>%
  ungroup() %>%
  tally()

counts %>%
  ungroup() %>%
  summarise(mean = mean(No..Individuals),
            sd = sd(No..Individuals),
            total = n(),
            se = sd/sqrt(total),
            max = max(No..Individuals))

samples %>%
  group_by(Extraction.ID) %>%
  summarise(sd = sd(Length_mm)) %>%
  filter(!is.na(sd)) %>%
  summarise(mean_sd = mean(sd),
            sd_sd = sd(sd),
            total = n(),
            se_sd = sd_sd/sqrt(total))

samples %>%
  filter(No..Individuals > 1) %>%
  summarise(mean = mean(Length_mm, na.rm = T),
            sd = sd(Length_mm, na.rm = T),
            total = n(),
            se = sd/sqrt(total))

# Model -------------------------------------------------------------------

## AB note: not sure if the model code is supposed to throw these warnings, but this is what I got (may be worth making a note in the code?)
## error message for the glmmTMB() lines:
# Warning message:
#   In (function (start, objective, gradient = NULL, hessian = NULL,  :
#                   NA/NaN function evaluation
## error message for the dredge() lines:
# Warning messages:
#   1: In (function (start, objective, gradient = NULL, hessian = NULL,  :
#                      NA/NaN function evaluation
#                    2: In (function (start, objective, gradient = NULL, hessian = NULL,  :
#                                       NA/NaN function evaluation


m <- glmmTMB(interactions ~ No..Individuals + (1|sample_str),
             data = counts, family = "genpois")

dredge(m)

m2 <- glmmTMB(interactions ~ 1 + (1|sample_str),
              data = counts, family = "genpois")


fit <- simulateResiduals(m2, plot = T)

a <- dredge(m)

a1 <- a[,c(3:8)]

# another warning:
# Warning message:
#   In `[.model.selection`(a, , c(3:8)) :
#   cannot recalculate "weights" on an incomplete object

a1 %>% 
  rename("No..Individuals" = "cond(No..Individuals)",) %>% 
  gt() %>% 
  fmt_number(
    columns = vars("No..Individuals", "logLik", "AICc", "delta"),
    decimals = 2) %>% 
  tab_header(
    title = "Model selection of number of individuals per sample model") 

# Visualize ---------------------------------------------------------------

counts %>%
  mutate(No..Individuals = factor(No..Individuals)) %>%
ggplot(aes(x = No..Individuals, y = interactions)) +
  geom_boxplot() +
  theme_bw()

counts %>%
  ggplot(aes(x = No..Individuals)) +
  geom_histogram() +
  theme_bw()


