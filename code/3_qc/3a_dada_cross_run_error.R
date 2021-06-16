#Dada2 Error Rate Check Across Runs
#Ana Miller-ter Kuile
#February 3, 2020
#This was drawn from both Happy Belly Informatics and the github issue about 
#combining runs for dada2
#https://astrobiomike.github.io/amplicon/dada2_workflow_ex
#https://github.com/benjjneb/dada2/issues/257

# Load packages -----------------------------------------------------------

library(here)
library(tidyverse)

# Import data -------------------------------------------------------------

error_f <- read.csv(here("data", 
                         "outputs", 
                         "3a_crossrun_error_checking",
                         "forward_errors.csv"))

error_r <- read.csv(here("data", 
                         "outputs", 
                         "3a_crossrun_error_checking",
                         "reverse_errors.csv"))


# Visualize ---------------------------------------------------------------

(forward_graph <- ggplot(error_f, 
                         aes(x = Qual, 
                             y = count, 
                             color = run)) +
   geom_point() + 
   scale_y_log10() + 
   facet_wrap(~Transition)) + 
  theme_bw()

(reverse_graph <- ggplot(error_r, 
                         aes(x = Qual, 
                             y = count, 
                             color = run)) +
    geom_point() + 
    scale_y_log10() + 
    facet_wrap(~Transition)) + 
  theme_bw() 










