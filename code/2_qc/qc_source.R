##################
# 2_qc source script
# Ana Miller-ter Kuile
# March 25, 2021
##################

# This script is the source script for all of the 2_qc scripts
# anything repeated across those scripts is contained in this
# script and called into each script 


# Load packages -----------------------------------------------------------

package.list <- c("here", 
                  "tidyverse", 
                  "gt", 
                  "vegan",
                  "glmmTMB",
                  "DHARMa",
                  "effects",
                  "emmeans",
                  "MuMIn",
                  "lattice",
                  "performance")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


