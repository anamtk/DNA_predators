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

# Load data -----------------------------------------------------------

library(here)

source(here("code", 
            "4_analyses",
            "4_source.R"))

# Body size model selection ---------------------------------------------------------

#is the prey size determined by some combination of predator identity
#and predator size?
m1 <- glmmTMB(mean_prey_log_mass_mg ~ pred_log_mass_mg*sample_str + (1|sample),
              data = size,
              REML = FALSE)

#interaction model: both slope and intercept vary by species
#mass + species model: slope invariant, only intercept varies by species
#mass model: only mass matters
#species model: only species matters
dredge(m1)

#best is mass + species model
m2 <- glmmTMB(mean_prey_log_mass_mg ~ pred_log_mass_mg + sample_str + (1|sample),
              data = size)

# Model Diagnostics -------------------------------------------------------

fit <- simulateResiduals(m2, plot = T)

#in this summary,
#intercept Estimate gives that base intercept for the model
#the intercept for each species is the sum of that and the 
#species intercept
summary(m2)
#the power law relationship is given by the estimate of 
#pred_log_mass_mg, which is sublinear with a relationship of
#y = a + x^0.34336

#rsquared
r.squaredGLMM(m2)

confint(m2) #from glmm package

#compare to each other
pairs(emmeans(m2, ~sample_str))

#look at level effects
plot(allEffects(m2))

# Summary Stats and Tables ------------------------------------------------

#how many total preds
size %>%
  distinct(sample) %>%
  summarise(total = n())

#how many interactions per predator individual
size %>%
  group_by(sample) %>%
  summarise(total = n()) %>%
  summarise(mean = mean(total),
            sd = sd(total),
            max = max(total))

# min and max size per predator species
size %>%
  group_by(sample_str) %>% 
  summarise(min = min(pred_mass_mg),
            max = max(pred_mass_mg))

#table of prey families
size %>%
  distinct(Class, Order, Family) %>%
  arrange(Class, Order, Family) %>% 
  gt() %>%
  tab_header(
    title = "Prey families from DNA diet data")

size %>%
  distinct(Class, Order, Family) %>%
  tally() 

#make a df of this data
prey_fams <- size %>%
  distinct(Class, Order, Family) %>%
  arrange(Class, Order, Family)

#export for supplementary table
write.csv(prey_fams, here("figures", 
                          "supp", 
                          "tables",
                          "prey_families.csv"))


# Model selection table

a <- dredge(m1)

a1 <- a[,c(3:9)]

a1 %>% 
  rename("log10 Predator mass" = "cond(pred_log_mass_mg)",
         "Predator species" = "cond(sample_str)",
         "log10 Predator mass*Predator species" = "cond(pred_log_mass_mg:sample_str)") %>% 
  gt() %>% 
  fmt_number(
    columns = vars("log10 Predator mass", "logLik", "AICc", "delta"),
    decimals = 2) %>% 
  tab_header(
    title = "Model selection of predator-prey size linear model") 

