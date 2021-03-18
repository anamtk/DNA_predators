
library(ape)
data(carnivora)
frm <- ~SuperFamily/Family/Genus/Species
tr <- as.phylo(frm, data = carnivora)
plot(tr)
Nnode(tr)
## compare with:
Nnode(as.phylo(frm, data = carnivora, collapse = FALSE))

pal_prey <- read.csv(here("figures", 
                          "supp", 
                          "tables", 
                          "prey_families.csv"))


pal_prey <- pal_prey %>%
  mutate(Class = as.factor(Class),
         Order = as.factor(Order),
         Family = as.factor(Family)) 
frm <- ~Class/Order/Family
tr <- as.phylo(frm, data = pal_prey, use.labels = T)
plot(tr)
Nnode(tr)