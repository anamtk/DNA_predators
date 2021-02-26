install.packages("ape")
library(ape)
data(carnivora)
frm <- ~SuperFamily/Family/Genus/Species
tr <- as.phylo(frm, data = carnivora)
plot(tr)
Nnode(tr)
## compare with:
Nnode(as.phylo(frm, data = carnivora, collapse = FALSE))

pal_prey <- read.csv(here("Drafts", 
                          "Figures", 
                          "Supp", 
                          "prey_families.csv"))


pal_prey <- pal_prey %>%
  mutate(Class = as.factor(Class),
         Order = as.factor(Order),
         Family = as.factor(Family),
         Phylum = as.factor(Phylum),
         Subphylum = as.factor(Subphylum)) 
frm <- ~Phylum/Subphylum/Class/Order/Family
tr <- as.phylo(frm, data = pal_prey, use.labels = T)
plot(tr)
Nnode(tr)

ggtree(tr, branch.length = F) +
  geom_tiplab() +
  xlim(0, 5)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")

library(ggtree)

ggtree(tr) 
?tiplabels

pal_prey <- read.tree(here("2_taxonomic_assignment",
                           "taxonomies",
                           "DADA_MAY.tre"))


ggtree(pal_prey, branch.length = F) +
  geom_tiplab() +
  xlim(0, 40)


library(igraph)
g <- make_full_bipartite_graph(2, 3)
plot(g)
