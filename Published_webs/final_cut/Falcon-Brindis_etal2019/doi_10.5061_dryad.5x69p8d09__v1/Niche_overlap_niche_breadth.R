library(spaa)

datos <- read.csv("Falcon-Brindis_general_table.csv") # Data conteining the spider abundance captured by each wasp species

datos <- datos[, -1]

data("datasample")
head(datasample)
  
niche.overlap(datos[,1:3], method = "pianka") # Could be a pairwise comparised.

niche.overlap(datos, method = "morisita") # Appropriate for cound data.

niche.overlap(datos, method = "pianka")

niche.overlap(datos, method = "schoener") # Multiply x 100 to obtain the percentage.

# Nothe that this result is not standardized, you should apply standarization. Follow 52.	Levins R. Evolution in changing environments: Some theoretical explorations. Princeton: Princeton University Press; 1968.

# Estimate niche breadth

niche.width(datos, method = "levins") # Abundant species
niche.width(datos, method = "shannon") # For rare species (uncommon)

