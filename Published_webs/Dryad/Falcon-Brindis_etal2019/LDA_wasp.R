
datos <- read.csv("Falcon-Brindis_Wasp_body_measurements.csv")
datos2 = datos[-1]

# Linnear Discriminant Analysis (LDA)
library(psych)
pairs.panels(datos[2:5],
             gap = 0,
             bg = c("red", "green", "blue") [datos$wasp],
             pch = 21)
library(tidyverse)
library(caret)
theme_set(theme_classic())
# Split the data into training (80%) and test set (20%)
set.seed(123)
training.samples <- datos2$wasp %>%
  createDataPartition(p = 0.8, list = FALSE)
train.data <- datos2[training.samples, ]
test.data <- datos2[-training.samples, ]

# Estimate preprocessing parameters
preproc.param <- train.data %>% 
  preProcess(method = c("center", "scale"))
# Transform the data using the estimated parameters
train.transformed <- preproc.param %>% predict(train.data)
test.transformed <- preproc.param %>% predict(test.data)

library(MASS)
# Fit the model
model <- lda(wasp~., data = train.transformed)
# Make predictions
predictions <- model %>% predict(test.transformed)
# Model accuracy
mean(predictions$class==test.transformed$wasp)

model

plot(model)
predictions <- model %>% predict(test.transformed)
names(predictions)

# Predicted classes
head(predictions$class, 6)
# Predicted probabilities of class memebership.
head(predictions$posterior, 6) 
# Linear discriminants
head(predictions$x, 3) 

lda.data <- cbind(train.transformed, predict(model)$x) # Gráfico
ggplot(lda.data, aes(LD1, LD2)) +
  geom_point(aes(color = wasp))