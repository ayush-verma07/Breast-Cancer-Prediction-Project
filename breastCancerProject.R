#All of the libraries we need
options(digits = 3)
library(matrixStats)
library(tidyverse)
library(caret)
library(dslabs)
data(brca)
library(gam)

#Section 1: Outlining all dimensions and properties
#number of samples
num_samples <- nrow(brca$x)
#number of predictors
num_predictors <- ncol(brca$x)
#proportion of malignant
prop_malignant <- mean(brca$y == "M")
#column number highest mean
highest_mean_col <- which.max(colMeans(brca$x))
#column number lowest standard deviation
lowest_sd_col <- which.min(colSds(brca$x))


#scaled the matrix using sweep
scaled_x <- sweep(brca$x, 2, colMeans(brca$x), "-")
scaled_x <- sweep(scaled_x, 2, colSds(brca$x), "/")

#standard deviation of the first column after scaling
sd_first_col <- sd(scaled_x[, 1])
#median value of the first column after scaling
median_first_col <- median(scaled_x[, 1])

#Section 2: do PCA
#scaled matrix
scaled_x <- sweep(brca$x, 2, colMeans(brca$x), "-")
scaled_x <- sweep(scaled_x, 2, colSds(brca$x), "/")

#do PCA on matrix
pca_result <- prcomp(scaled_x, center = TRUE, scale. = TRUE)

#proportion of variance explained by the first principal component
prop_variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
prop_variance_explained_pc1 <- prop_variance_explained[1]

#number of principal components required to explain at least 90% of the variance
cumulative_variance_explained <- cumsum(prop_variance_explained)
num_pcs_90_var <- which(cumulative_variance_explained >= 0.90)[1]

#plot the first two principal components with color representing tumor type
pca_data <- data.frame(pca_result$x[, 1:2], Type = brca$y)
ggplot(pca_data, aes(x = PC1, y = PC2, color = Type)) +
  geom_point() +
  labs(title = "PCA of BRCA Dataset",
       x = "Principal Component 1",
       y = "Principal Component 2")

#create a boxplot of the first 10 PCs grouped by tumor type
pca_data_10pcs <- data.frame(pca_result$x[, 1:10], Type = brca$y)
pca_data_10pcs_melted <- pca_data_10pcs %>%
  pivot_longer(cols = -Type, names_to = "PC", values_to = "Value")

ggplot(pca_data_10pcs_melted, aes(x = PC, y = Value, fill = Type)) +
  geom_boxplot() +
  labs(title = "Boxplot of the First 10 Principal Components",
       x = "Principal Component",
       y = "Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#identify PCs with no overlap in IQRs
iqr_diff <- pca_data_10pcs %>%
  group_by(Type) %>%
  summarise(across(starts_with("PC"), ~ IQR(.))) %>%
  pivot_longer(cols = -Type, names_to = "PC", values_to = "IQR") %>%
  spread(Type, IQR) %>%
  mutate(no_overlap = (benign + malignant < abs(median(brca$x[, 1]))))

#which PCS dont overlap IQRs for benign and malignant
no_overlap_pcs <- iqr_diff %>% filter(no_overlap) %>% pull(PC)


#Section 3
#create data partition
set.seed(1)
test_index <- createDataPartition(brca$y, times = 1, p = 0.2, list = FALSE)
test_x <- x_scaled[test_index,]
test_y <- brca$y[test_index]
train_x <- x_scaled[-test_index,]
train_y <- brca$y[-test_index]

#get the proportion of benign sasmples in train and test sets
mean(train_y=="B")
mean(test_y=="B")

#log reg. model
set.seed(1)
log_reg_model <- train(train_x, train_y, method = "glm", family = binomial)
log_reg_predictions <- predict(log_reg_model, test_x)
log_reg_accuracy <- mean(log_reg_predictions == test_y)

#loess model
set.seed(5)
if (!requireNamespace("gam", quietly = TRUE)) {
  install.packages("gam")
}
library(gam)
loess_model <- train(train_x, train_y, method = "gamLoess")

#generate predictions on the test set
loess_predictions <- predict(loess_model, test_x)

#calculate and print the accuracy of the loess model
loess_accuracy <- mean(loess_predictions == test_y)

#Section 4
#train k-nearest neighbors model
set.seed(7)
knn_grid <- expand.grid(k = seq(3, 21, 2))
knn_model <- train(train_x, train_y, method = "knn", tuneGrid = knn_grid)

#generate predictions and calculate accuracy
knn_predictions <- predict(knn_model, test_x)
knn_accuracy <- mean(knn_predictions == test_y)
final_k <- knn_model$bestTune$k

#train random forest model
set.seed(9)
rf_grid <- expand.grid(mtry = c(3, 5, 7, 9))
rf_model <- train(train_x, train_y, method = "rf", tuneGrid = rf_grid, importance = TRUE)

#generate predictions and calculate accuracy
rf_predictions <- predict(rf_model, test_x)
rf_accuracy <- mean(rf_predictions == test_y)
best_mtry <- rf_model$bestTune$mtry
importance <- varImp(rf_model, scale = FALSE)
most_important_variable <- rownames(importance$importance)[which.max(importance$importance$Overall)]

#extract the top 10 most important variables
top_10_importance <- head(importance$importance[order(-importance$importance$Overall), ], 10)

#determine the feature set
feature_set <- function(variable) {
  if (grepl("mean", variable)) {
    return("mean values")
  } else if (grepl("se", variable)) {
    return("standard errors")
  } else if (grepl("worst", variable)) {
    return("worst values")
  } else {
    return("other")
  }
}

top_10_feature_set <- table(sapply(rownames(top_10_importance), feature_set))
most_important_feature_set <- names(top_10_feature_set)[which.max(top_10_feature_set)]

#generate predictions for all models
log_reg_predictions <- predict(log_reg_model, test_x)
loess_predictions <- predict(loess_model, test_x)
knn_predictions <- predict(knn_model, test_x)
rf_predictions <- predict(rf_model, test_x)

#combine predictions into a data frame
ensemble_predictions <- data.frame(log_reg_predictions, loess_predictions, knn_predictions, rf_predictions)

#create a majority vote prediction
majority_vote <- apply(ensemble_predictions, 1, function(row) {
  ifelse(mean(row == "M") > 0.5, "M", "B")
})

#calculate the accuracy of the ensemble
ensemble_accuracy <- mean(majority_vote == test_y)

#table of accuracies
accuracies <- data.frame(
  Model = c("Logistic Regression", "Loess", "kNN", "Random Forest", "Ensemble"),
  Accuracy = c(log_reg_accuracy, loess_accuracy, knn_accuracy, rf_accuracy, ensemble_accuracy)
)

#print the table
print(accuracies)

#model with the highest accuracy
best_model <- accuracies$Model[which.max(accuracies$Accuracy)]

#combine everything to get the ensemble accuracy and compare to others
set.seed(9)
rf_grid <- expand.grid(mtry = c(3, 5, 7, 9))
rf_model <- train(train_x, train_y, method = "rf", tuneGrid = rf_grid, importance = TRUE)
rf_predictions <- predict(rf_model, test_x)
rf_accuracy <- mean(rf_predictions == test_y)
importance <- varImp(rf_model, scale = FALSE)
most_important_variable <- rownames(importance$importance)[which.max(importance$importance$Overall)]

ensemble_predictions <- data.frame(
  log_reg_predictions,
  loess_predictions,
  knn_predictions,
  rf_predictions
)

majority_vote <- apply(ensemble_predictions, 1, function(row) {
  ifelse(mean(row == "M") > 0.5, "M", "B")
})

ensemble_accuracy <- mean(majority_vote == test_y)
