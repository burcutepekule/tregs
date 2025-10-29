rm(list=ls())
library(dplyr)
library(tidyr)
library(cowplot)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(grid)
library(zoo)
library(slider)
library(purrr)
library(tibble)
library(RcppRoll)
library(MASS)
library(ggrepel)
library(randomForest)
library(vip)
library(pdp)

# Using your existing data preparation
df_short_merged_with_params = readRDS('/Users/burcutepekule/Desktop/tregs/df_short_merged_with_params.rds')

param_cols = colnames(df_short_merged_with_params)[8:31]

# Prepare data (same as your LDA prep)
shap_data = df_short_merged_with_params %>%
  dplyr::select(class_code, all_of(param_cols)) %>%
  filter(!is.na(class_code))

# Scale features
features = shap_data %>% dplyr::select(-class_code)
features_scaled = as.data.frame(scale(features))

# Split data for SHAP (need train/test split)
set.seed(123)
train_idx = sample(1:nrow(shap_data), 0.8 * nrow(shap_data))
train_data = features_scaled[train_idx, ]
train_labels = shap_data$class_code[train_idx]
test_data = features_scaled[-train_idx, ]
test_labels = shap_data$class_code[-train_idx]

# ============================================
# OPTION 1: Random Forest + TreeSHAP (Faster)
# ============================================

# Train Random Forest
rf_model = randomForest(
  x = train_data,
  y = as.factor(train_labels),
  ntree = 500,
  importance = TRUE,
  localImp = TRUE  # This enables local importance
)

# Extract importance directly from the model
imp_matrix = rf_model$importance

# Create importance dataframe using the direct access
importance_df = data.frame(
  feature = rownames(imp_matrix),
  MeanDecreaseGini = imp_matrix[, "MeanDecreaseGini"],
  MeanDecreaseAccuracy = imp_matrix[, "MeanDecreaseAccuracy"]
) %>%
  arrange(desc(MeanDecreaseGini))

# Plot RF feature importance
p_rf_importance = ggplot(importance_df[1:15,], aes(x = reorder(feature, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Random Forest Feature Importance",
       x = "Feature", y = "Mean Decrease in Gini") +
  theme_minimal()

print(p_rf_importance)

# ============================================
# LOCAL FEATURE IMPORTANCE & INTERPRETABILITY
# Using packages you already have
# ============================================

# 1. LIME-style local explanations (simplified version)
# This shows which features drive individual predictions

explain_single_prediction <- function(rf_model, instance, train_data, n_features = 10) {
  # Get original prediction
  orig_pred <- predict(rf_model, instance, type = "prob")
  orig_class <- predict(rf_model, instance)
  
  # For each feature, measure impact of changing it
  feature_impacts <- data.frame(
    feature = colnames(train_data),
    impact = numeric(ncol(train_data))
  )
  
  for(i in 1:ncol(train_data)) {
    # Replace feature with mean value
    modified_instance <- instance
    modified_instance[, i] <- mean(train_data[, i])
    
    # Get new prediction
    new_pred <- predict(rf_model, modified_instance, type = "prob")
    
    # Measure change in probability for predicted class
    class_col <- which(colnames(orig_pred) == as.character(orig_class))
    feature_impacts$impact[i] <- orig_pred[, class_col] - new_pred[, class_col]
  }
  
  feature_impacts <- feature_impacts %>%
    arrange(desc(abs(impact))) %>%
    head(n_features)
  
  return(list(
    predicted_class = as.character(orig_class),
    probability = orig_pred[, class_col],
    feature_impacts = feature_impacts
  ))
}

# 2. Analyze predictions for each class
# Let's look at examples from each class
unique_classes <- unique(test_labels)
example_explanations <- list()

for(class_label in unique_classes[1:min(5, length(unique_classes))]) {
  # Find examples of this class
  class_indices <- which(test_labels == class_label)
  
  if(length(class_indices) > 0) {
    # Take first example
    idx <- class_indices[1]
    
    explanation <- explain_single_prediction(
      rf_model, 
      test_data[idx, , drop = FALSE], 
      train_data
    )
    
    example_explanations[[class_label]] <- data.frame(
      class = class_label,
      sample_idx = idx,
      explanation$feature_impacts
    )
  }
}

# Combine all explanations
all_explanations <- do.call(rbind, example_explanations)

# Visualize feature impacts for different classes
p_local_importance <- ggplot(all_explanations, 
                             aes(x = impact, y = reorder(feature, abs(impact)))) +
  geom_bar(stat = "identity", aes(fill = impact > 0)) +
  facet_wrap(~ class, scales = "free_x") +
  scale_fill_manual(values = c("FALSE" = "red", "TRUE" = "blue"),
                    labels = c("Decreases prob", "Increases prob"),
                    name = "Effect") +
  labs(title = "Feature Impacts on Individual Predictions by Class",
       x = "Impact on Prediction Probability", y = "Feature") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 7))

print(p_local_importance)


