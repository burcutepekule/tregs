# Clear workspace
rm(list=ls())

# Load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(randomForest)
library(MASS)  # For LDA comparison if needed

# Load your data
df_short_merged_with_params = readRDS('/Users/burcutepekule/Desktop/tregs/df_short_merged_with_params.rds')
df_short_merged_with_params = df_short_merged_with_params %>% mutate(class_code = case_when(
  class_code == "0000" ~ "Drift", #"No effect",
  class_code == "1000" ~ "Favorable", #Sterile better only",
  class_code == "0100" ~ "Favorable", #"Pathogenic better only",
  class_code == "0010" ~ "Unfavorable", #"Sterile worse only",
  class_code == "0001" ~ "Unfavorable", #"Pathogenic worse only",
  class_code == "1100" ~ "Favorable", #"Both better",
  class_code == "0011" ~ "Unfavorable", #"Both worse",
  class_code == "1001" ~ "Unfavorable", # "Sterile better, pathogenic worse",
  class_code == "0110" ~ "Unfavorable", # "Pathogenic better, sterile worse",
  TRUE ~ "Other"
))

# Get parameter columns (columns 8-31)
param_cols = colnames(df_short_merged_with_params)[8:31]

# Check the data
dim(df_short_merged_with_params)
table(df_short_merged_with_params$class_code)

# Prepare data for modeling
model_data = df_short_merged_with_params %>%
  dplyr::select(class_code, all_of(param_cols)) %>%
  filter(!is.na(class_code))

# Scale features
features = model_data %>% dplyr::select(-class_code)
features_scaled = as.data.frame(scale(features))

# Split into train/test (80/20)
set.seed(123)
n_samples = nrow(model_data)
train_idx = sample(1:n_samples, size = 0.8 * n_samples)

train_data = features_scaled[train_idx, ]
train_labels = model_data$class_code[train_idx]
test_data = features_scaled[-train_idx, ]
test_labels = model_data$class_code[-train_idx]

# Check the splits
cat("Training samples:", nrow(train_data), "\n")
cat("Test samples:", nrow(test_data), "\n")
cat("Classes in training:", unique(train_labels), "\n")
cat("Classes in test:", unique(test_labels), "\n")

# Train Random Forest with importance=TRUE
cat("Training Random Forest model...\n")

table(train_labels)
library(caret)
rf_model = randomForest(
  x = train_data,
  y = as.factor(train_labels),
  ntree = 500,
  importance = TRUE,
  localImp = TRUE,
  sampsize = c("Drift" = 35, "Unfavorable" = 35, "Favorable" = 35)  # Sample equally from each class
)

# Check the model
print(rf_model)

# Quick check - did importance get calculated?
dim(rf_model$importance)
colnames(rf_model$importance)

# Function to calculate SHAP-like values manually
calculate_shap_manual <- function(rf_model, X_test, X_train, n_samples = 50, target_class = "0000") {
  
  # Use subset for speed
  n_test = min(nrow(X_test), n_samples)
  X_test_subset = X_test[1:n_test, ]
  
  n_features = ncol(X_test)
  shap_values = matrix(0, nrow = n_test, ncol = n_features)
  colnames(shap_values) = colnames(X_test)
  
  # Baseline prediction (mean prediction on training data)
  baseline = mean(predict(rf_model, X_train, type = "prob")[, target_class])
  
  cat("Calculating SHAP values for", n_test, "samples...\n")
  
  for(i in 1:n_test) {
    if(i %% 10 == 0) cat("  Processing sample", i, "/", n_test, "\n")
    
    instance = X_test_subset[i, , drop = FALSE]
    
    for(j in 1:n_features) {
      # Create comparison datasets
      X_with = X_train[sample(1:nrow(X_train), 100), ]  # Sample for speed
      X_with[, j] = rep(instance[, j], nrow(X_with))
      
      X_without = X_with
      X_without[, j] = rep(mean(X_train[, j]), nrow(X_without))
      
      # Calculate marginal contribution
      pred_with = mean(predict(rf_model, X_with, type = "prob")[, target_class])
      pred_without = mean(predict(rf_model, X_without, type = "prob")[, target_class])
      
      shap_values[i, j] = pred_with - pred_without
    }
  }
  
  return(list(
    shap_values = shap_values,
    baseline = baseline
  ))
}

# Get all class names
unique_classes = as.character(unique(train_labels))
print(unique_classes)  # Should be "Drift", "Favorable", "Unfavorable"

# Store SHAP results for all classes
all_shap_results = list()
all_shap_importance = list()

# Calculate SHAP for each class
for(target_class in unique_classes) {
  cat("\n========================================\n")
  cat("Calculating SHAP for class:", target_class, "\n")
  
  # Calculate SHAP values
  shap_result = calculate_shap_manual(
    rf_model = rf_model,
    X_test = test_data,
    X_train = train_data,
    n_samples = 20,  # You can increase this for more accuracy
    target_class = target_class
  )
  
  # Store results
  all_shap_results[[target_class]] = shap_result
  
  # Calculate importance for this class
  shap_importance = data.frame(
    feature = colnames(shap_result$shap_values),
    mean_abs_shap = colMeans(abs(shap_result$shap_values)),
    mean_shap = colMeans(shap_result$shap_values),
    class = target_class
  ) %>%
    arrange(desc(mean_abs_shap))
  
  all_shap_importance[[target_class]] = shap_importance
  
  cat("Baseline for", target_class, ":", shap_result$baseline, "\n")
  cat("Top 3 features:", head(shap_importance$feature, 3), "\n")
}

# Combine all importance scores
combined_importance = do.call(rbind, all_shap_importance)

cat("\n========================================\n")
cat("SHAP calculation complete for all classes!\n")

top_features_overall = combined_importance %>%
  group_by(feature) %>%
  summarise(avg_importance = mean(mean_abs_shap)) %>%
  arrange(desc(avg_importance)) %>%
  head(10) %>%
  pull(feature)

# 1. Create a plot showing signed (directional) SHAP values
importance_signed = combined_importance %>%
  filter(feature %in% top_features_overall)

p_signed = ggplot(importance_signed, 
                  aes(x = reorder(feature, abs(mean_shap)), y = mean_shap, fill = class)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_fill_manual(values = c("Drift" = "#AFAAB9", 
                               "Favorable" = "lightgreen", 
                               "Unfavorable" = "coral")) +
  labs(title = "Signed SHAP Values: How Features Push Toward Each Class",
       subtitle = "Positive = increases probability, Negative = decreases probability",
       x = "Feature", y = "Mean SHAP value",
       fill = "Class") +
  theme_minimal()

print(p_signed)

# 2. Create a diverging heatmap for signed values
p_heatmap_signed = ggplot(importance_signed, 
                          aes(x = class, y = feature, fill = mean_shap)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0,
                       name = "Mean SHAP\n(signed)") +
  geom_text(aes(label = round(mean_shap, 3)), size = 3) +
  labs(title = "Directional SHAP Values by Class",
       subtitle = "Red = feature increases class probability, Blue = decreases",
       x = "Class", y = "Feature") +
  theme_minimal()

print(p_heatmap_signed)



