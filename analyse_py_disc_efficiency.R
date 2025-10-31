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
library(caret)

df_short_merged_with_params = readRDS('/Users/burcutepekule/Desktop/tregs/df_short_merged_with_params.rds')


table(df_short_merged_with_params$class_description)
df_short_merged_with_params = df_short_merged_with_params %>% dplyr::filter(class_description!="Sterile worse only")
df_short_merged_with_params = df_short_merged_with_params %>% dplyr::filter(class_description!="Sterile better, pathogenic worse")
length(unique(df_short_merged_with_params$param_set_id))
dim(df_short_merged_with_params)

table(df_short_merged_with_params$class_description)

# Upsampling
df_short_merged_with_params = upSample(x = df_short_merged_with_params[, -which(names(df_short_merged_with_params) == "class_description")],
                        y = as.factor(df_short_merged_with_params$class_description)) # class_description changed to Class
colnames(df_short_merged_with_params)[31] = 'class_description'
table(df_short_merged_with_params$class_description)

param_cols = colnames(df_short_merged_with_params)[7:30]
# Prepare data
lda_data = df_short_merged_with_params %>%
  dplyr::select(class_description, all_of(param_cols)) %>%
  filter(!is.na(class_description))  # Remove any NA classes

features = lda_data %>% dplyr::select(-class_description)
# features = as.data.frame(scale(features))

library(rpart)
library(rpart.plot)
library(partykit)

# Prepare data
tree_data = cbind(class_description = lda_data$class_description, features)

# Create a new factor with shorter, clear labels
tree_data$class_description = factor(tree_data$class_description,
                                           levels = c("Both better", 
                                                      "Pathogenic better only",
                                                      "Sterile better only",
                                                      "No effect",
                                                      "Pathogenic worse only"
                                           ),
                                           labels = c("P+S+",     
                                                      "P+",       
                                                      "S+",      
                                                      "P0S0",
                                                      "P-"))     

# Fit decision tree
tree_model = rpart(class_description ~ ., 
                   data = tree_data,
                   method = "class",
                   control = rpart.control(cp = 0.03, # complexity parameter (smaller, more complex)
                                           minsplit = 30, # min observations for split, 20 default
                                           maxdepth = 4)) # max tree depth

# First, check the order of your classes
classes = levels(as.factor(tree_data$class_description))
print(classes)  # See the order

palette_list = list("Greens", "Grays", "Purples","Reds", "Blues", "Oranges")

# Visualize tree with rules
pdf("tree_model.pdf", width = 24, height = 8)
rpart.plot(tree_model, 
           type = 4, 
           extra = 104, # show probability and fraction
           under = TRUE,
           cex = 0.8,
           box.palette = palette_list)
dev.off()

# Get text rules
print(tree_model)

# Convert to party object and plot
party_tree = as.party(tree_model)

pdf("decision_tree.pdf", width = 22, height = 10)
plot(party_tree, tp_args = list(rot = 45))  # Now rotation should work fine with short labels
dev.off()

# Extract rules as text
library(rattle)
asRules(tree_model) 

# Variable importance
tree_model$variable.importance


# # Now cross-validation should work
# library(MLmetrics)
# library(caret)
# 
# # Set up cross-validation
# train_control <- trainControl(
#   method = "cv",           
#   number = 10,             
#   classProbs = TRUE,       
#   summaryFunction = multiClassSummary
# )
# 
# # Define hyperparameter grid
# hyper_grid <- expand.grid(
#   cp = seq(0.001, 0.05, by = 0.005)
# )
# 
# # Train model with CV
# set.seed(123)
# cv_model <- train(
#   class_description ~ .,
#   data = tree_data_clean,
#   method = "rpart",
#   trControl = train_control,
#   tuneGrid = hyper_grid,
#   metric = "Accuracy"
# )
# 
# # View best parameters
# print(cv_model$bestTune)
# plot(cv_model)