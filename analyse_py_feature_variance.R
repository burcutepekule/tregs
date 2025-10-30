rm(list=ls())
library(dplyr)
library(tidyr)
source("/Users/burcutepekule/Dropbox/Treg_problem_v2/MISC/PLOT_FUNCTIONS.R")
scenarios      = read_csv('/Users/burcutepekule/Desktop/tregs/mass_sim_results/scenarios.csv', show_col_types = FALSE)

params           = read_csv('/Users/burcutepekule/Desktop/tregs/mass_sim_results/sampled_parameters.csv', show_col_types = FALSE)
df_short_merged_with_params = readRDS('/Users/burcutepekule/Desktop/tregs/df_short_merged_with_params.rds')
df_short_merged_with_params_keep = df_short_merged_with_params
param_cols = colnames(df_short_merged_with_params_keep)[8:31]

# Prepare data
lda_data = df_short_merged_with_params %>%
  dplyr::select(class_description, all_of(param_cols)) %>%
  filter(!is.na(class_description))  # Remove any NA classes

features = lda_data %>% dplyr::select(-class_description)
features = as.data.frame(scale(features)) #Each feature has mean = 0 and standard deviation = 1, 
# Features with different units/scales become comparable
# Features with larger natural variance don't dominate the analysis

# bind back
df_short_merged_with_params = cbind(class_description = lda_data$class_description, features)
table(df_short_merged_with_params$class_description)

# Get class frequencies first
class_frequencies = df_short_merged_with_params %>%
  filter(!is.na(class_description)) %>%
  count(class_description) %>%
  mutate(freq = n / sum(n))

# Calculate your original summary
feature_summary = df_short_merged_with_params %>%
  dplyr::select(class_description, all_of(param_cols)) %>%
  filter(!is.na(class_description)) %>%
  group_by(class_description) %>%
  summarise(across(all_of(param_cols), 
                   list(mean = mean),
                   .names = "{.col}_{.fn}"))
  # summarise(across(all_of(param_cols), 
  #                  list(min = min, 
  #                       mean = mean, 
  #                       max = max),
  #                  .names = "{.col}_{.fn}"))

feature_summary_orig = df_short_merged_with_params_keep %>%
  dplyr::select(class_description, all_of(param_cols)) %>%
  filter(!is.na(class_description)) %>%
  group_by(class_description) %>%
  summarise(across(all_of(param_cols), 
                   list(mean = mean),
                   .names = "{.col}_{.fn}"))

feature_summary_long = feature_summary %>%
  pivot_longer(cols = -class_description, 
               names_to = "feature_stat", 
               values_to = "value_scaled")

feature_summary_orig_long = feature_summary_orig %>%
  pivot_longer(cols = -class_description, 
               names_to = "feature_stat", 
               values_to = "value")

# Calculate variance of each column (across the class means, mins, maxes)
variance_row = feature_summary %>%
  dplyr::select(-class_description) %>%
  summarise(across(everything(), var, na.rm = TRUE))

variance_row_long = variance_row %>%
  pivot_longer(cols = everything(),  
               names_to = "feature_stat", 
               values_to = "variance")

variance_row_long$variance=as.numeric(variance_row_long$variance)
feature_summary = merge(feature_summary_long, variance_row_long, by='feature_stat')
feature_summary = merge(feature_summary, feature_summary_orig_long, by=c('feature_stat','class_description'))

