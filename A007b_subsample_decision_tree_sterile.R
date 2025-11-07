rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)  # For read_csv
library(stringr)
library(zoo)
library(mgcv)

inj_type             = 'sterile'
ss_start_threshold   = 250
t_max                = 500
tol_in               = 25*0.25

# source('~/Desktop/tregs/A006_decision_tree_better.R')
# source('~/Desktop/tregs/A006_decision_tree_worse.R')
# source('~/Desktop/tregs/A006_decision_tree_drift.R')

source("/Users/burcutepekule/Dropbox/Treg_problem_v2/MISC/PLOT_FUNCTIONS.R")
df_raw    = readRDS('/Users/burcutepekule/Desktop/tregs/all_comparison_results_0.rds')
df_params = read_csv('/Users/burcutepekule/Desktop/tregs/mass_sim_results/sampled_parameters.csv', show_col_types = FALSE)

df_model = df_raw %>% dplyr::filter(comparison=='Treg_OFF_ON' & injury_type==inj_type)
df_model = inner_join(df_model, df_params, by='param_set_id')

#----- filter based on ss_start, it cannot be too large otherwise not much to compare!
param_id_all_below = df_model %>%
  dplyr::group_by(param_set_id) %>%
  dplyr::summarise(all_below = all(ss_start < ss_start_threshold), .groups = "drop") %>%
  dplyr::filter(all_below) %>%
  dplyr::pull(param_set_id)
df_model = df_model %>% dplyr::filter(param_set_id %in% param_id_all_below)

#----- filter based on replicate_id, less than 10 means incomplete!
param_id_all_complete = df_model %>%
  dplyr::group_by(param_set_id, comparison, injury_type) %>%
  dplyr::summarise(all_complete = (n_distinct(replicate_id) == 10), .groups = "drop") %>%
  dplyr::group_by(param_set_id) %>%
  dplyr::summarise(all_complete = all(all_complete), .groups = "drop") %>%
  dplyr::filter(all_complete) %>%
  dplyr::pull(param_set_id)
df_model = df_model %>% dplyr::filter(param_set_id %in% param_id_all_complete)

#------- Effect size based on cohens d -----------------------------------------
df_model = df_model %>% dplyr::mutate(abs_cohens_d = abs(cohens_d))
df_model = df_model %>% dplyr::mutate(effect_size = case_when(
  abs_cohens_d < 0.2 ~ "Negligible",
  abs_cohens_d < 0.5 & abs_cohens_d>= 0.2  ~ "Small",
  abs_cohens_d < 0.8 & abs_cohens_d>= 0.5 ~ "Medium",
  TRUE ~ "Large"
))

df_model$tol = tol_in
df_summary = df_model %>%
  # dplyr::group_by(param_set_id, injury_type, comparison) %>%
  dplyr::group_by(param_set_id, comparison) %>%
  dplyr::summarise(
    n_better = sum(effect_size %in% c('Large','Medium') & mean_diff > tol, na.rm = TRUE),
    n_drift  = sum((mean_diff <= tol & mean_diff >= -1*tol), na.rm = TRUE),
    n_worse  = sum(effect_size %in% c('Large','Medium') & mean_diff < -1*tol, na.rm = TRUE),
    .groups = "drop"
  )

df_summary = inner_join(df_summary %>% dplyr::select(param_set_id, n_better, n_drift, n_worse), df_model, by='param_set_id')

param_names = c(
  "th_ROS_microbe",
  "th_ROS_epith_recover",
  "epith_recovery_chance",
  "rat_com_pat_threshold",
  "diffusion_speed_DAMPs",
  "diffusion_speed_SAMPs",
  "diffusion_speed_ROS",
  "add_ROS",
  "add_DAMPs",
  "add_SAMPs",
  "ros_decay",
  "DAMPs_decay",
  "SAMPs_decay",
  "activation_threshold_DAMPs",
  "activation_threshold_SAMPs",
  "activity_engulf_M0_baseline",
  "activity_engulf_M1_baseline",
  "activity_engulf_M2_baseline",
  "activity_ROS_M1_baseline",
  "rate_leak_commensal_injury",
  "rate_leak_pathogen_injury",
  "rate_leak_commensal_baseline",
  "active_age_limit",
  "treg_discrimination_efficiency"
)

df_summary = df_summary %>%
  dplyr::select(n_better, n_drift, n_worse, all_of(param_names)) %>%
  distinct()

library(ggplot2)
library(ggdendro)
library(factoextra)
library(cluster)
library(partykit)

# Create outcome proportions
df_clustering = df_summary %>%
  mutate(
    total = n_better + n_drift + n_worse,
    frac_better = div0(n_better,total),
    frac_drift = div0(n_drift,total),
    frac_worse = div0(n_worse,total)
    )

unique(df_clustering$n_better+df_clustering$n_drift+df_clustering$n_worse) # all 10

ctree_model_better = readRDS('ctree_model_better.rds')
ctree_model_worse  = readRDS('ctree_model_worse.rds')
ctree_model_drift  = readRDS('ctree_model_drift.rds')

# ------------------------------------------------------------------------------
library(dplyr)
library(lhs)

# Step 0: Define total samples first!
n_total      = 10000
# n_per_category = n_total / 3  # 3000 each
n_per_better = 5000
n_per_worse  = 5000
n_per_drift  = n_total-(n_per_better+n_per_worse)

# Step 1: Get node assignments and calculate outcome rates for each node
df_with_nodes = df_clustering %>%
  mutate(
    node_better = predict(ctree_model_better, type = "node"),
    node_worse = predict(ctree_model_worse, type = "node"),
    node_drift = predict(ctree_model_drift, type = "node")
  )

# Calculate the RATE of better outcomes in each node of the better tree
node_stats_better = df_with_nodes %>%
  group_by(node_better) %>%
  summarise(
    n_orig = n(),
    mean_better = mean(frac_better),     
    mean_worse = mean(frac_worse),
    mean_drift = mean(frac_drift),
    tree_type  = 'better',
    .groups = "drop"
  )
colnames(node_stats_better)[1] = 'node_index'

node_stats_worse = df_with_nodes %>%
  group_by(node_worse) %>%
  summarise(
    n_orig = n(),
    mean_better = mean(frac_better),     
    mean_worse = mean(frac_worse),
    mean_drift = mean(frac_drift),
    tree_type  = 'worse',
    .groups = "drop"
  )
colnames(node_stats_worse)[1] = 'node_index'

node_stats_drift = df_with_nodes %>%
  group_by(node_drift) %>%
  summarise(
    n_orig = n(),
    mean_better = mean(frac_better),
    mean_worse = mean(frac_worse),
    mean_drift = mean(frac_drift),      
    tree_type  = 'drift',
    .groups = "drop"
  )
colnames(node_stats_drift)[1] = 'node_index'

node_stats = rbind(node_stats_drift, node_stats_worse, node_stats_better)

#-------------------------------------------------------------------------------
library(lhs)  # For Latin Hypercube Sampling
library(dplyr)
library(tidyr)

# Step 1: Function to extract parameter bounds for a specific node
get_node_bounds = function(data_with_nodes, node_col, node_id, param_names) {
  # Filter data to this specific node
  node_data = data_with_nodes %>%
    filter(!!sym(node_col) == node_id)
  
  if (nrow(node_data) == 0) {
    warning(paste("No data found for node", node_id))
    return(NULL)
  }
  
  # Get min/max for each parameter
  # Use 5th and 95th percentiles to avoid extreme outliers
  bounds = node_data %>%
    select(all_of(param_names)) %>%
    summarise(across(everything(), list(
      min = ~quantile(., 0.05, na.rm = TRUE),
      max = ~quantile(., 0.95, na.rm = TRUE)
    ))) %>%
    pivot_longer(everything(), 
                 names_to = "param_stat", 
                 values_to = "value") %>%
    separate(param_stat, into = c("parameter", "stat"), sep = "_(?=[^_]+$)") %>%
    pivot_wider(names_from = stat, values_from = value)
  
  return(bounds)
}

# Step 2: Function to generate LHS samples within parameter bounds
generate_lhs_samples = function(bounds, n_samples, param_names) {
  if (is.null(bounds) || n_samples == 0) {
    return(NULL)
  }
  
  n_params = length(param_names)
  
  # Generate LHS in [0,1]^n_params space
  # Each column is one parameter, each row is one parameter set
  lhs_unit = randomLHS(n_samples, n_params)
  
  # Transform from [0,1] to actual parameter ranges
  param_samples = as.data.frame(lhs_unit)
  names(param_samples) = param_names
  
  # Scale each parameter to its actual range
  for (param in param_names) {
    param_min = bounds$min[bounds$parameter == param]
    param_max = bounds$max[bounds$parameter == param]
    
    # Transform: value_actual = value_[0,1] * (max - min) + min
    param_samples[[param]] = param_samples[[param]] * (param_max - param_min) + param_min
  }
  
  return(param_samples)
}

# Step 3: Generate samples for BETTER category
cat("Generating samples for BETTER category...\n")
better_samples_list = list()

for (i in 1:nrow(node_stats_better)) {
  node_id = node_stats_better$node_better[i]
  n_samples = node_stats_better$n_to_sample[i]
  
  if (n_samples > 0) {
    # Get parameter bounds for this node
    bounds = get_node_bounds(df_with_nodes, "node_better", node_id, param_names)
    
    # Generate LHS samples within these bounds
    samples = generate_lhs_samples(bounds, n_samples, param_names)
    
    if (!is.null(samples)) {
      # Add metadata
      samples$target_category = "better"
      samples$source_node = node_id
      samples$expected_better_rate = node_stats_better$mean_better[i]
      
      better_samples_list[[i]] = samples
      
      cat(sprintf("  Node %d: generated %d samples (expected better rate: %.1f%%)\n", 
                  node_id, n_samples, node_stats_better$mean_better[i] * 100))
    }
  }
}

better_samples = bind_rows(better_samples_list)
cat(sprintf("Total BETTER samples: %d\n\n", nrow(better_samples)))

# Step 4: Generate samples for WORSE category
cat("Generating samples for WORSE category...\n")
worse_samples_list = list()

for (i in 1:nrow(node_stats_worse)) {
  node_id = node_stats_worse$node_worse[i]
  n_samples = node_stats_worse$n_to_sample[i]
  
  if (n_samples > 0) {
    bounds = get_node_bounds(df_with_nodes, "node_worse", node_id, param_names)
    samples = generate_lhs_samples(bounds, n_samples, param_names)
    
    if (!is.null(samples)) {
      samples$target_category = "worse"
      samples$source_node = node_id
      samples$expected_worse_rate = node_stats_worse$mean_worse[i]
      
      worse_samples_list[[i]] = samples
      
      cat(sprintf("  Node %d: generated %d samples (expected worse rate: %.1f%%)\n", 
                  node_id, n_samples, node_stats_worse$mean_worse[i] * 100))
    }
  }
}

worse_samples = bind_rows(worse_samples_list)
cat(sprintf("Total WORSE samples: %d\n\n", nrow(worse_samples)))

# Step 5: Generate samples for DRIFT category
cat("Generating samples for DRIFT category...\n")
drift_samples_list = list()

for (i in 1:nrow(node_stats_drift)) {
  node_id = node_stats_drift$node_drift[i]
  n_samples = node_stats_drift$n_to_sample[i]
  
  if (n_samples > 0) {
    bounds = get_node_bounds(df_with_nodes, "node_drift", node_id, param_names)
    samples = generate_lhs_samples(bounds, n_samples, param_names)
    
    if (!is.null(samples)) {
      samples$target_category = "drift"
      samples$source_node = node_id
      samples$expected_drift_rate = node_stats_drift$mean_drift[i]
      
      drift_samples_list[[i]] = samples
      
      cat(sprintf("  Node %d: generated %d samples (expected drift rate: %.1f%%)\n", 
                  node_id, n_samples, node_stats_drift$mean_drift[i] * 100))
    }
  }
}

drift_samples = bind_rows(drift_samples_list)
cat(sprintf("Total DRIFT samples: %d\n\n", nrow(drift_samples)))

# Step 6: Combine all samples
balanced_lhs_dataset = bind_rows(better_samples, worse_samples, drift_samples)

# Shuffle the rows 
balanced_lhs_dataset = balanced_lhs_dataset %>% slice_sample(n = nrow(.), replace = FALSE)

cat("=== FINAL SUMMARY ===\n")
cat(sprintf("Total parameter sets generated: %d\n", nrow(balanced_lhs_dataset)))
print(table(balanced_lhs_dataset$target_category))

# Step 7: Verify parameter ranges are reasonable
cat("\nParameter range check:\n")
print(summary(balanced_lhs_dataset %>% select(all_of(param_names))))

# Step 8: Export for your simulation
balanced_lhs_dataset$param_set_id = 0:(dim(balanced_lhs_dataset)[1]-1)
balanced_lhs_dataset = balanced_lhs_dataset[c('param_set_id',param_names)]
write.csv(balanced_lhs_dataset, "balanced_lhs_parameters_for_simulation.csv")

cat("\nDataset exported to: balanced_lhs_parameters_for_simulation.csv\n")

#----------------------------- Plots -------------------------------------------
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)  # For arranging multiple plots

# Prepare the data for comparison
original_data = df_clustering %>%
  select(all_of(param_names)) %>%
  mutate(source = "Original")

lhs_data = balanced_lhs_dataset %>%
  select(all_of(param_names)) %>%
  mutate(source = "LHS Generated")

# Combine
comparison_data = bind_rows(original_data, lhs_data)

# Pivot to long format for faceting
comparison_long = comparison_data %>%
  pivot_longer(cols = all_of(param_names), 
               names_to = "parameter", 
               values_to = "value")

# Create the overlaid density plots
density_plot = ggplot(comparison_long, aes(x = value, fill = source, color = source)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  facet_wrap(~parameter, scales = "free", ncol = 4) +
  scale_fill_manual(values = c("Original" = "gray60", "LHS Generated" = "steelblue"),
                    name = "Dataset") +
  scale_color_manual(values = c("Original" = "gray30", "LHS Generated" = "darkblue"),
                     name = "Dataset") +
  labs(title = "Parameter Distributions: Original vs LHS Generated",
       subtitle = "Overlaid density plots for all 24 parameters",
       x = "Parameter Value",
       y = "Density") +
  theme_minimal() +
  theme(
    legend.position = "top",
    strip.text = element_text(size = 8),
    axis.text = element_text(size = 7)
  )

print(density_plot)

# Save as high-res image
ggsave("parameter_distributions_comparison.png", 
       density_plot, 
       width = 16, 
       height = 12, 
       dpi = 300)
