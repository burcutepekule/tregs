rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(lhs)
library(readr)

# Define the original parameter bounds (matching mass_simulation_LHS.py)
param_bounds = list(
  th_ROS_microbe = c(0, 1),
  th_ROS_epith_recover = c(0, 1),
  epith_recovery_chance = c(0, 1),
  rat_com_pat_threshold = c(0.5, 1),
  diffusion_speed_DAMPs = c(0, 0.12),
  diffusion_speed_SAMPs = c(0, 0.12),
  diffusion_speed_ROS = c(0, 0.12),
  add_ROS = c(0, 1),
  add_DAMPs = c(0, 1),
  add_SAMPs = c(0, 1),
  ros_decay = c(0, 1),
  DAMPs_decay = c(0, 1),
  SAMPs_decay = c(0, 1),
  activation_threshold_DAMPs = c(0, 1),
  activation_threshold_SAMPs = c(0, 1),
  activity_engulf_M0_baseline = c(0, 0.5),
  activity_engulf_M1_baseline = c(0, 0.5),
  activity_engulf_M2_baseline = c(0, 0.5),
  activity_ROS_M1_baseline = c(0, 0.5),
  rate_leak_commensal_injury = c(0.5, 1),
  rate_leak_pathogen_injury = c(0.5, 1),
  rate_leak_commensal_baseline = c(0, 0.25),
  active_age_limit = c(3, 30),  # discrete parameter, will be rounded
  treg_discrimination_efficiency = c(0, 1)
)

param_names = names(param_bounds)

# Set parameters
set.seed(123)
n_samples = 10000  # Total number of LHS samples to generate

# Generate LHS samples from original bounds
cat("Generating", n_samples, "LHS samples from original parameter bounds...\n")

n_params = length(param_names)
lhs_unit = randomLHS(n_samples, n_params)

# Create dataframe for samples
lhs_samples = as.data.frame(lhs_unit)
names(lhs_samples) = param_names

# Scale each parameter to its original bounds
for (param in param_names) {
  param_min = param_bounds[[param]][1]
  param_max = param_bounds[[param]][2]
  lhs_samples[[param]] = lhs_samples[[param]] * (param_max - param_min) + param_min
}

# Round the discrete parameter
lhs_samples$active_age_limit = round(lhs_samples$active_age_limit)

# Add parameter set ID
lhs_samples$param_set_id = 0:(nrow(lhs_samples) - 1)

# Reorder columns to put param_set_id first
lhs_samples = lhs_samples[c('param_set_id', param_names)]

# Export
output_file = "original_lhs_parameters.csv"
write.csv(lhs_samples, output_file, row.names = FALSE)
cat("\nDataset saved to:", output_file, "\n")
cat("Total parameter sets:", nrow(lhs_samples), "\n")

# # Optional: Compare with existing data if available
# # Load previous data for comparison (from A006 script)
# df_params_file = '/Users/burcutepekule/Desktop/tregs/mass_sim_results/sampled_parameters.csv'
# if (file.exists(df_params_file)) {
#   df_params = read_csv(df_params_file, show_col_types = FALSE)
# 
#   # Prepare the data for comparison
#   original_data = df_params %>%
#     select(all_of(param_names)) %>%
#     mutate(source = "Original Python LHS")
# 
#   lhs_data = lhs_samples %>%
#     select(all_of(param_names)) %>%
#     mutate(source = "R LHS Original Bounds")
# 
#   # Combine
#   comparison_data = bind_rows(original_data, lhs_data)
# 
#   # Pivot to long format for faceting
#   comparison_long = comparison_data %>%
#     pivot_longer(cols = all_of(param_names),
#                  names_to = "parameter",
#                  values_to = "value")
# 
#   # Create the overlaid density plots
#   density_plot = ggplot(comparison_long, aes(x = value, fill = source, color = source)) +
#     geom_density(alpha = 0.4, linewidth = 0.8) +
#     facet_wrap(~parameter, scales = "free", ncol = 4) +
#     scale_fill_manual(values = c("Original Python LHS" = "gray60",
#                                   "R LHS Original Bounds" = "steelblue"),
#                       name = "Dataset") +
#     scale_color_manual(values = c("Original Python LHS" = "gray30",
#                                    "R LHS Original Bounds" = "darkblue"),
#                        name = "Dataset") +
#     labs(title = "Parameter Distributions: Python LHS vs R LHS (Original Bounds)",
#          subtitle = "Overlaid density plots for all 24 parameters",
#          x = "Parameter Value",
#          y = "Density") +
#     theme_minimal() +
#     theme(
#       legend.position = "top",
#       strip.text = element_text(size = 8),
#       axis.text = element_text(size = 7)
#     )
# 
#   print(density_plot)
# 
#   # Save as high-res image
#   ggsave("parameter_distributions_original_bounds_comparison.png",
#          density_plot,
#          width = 16,
#          height = 12,
#          dpi = 300)
# 
#   cat("\nComparison plot saved to: parameter_distributions_original_bounds_comparison.png\n")
# } else {
#   cat("\nOriginal Python LHS file not found, skipping comparison plot.\n")
# }
# 
# # Summary statistics
# cat("\n=== SUMMARY STATISTICS ===\n")
# summary_stats = lhs_samples %>%
#   select(all_of(param_names)) %>%
#   summarise(across(everything(), list(
#     min = ~min(.),
#     max = ~max(.),
#     mean = ~mean(.),
#     median = ~median(.)
#   )))
# 
# cat("\nParameter ranges:\n")
# for (param in param_names) {
#   min_val = param_bounds[[param]][1]
#   max_val = param_bounds[[param]][2]
#   actual_min = min(lhs_samples[[param]])
#   actual_max = max(lhs_samples[[param]])
#   cat(sprintf("  %-35s Expected: [%.3f, %.3f]  Actual: [%.3f, %.3f]\n",
#               param, min_val, max_val, actual_min, actual_max))
# }
