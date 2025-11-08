rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)
library(partykit)

source("/Users/burcutepekule/Dropbox/Treg_problem_v2/MISC/PLOT_FUNCTIONS.R")
df_raw    = readRDS('/Users/burcutepekule/Desktop/tregs/all_comparison_results_0.rds')
df_params = read_csv('/Users/burcutepekule/Desktop/tregs/mass_sim_results/sampled_parameters.csv', show_col_types = FALSE)

inj_type             = 'sterile'
ss_start_threshold   = 450
t_max                = 500
tol_in               = 25*0.25

# [Lines 22-123: Same filtering and preprocessing as original...]
df_model = df_raw %>% dplyr::filter(comparison=='Treg_OFF_ON' & injury_type==inj_type)
df_model = inner_join(df_model, df_params, by='param_set_id')

param_id_all_below = df_model %>%
  dplyr::group_by(param_set_id) %>%
  dplyr::summarise(all_below = all(ss_start < ss_start_threshold), .groups = "drop") %>%
  dplyr::filter(all_below) %>%
  dplyr::pull(param_set_id)
df_model = df_model %>% dplyr::filter(param_set_id %in% param_id_all_below)

df_model = df_model %>% dplyr::mutate(abs_cohens_d = abs(cohens_d))
df_model = df_model %>% dplyr::mutate(effect_size = case_when(
  abs_cohens_d < 0.2 ~ "Negligible",
  abs_cohens_d < 0.5 & abs_cohens_d>= 0.2  ~ "Small",
  abs_cohens_d < 0.8 & abs_cohens_d>= 0.5 ~ "Medium",
  TRUE ~ "Large"
))

df_model$tol = tol_in
df_summary = df_model %>%
  dplyr::group_by(param_set_id, injury_type, comparison) %>%
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

# Create outcome proportions and dominant outcome
df_clustering = df_summary %>%
  mutate(
    total = n_better + n_drift + n_worse,
    frac_better = div0(n_better,total),
    frac_drift = div0(n_drift,total),
    frac_worse = div0(n_worse,total),
    dominant_outcome = case_when(
      frac_better > frac_worse & frac_better > frac_drift ~ "better",
      frac_worse > frac_better & frac_worse > frac_drift ~ "worse",
      TRUE ~ "drift"
    ),
    dominant_outcome = factor(dominant_outcome, levels = c("better", "drift", "worse"))
  )

cat("\n=== ORIGINAL OUTCOME DISTRIBUTION ===\n")
print(table(df_clustering$dominant_outcome))
cat("\nProportions:\n")
print(100 * prop.table(table(df_clustering$dominant_outcome)))

# === MODIFIED APPROACH: Direct resampling from original parameter sets ===
# Instead of LHS within bounds, directly resample from existing param sets by outcome

set.seed(123)

# Separate parameter sets by dominant outcome
params_better = df_clustering %>% filter(dominant_outcome == "better")
params_worse = df_clustering %>% filter(dominant_outcome == "worse")
params_drift = df_clustering %>% filter(dominant_outcome == "drift")

cat("\n=== AVAILABLE PARAMETER SETS ===\n")
cat("Better:", nrow(params_better), "\n")
cat("Worse:", nrow(params_worse), "\n")
cat("Drift:", nrow(params_drift), "\n")

# Define target counts (aiming for 1/3 each)
target_total = 10000  # Adjust as needed
target_per_category = target_total / 3

# Sample with replacement if needed to reach target
sample_with_target = function(df, target_n, param_cols) {
  if (nrow(df) == 0) {
    cat("WARNING: No samples available for this category!\n")
    return(NULL)
  }

  # Sample with replacement if target > available
  replace_flag = (target_n > nrow(df))

  sampled_indices = sample(1:nrow(df), size = target_n, replace = replace_flag)
  sampled = df[sampled_indices, param_cols]

  if (replace_flag) {
    cat("  (Sampled with replacement due to insufficient original samples)\n")
  }

  return(sampled)
}

cat("\n=== GENERATING BALANCED DATASET ===\n")

# Sample each category
cat("Sampling better...\n")
balanced_better = sample_with_target(params_better, round(target_per_category), param_names)
if (!is.null(balanced_better)) balanced_better$target_category = "better"

cat("Sampling worse...\n")
balanced_worse = sample_with_target(params_worse, round(target_per_category), param_names)
if (!is.null(balanced_worse)) balanced_worse$target_category = "worse"

cat("Sampling drift...\n")
balanced_drift = sample_with_target(params_drift, round(target_per_category), param_names)
if (!is.null(balanced_drift)) balanced_drift$target_category = "drift"

# Combine
balanced_lhs_dataset = bind_rows(balanced_better, balanced_worse, balanced_drift)

# Shuffle
balanced_lhs_dataset = balanced_lhs_dataset %>%
  slice_sample(n = nrow(.))

# Add param_set_id but KEEP target_category for validation
balanced_lhs_dataset$param_set_id = 0:(nrow(balanced_lhs_dataset)-1)

cat("\n=== FINAL BALANCED DATASET ===\n")
cat("Total parameter sets:", nrow(balanced_lhs_dataset), "\n")
print(table(balanced_lhs_dataset$target_category))
cat("\nProportions:\n")
print(100 * prop.table(table(balanced_lhs_dataset$target_category)))

# Save WITH target_category for validation
write.csv(balanced_lhs_dataset, "balanced_lhs_parameters_DIRECT.csv", row.names = FALSE)
cat("\nDataset saved to: balanced_lhs_parameters_DIRECT.csv\n")
cat("(Note: target_category column is INCLUDED for validation)\n")
