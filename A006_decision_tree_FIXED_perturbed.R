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

# [Same preprocessing as before...]
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

# === PERTURBED SAMPLING APPROACH ===
# Sample from original param sets and add small perturbations

# Define parameter bounds for perturbation clipping
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
  active_age_limit = c(3, 30),
  treg_discrimination_efficiency = c(0, 1)
)

# Function to add perturbations
add_perturbation = function(df, param_cols, perturbation_sd = 0.05, n_copies = 5) {
  # For each original param set, create n_copies perturbed versions
  all_perturbed = list()

  for (i in 1:nrow(df)) {
    original_row = df[i, param_cols]

    for (j in 1:n_copies) {
      perturbed = original_row

      for (param in param_cols) {
        original_val = original_row[[param]]
        bounds = param_bounds[[param]]

        # Calculate perturbation relative to parameter range
        param_range = bounds[2] - bounds[1]
        noise = rnorm(1, mean = 0, sd = perturbation_sd * param_range)

        # Add noise and clip to bounds
        new_val = original_val + noise
        new_val = max(bounds[1], min(bounds[2], new_val))

        # Special handling for integer parameter
        if (param == "active_age_limit") {
          new_val = round(new_val)
        }

        perturbed[[param]] = new_val
      }

      all_perturbed[[length(all_perturbed) + 1]] = perturbed
    }
  }

  return(bind_rows(all_perturbed))
}

set.seed(123)

# Separate by outcome
params_better = df_clustering %>% filter(dominant_outcome == "better")
params_worse = df_clustering %>% filter(dominant_outcome == "worse")
params_drift = df_clustering %>% filter(dominant_outcome == "drift")

cat("\n=== GENERATING PERTURBED SAMPLES ===\n")
cat("Original counts - Better:", nrow(params_better),
    "Worse:", nrow(params_worse),
    "Drift:", nrow(params_drift), "\n")

# Target: ~3333 of each
target_per_category = 3333

# Calculate how many originals to sample
n_copies_per_original = 5  # Each original will generate 5 perturbed copies
n_originals_needed = ceiling(target_per_category / n_copies_per_original)

cat("\nPerturbation settings:\n")
cat("  Copies per original:", n_copies_per_original, "\n")
cat("  Perturbation SD: 5% of parameter range\n")

# Sample originals with replacement if needed
sample_and_perturb = function(df, n_needed, param_cols, n_copies) {
  if (nrow(df) == 0) return(NULL)

  # Sample originals
  replace_flag = (n_needed > nrow(df))
  sampled_df = df[sample(1:nrow(df), size = n_needed, replace = replace_flag), ]

  # Add perturbations
  perturbed = add_perturbation(sampled_df, param_cols,
                               perturbation_sd = 0.05,
                               n_copies = n_copies)

  return(perturbed)
}

cat("\nGenerating better samples...\n")
balanced_better = sample_and_perturb(params_better, n_originals_needed, param_names, n_copies_per_original)
if (!is.null(balanced_better)) balanced_better$target_category = "better"

cat("Generating worse samples...\n")
balanced_worse = sample_and_perturb(params_worse, n_originals_needed, param_names, n_copies_per_original)
if (!is.null(balanced_worse)) balanced_worse$target_category = "worse"

cat("Generating drift samples...\n")
balanced_drift = sample_and_perturb(params_drift, n_originals_needed, param_names, n_copies_per_original)
if (!is.null(balanced_drift)) balanced_drift$target_category = "drift"

# Combine and trim to exact target
balanced_lhs_dataset = bind_rows(balanced_better, balanced_worse, balanced_drift)

# Shuffle
balanced_lhs_dataset = balanced_lhs_dataset %>%
  slice_sample(n = nrow(.))

balanced_lhs_dataset$param_set_id = 0:(nrow(balanced_lhs_dataset)-1)

cat("\n=== FINAL BALANCED DATASET ===\n")
cat("Total parameter sets:", nrow(balanced_lhs_dataset), "\n")
print(table(balanced_lhs_dataset$target_category))
cat("\nProportions:\n")
print(100 * prop.table(table(balanced_lhs_dataset$target_category)))

# Save WITH target_category
write.csv(balanced_lhs_dataset, "balanced_lhs_parameters_PERTURBED_with_target_cat.csv", row.names = FALSE)
# Save WITHOUT target_category
balanced_lhs_dataset = balanced_lhs_dataset[c('param_set_id',param_names)]
write.csv(balanced_lhs_dataset, "balanced_lhs_parameters_PERTURBED.csv", row.names = FALSE)

