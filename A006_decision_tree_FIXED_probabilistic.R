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

# [Same preprocessing as original...]
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

# Create outcome proportions - KEEP THE FRACTIONS, don't discretize yet!
df_clustering = df_summary %>%
  mutate(
    total = n_better + n_drift + n_worse,
    frac_better = n_better / total,
    frac_drift = n_drift / total,
    frac_worse = n_worse / total,
    # Also keep dominant for tree building
    dominant_outcome = case_when(
      frac_better > frac_worse & frac_better > frac_drift ~ "better",
      frac_worse > frac_better & frac_worse > frac_drift ~ "worse",
      TRUE ~ "drift"
    ),
    dominant_outcome = factor(dominant_outcome, levels = c("better", "drift", "worse"))
  )

cat("\n=== PARAMETER SET OUTCOME DISTRIBUTIONS ===\n")
cat("Example parameter sets and their outcome probabilities:\n")
print(head(df_clustering %>% select(frac_better, frac_drift, frac_worse, dominant_outcome), 10))

# Build decision tree on dominant outcome
ctree_model = ctree(dominant_outcome ~ .,
                    data = df_clustering %>% select(dominant_outcome, all_of(param_names)),
                    control = ctree_control(
                      mincriterion = 0.90,
                      minsplit = 5,
                      minbucket = 5,
                      maxdepth = 5,
                      testtype = "Univariate"
                    ))

df_clustering = df_clustering %>% mutate(ctree_node = predict(ctree_model, type = "node"))

# === KEY DIFFERENCE: Use FRACTIONAL probabilities, not dominant outcome ===

# For each node, calculate the AVERAGE outcome probabilities across all param sets in that node
node_summary = df_clustering %>%
  group_by(ctree_node) %>%
  summarise(
    n_param_sets = n(),
    # Average probabilities across all param sets in this node
    prob_better = mean(frac_better),
    prob_worse = mean(frac_worse),
    prob_drift = mean(frac_drift),
    .groups = "drop"
  )

cat("\n=== NODE SUMMARY (Probabilistic) ===\n")
print(node_summary)

# Run simulated annealing to balance (same as original)
objective_function = function(random_integers, node_summary) {
  node_summary_temp = node_summary
  node_summary_temp$pop_size = round(random_integers)

  node_summary_temp = node_summary_temp %>%
    mutate(
      n_for_better = round(prob_better * pop_size),
      n_for_worse = round(prob_worse * pop_size),
      n_for_drift = round(prob_drift * pop_size)
    )

  all_sum = sum(node_summary_temp$n_for_better +
                  node_summary_temp$n_for_worse +
                  node_summary_temp$n_for_drift)

  better_rat = sum(node_summary_temp$n_for_better) / all_sum
  worse_rat = sum(node_summary_temp$n_for_worse) / all_sum
  drift_rat = sum(node_summary_temp$n_for_drift) / all_sum

  diff_rat = rep(1/3, 3) - c(better_rat, worse_rat, drift_rat)
  diff_rat_euc = sqrt(sum(diff_rat^2))

  return(diff_rat_euc)
}

M = 10000
set.seed(42)
current_solution = sample(1:M, size = nrow(node_summary), replace = TRUE)
current_score = objective_function(current_solution, node_summary)

best_solution = current_solution
best_score = current_score

temperature  = 100
cooling_rate = 0.9

for (iter in 1:100000) {
  new_solution = current_solution
  positions_to_change = sample(1:length(new_solution), size = sample(1:3, 1))
  new_solution[positions_to_change] = sample(0:M, size = length(positions_to_change), replace = TRUE)

  new_score = objective_function(new_solution, node_summary)

  if (new_score < current_score || runif(1) < exp((current_score - new_score) / temperature)) {
    current_solution = new_solution
    current_score = new_score

    if (current_score < best_score) {
      best_solution = current_solution
      best_score = current_score
    }
  }

  temperature = temperature * cooling_rate

  if (iter %% 1000 == 0) {
    print(paste("Iteration:", iter, "Best score:", best_score))
  }
}

node_allocations = node_summary
node_allocations$pop_size = best_solution
node_allocations = node_allocations %>%
  mutate(
    n_for_better = round(prob_better * pop_size),
    n_for_worse = round(prob_worse * pop_size),
    n_for_drift = round(prob_drift * pop_size)
  )

all_sum    = sum(node_allocations$n_for_better + node_allocations$n_for_worse + node_allocations$n_for_drift)
better_rat = sum(node_allocations$n_for_better) / all_sum
worse_rat  = sum(node_allocations$n_for_worse) / all_sum
drift_rat  = sum(node_allocations$n_for_drift) / all_sum

cat("\n=== EXPECTED OUTCOME PROPORTIONS ===\n")
print(100 * round(c(better_rat, worse_rat, drift_rat), 3))

# === NEW SAMPLING STRATEGY: Mixture-based sampling ===
# Instead of LHS within bounds, use a mixture of direct resampling and perturbation
# weighted by each parameter set's outcome probabilities

cat("\n=== GENERATING BALANCED DATASET (Mixture Approach) ===\n")

all_samples = list()

for (i in 1:nrow(node_allocations)) {
  node_id = node_allocations$ctree_node[i]
  n_better = node_allocations$n_for_better[i]
  n_worse = node_allocations$n_for_worse[i]
  n_drift = node_allocations$n_for_drift[i]

  # Get all parameter sets in this node
  node_param_sets = df_clustering %>% filter(ctree_node == node_id)

  if (nrow(node_param_sets) == 0) next

  # === Sample for BETTER outcomes ===
  if (n_better > 0) {
    # Sample parameter sets weighted by their frac_better
    weights = node_param_sets$frac_better
    if (sum(weights) > 0) {
      # Sample with replacement, weighted by probability of producing "better"
      sampled_indices = sample(1:nrow(node_param_sets),
                              size = n_better,
                              replace = TRUE,
                              prob = weights)

      samples = node_param_sets[sampled_indices, param_names]
      samples$expected_category = "better"
      samples$source_node = node_id
      samples$source_prob_better = node_param_sets$frac_better[sampled_indices]

      all_samples[[length(all_samples) + 1]] = samples
    }
  }

  # === Sample for WORSE outcomes ===
  if (n_worse > 0) {
    weights = node_param_sets$frac_worse
    if (sum(weights) > 0) {
      sampled_indices = sample(1:nrow(node_param_sets),
                              size = n_worse,
                              replace = TRUE,
                              prob = weights)

      samples = node_param_sets[sampled_indices, param_names]
      samples$expected_category = "worse"
      samples$source_node = node_id
      samples$source_prob_worse = node_param_sets$frac_worse[sampled_indices]

      all_samples[[length(all_samples) + 1]] = samples
    }
  }

  # === Sample for DRIFT outcomes ===
  if (n_drift > 0) {
    weights = node_param_sets$frac_drift
    if (sum(weights) > 0) {
      sampled_indices = sample(1:nrow(node_param_sets),
                              size = n_drift,
                              replace = TRUE,
                              prob = weights)

      samples = node_param_sets[sampled_indices, param_names]
      samples$expected_category = "drift"
      samples$source_node = node_id
      samples$source_prob_drift = node_param_sets$frac_drift[sampled_indices]

      all_samples[[length(all_samples) + 1]] = samples
    }
  }

  cat(sprintf("  Node %d: sampled %d sets (better:%d, worse:%d, drift:%d)\n",
              node_id, n_better + n_worse + n_drift, n_better, n_worse, n_drift))
}

# Combine all samples
balanced_lhs_dataset = bind_rows(all_samples)

# Shuffle
set.seed(123)
balanced_lhs_dataset = balanced_lhs_dataset %>%
  slice_sample(n = nrow(.))

balanced_lhs_dataset$param_set_id = 0:(nrow(balanced_lhs_dataset)-1)

cat("\n=== FINAL BALANCED DATASET ===\n")
cat("Total parameter sets:", nrow(balanced_lhs_dataset), "\n")
cat("\nExpected category distribution:\n")
print(table(balanced_lhs_dataset$expected_category))
cat("\nProportions:\n")
print(100 * prop.table(table(balanced_lhs_dataset$expected_category)))

# # Save with metadata
# write.csv(balanced_lhs_dataset, "balanced_lhs_parameters_PROBABILISTIC.csv", row.names = FALSE)

# Save WITH target_category
write.csv(balanced_lhs_dataset, "balanced_lhs_parameters_PERTURBED_with_target_cat.csv", row.names = FALSE)
# Save WITHOUT target_category
balanced_lhs_dataset = balanced_lhs_dataset[c('param_set_id',param_names)]
write.csv(balanced_lhs_dataset, "balanced_lhs_parameters_PERTURBED.csv", row.names = FALSE)

# cat("\n=== KEY INSIGHT ===\n")
# cat("This approach:\n")
# cat("1. Preserves the STOCHASTIC nature of outcomes\n")
# cat("2. Samples parameter sets WEIGHTED by their probability of producing each outcome\n")
# cat("3. Parameter sets with higher frac_worse are more likely to be sampled for 'worse'\n")
# cat("4. Expected category is probabilistic, not deterministic\n")
# cat("\nDataset saved to: balanced_lhs_parameters_PROBABILISTIC.csv\n")
