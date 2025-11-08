#!/usr/bin/env Rscript
# Quick validation script to compare original vs fixed approaches

library(dplyr)
library(readr)

cat("\n=== VALIDATION: Original vs Fixed Approach ===\n\n")

# Check if original balanced parameters exist
if (file.exists("balanced_lhs_parameters.csv")) {
  cat("ORIGINAL APPROACH (balanced_lhs_parameters.csv):\n")
  original_params = read_csv("balanced_lhs_parameters.csv", show_col_types = FALSE)
  cat("  Total parameter sets:", nrow(original_params), "\n")
  cat("  Columns:", paste(names(original_params), collapse=", "), "\n")

  if ("target_category" %in% names(original_params)) {
    cat("\n  Expected distribution:\n")
    print(table(original_params$target_category))
    cat("\n  Percentages:\n")
    print(100 * prop.table(table(original_params$target_category)))
  } else {
    cat("\n  ⚠️  WARNING: 'target_category' column NOT found!\n")
    cat("     This means you cannot validate expected vs actual outcomes.\n")
  }
} else {
  cat("Original balanced_lhs_parameters.csv not found.\n")
}

cat("\n" , rep("=", 60), "\n\n", sep="")

# Check if DIRECT resampling results exist
if (file.exists("balanced_lhs_parameters_DIRECT.csv")) {
  cat("FIXED APPROACH - DIRECT RESAMPLING:\n")
  direct_params = read_csv("balanced_lhs_parameters_DIRECT.csv", show_col_types = FALSE)
  cat("  Total parameter sets:", nrow(direct_params), "\n")

  if ("target_category" %in% names(direct_params)) {
    cat("\n  Expected distribution:\n")
    print(table(direct_params$target_category))
    cat("\n  Percentages:\n")
    print(100 * prop.table(table(direct_params$target_category)))
    cat("\n  ✓ target_category column IS present for validation\n")
  }
} else {
  cat("FIXED APPROACH - DIRECT RESAMPLING: Not generated yet\n")
  cat("  Run: Rscript A006_decision_tree_FIXED_direct_resample.R\n")
}

cat("\n" , rep("=", 60), "\n\n", sep="")

# Check if PERTURBED results exist
if (file.exists("balanced_lhs_parameters_PERTURBED.csv")) {
  cat("FIXED APPROACH - PERTURBED SAMPLING:\n")
  perturbed_params = read_csv("balanced_lhs_parameters_PERTURBED.csv", show_col_types = FALSE)
  cat("  Total parameter sets:", nrow(perturbed_params), "\n")

  if ("target_category" %in% names(perturbed_params)) {
    cat("\n  Expected distribution:\n")
    print(table(perturbed_params$target_category))
    cat("\n  Percentages:\n")
    print(100 * prop.table(table(perturbed_params$target_category)))
    cat("\n  ✓ target_category column IS present for validation\n")
  }
} else {
  cat("FIXED APPROACH - PERTURBED SAMPLING: Not generated yet\n")
  cat("  Run: Rscript A006_decision_tree_FIXED_perturbed.R\n")
}

cat("\n" , rep("=", 60), "\n\n", sep="")

# If actual results exist, compare expected vs actual
if (file.exists("all_comparison_results_0_pre.rds")) {
  cat("ACTUAL SIMULATION RESULTS (Presampled):\n")

  df_raw = readRDS('all_comparison_results_0_pre.rds')
  params_file = NULL

  # Try to load the parameters file used
  if (file.exists("mass_sim_results_presampled/loaded_parameters.csv")) {
    params_file = "mass_sim_results_presampled/loaded_parameters.csv"
  }

  if (!is.null(params_file)) {
    df_params = read_csv(params_file, show_col_types = FALSE)

    # Apply same filtering as original script
    inj_type = 'sterile'
    ss_start_threshold = 450
    tol_in = 25 * 0.25

    df_raw = df_raw %>% dplyr::filter(comparison=='Treg_OFF_ON' & injury_type==inj_type)

    param_id_all_below = df_raw %>%
      dplyr::group_by(param_set_id) %>%
      dplyr::summarise(all_below = all(ss_start < ss_start_threshold), .groups = "drop") %>%
      dplyr::filter(all_below) %>%
      dplyr::pull(param_set_id)

    df_raw = df_raw %>% dplyr::filter(param_set_id %in% param_id_all_below)

    df_raw = df_raw %>%
      dplyr::mutate(
        abs_cohens_d = abs(cohens_d),
        effect_size = case_when(
          abs_cohens_d < 0.2 ~ "Negligible",
          abs_cohens_d < 0.5 & abs_cohens_d >= 0.2 ~ "Small",
          abs_cohens_d < 0.8 & abs_cohens_d >= 0.5 ~ "Medium",
          TRUE ~ "Large"
        ),
        outcome_class = case_when(
          effect_size %in% c('Large','Medium') & mean_diff > tol_in ~ "better",
          effect_size %in% c('Large','Medium') & mean_diff < -1*tol_in ~ "worse",
          TRUE ~ "drift"
        )
      )

    cat("\n  Actual outcome distribution:\n")
    actual_counts = table(df_raw$outcome_class)
    print(actual_counts)
    cat("\n  Percentages:\n")
    actual_pct = 100 * prop.table(actual_counts)
    print(actual_pct)

    # If target_category exists in params, show confusion matrix
    if ("target_category" %in% names(df_params)) {
      cat("\n  CONFUSION MATRIX (Expected vs Actual):\n")
      df_with_expected = df_raw %>%
        left_join(df_params %>% select(param_set_id, target_category), by = "param_set_id")

      confusion = table(Expected = df_with_expected$target_category,
                       Actual = df_with_expected$outcome_class)
      print(confusion)

      # Calculate accuracy
      accuracy_total = sum(diag(confusion)) / sum(confusion)
      cat("\n  Overall accuracy:", round(100 * accuracy_total, 2), "%\n")

      # Accuracy by category
      cat("\n  Accuracy by expected category:\n")
      for (cat in rownames(confusion)) {
        cat_accuracy = confusion[cat, cat] / sum(confusion[cat, ])
        cat("   ", cat, ":", round(100 * cat_accuracy, 2), "%\n")
      }
    }
  }
} else {
  cat("ACTUAL SIMULATION RESULTS: Not available yet\n")
  cat("  Run simulations first, then process with A008 and A009 scripts\n")
}

cat("\n" , rep("=", 60), "\n\n", sep="")

cat("RECOMMENDATIONS:\n")
cat("1. Use the FIXED approach (direct resampling or perturbed)\n")
cat("2. Always keep target_category column for validation\n")
cat("3. After running simulations, compare expected vs actual using this script\n")
cat("4. If accuracy is low, investigate which parameters are causing misclassification\n\n")
