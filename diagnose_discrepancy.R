rm(list=ls())
library(dplyr)
library(tidyr)
library(readr)

# Load the original decision tree results
df_raw_original = readRDS('all_comparison_results_0.rds')
df_params_original = read_csv('mass_sim_results/sampled_parameters.csv', show_col_types = FALSE)

# Load the presampled results
df_raw_presampled = readRDS('all_comparison_results_0_pre.rds')
df_params_presampled = read_csv('mass_sim_results_presampled/loaded_parameters.csv', show_col_types = FALSE)
balanced_params = read_csv('balanced_lhs_parameters.csv', show_col_types = FALSE)

# Settings from the original script
inj_type = 'sterile'
ss_start_threshold = 450
tol_in = 25 * 0.25  # 6.25

# Process original data
df_original = df_raw_original %>%
  filter(comparison=='Treg_OFF_ON' & injury_type==inj_type)

param_id_all_below_orig = df_original %>%
  group_by(param_set_id) %>%
  summarise(all_below = all(ss_start < ss_start_threshold), .groups = "drop") %>%
  filter(all_below) %>%
  pull(param_set_id)

df_original = df_original %>%
  filter(param_set_id %in% param_id_all_below_orig) %>%
  inner_join(df_params_original, by='param_set_id')

# Process presampled data
df_presampled = df_raw_presampled %>%
  filter(comparison=='Treg_OFF_ON' & injury_type==inj_type)

param_id_all_below_pre = df_presampled %>%
  group_by(param_set_id) %>%
  summarise(all_below = all(ss_start < ss_start_threshold), .groups = "drop") %>%
  filter(all_below) %>%
  pull(param_set_id)

df_presampled = df_presampled %>%
  filter(param_set_id %in% param_id_all_below_pre) %>%
  inner_join(df_params_presampled, by='param_set_id')

# Calculate effect sizes and classifications for both
classify_outcomes = function(df, tol) {
  df %>%
    mutate(
      abs_cohens_d = abs(cohens_d),
      effect_size = case_when(
        abs_cohens_d < 0.2 ~ "Negligible",
        abs_cohens_d < 0.5 & abs_cohens_d >= 0.2 ~ "Small",
        abs_cohens_d < 0.8 & abs_cohens_d >= 0.5 ~ "Medium",
        TRUE ~ "Large"
      ),
      outcome_class = case_when(
        effect_size %in% c('Large','Medium') & mean_diff > tol ~ "better",
        effect_size %in% c('Large','Medium') & mean_diff < -1*tol ~ "worse",
        TRUE ~ "drift"
      )
    )
}

df_original_classified = classify_outcomes(df_original, tol_in)
df_presampled_classified = classify_outcomes(df_presampled, tol_in)

# Summary statistics
cat("\n=== ORIGINAL DATA (from mass_simulation_LHS.py) ===\n")
cat("Total observations:", nrow(df_original_classified), "\n")
cat("Unique param sets:", n_distinct(df_original_classified$param_set_id), "\n\n")

orig_summary = df_original_classified %>%
  group_by(outcome_class) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(pct = 100 * n / sum(n))
print(orig_summary)

cat("\n=== PRESAMPLED DATA (from mass_simulation_LHS_presampled.py) ===\n")
cat("Total observations:", nrow(df_presampled_classified), "\n")
cat("Unique param sets:", n_distinct(df_presampled_classified$param_set_id), "\n\n")

pre_summary = df_presampled_classified %>%
  group_by(outcome_class) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(pct = 100 * n / sum(n))
print(pre_summary)

# Compare distributions of key metrics
cat("\n=== COMPARISON OF KEY METRICS ===\n")

cat("\nMean diff distribution:\n")
cat("Original - Mean:", mean(df_original_classified$mean_diff),
    "SD:", sd(df_original_classified$mean_diff), "\n")
cat("Presampled - Mean:", mean(df_presampled_classified$mean_diff),
    "SD:", sd(df_presampled_classified$mean_diff), "\n")

cat("\nCohen's d distribution:\n")
cat("Original - Mean:", mean(abs(df_original_classified$cohens_d)),
    "SD:", sd(abs(df_original_classified$cohens_d)), "\n")
cat("Presampled - Mean:", mean(abs(df_presampled_classified$cohens_d)),
    "SD:", sd(abs(df_presampled_classified$cohens_d)), "\n")

cat("\nEffect size distribution:\n")
print(table(df_original_classified$effect_size) / nrow(df_original_classified))
cat("\n")
print(table(df_presampled_classified$effect_size) / nrow(df_presampled_classified))

# Check if balanced_params has target_category info
if("target_category" %in% names(balanced_params)) {
  cat("\n=== EXPECTED vs ACTUAL (from balanced_lhs_parameters.csv) ===\n")
  expected = balanced_params %>%
    count(target_category) %>%
    mutate(pct = 100 * n / sum(n))
  print(expected)

  # Match param sets with their expected categories
  df_with_expected = df_presampled_classified %>%
    left_join(balanced_params %>% select(param_set_id, target_category),
              by = "param_set_id")

  if(nrow(df_with_expected) > 0 && "target_category" %in% names(df_with_expected)) {
    cat("\n=== CONFUSION MATRIX: Expected vs Actual ===\n")
    confusion = df_with_expected %>%
      group_by(target_category, outcome_class) %>%
      summarise(n = n(), .groups = "drop") %>%
      pivot_wider(names_from = outcome_class, values_from = n, values_fill = 0)
    print(confusion)

    # Calculate accuracy by category
    cat("\n=== PREDICTION ACCURACY BY CATEGORY ===\n")
    accuracy = df_with_expected %>%
      mutate(correct = (target_category == outcome_class)) %>%
      group_by(target_category) %>%
      summarise(
        total = n(),
        correct = sum(correct),
        accuracy = 100 * sum(correct) / n(),
        .groups = "drop"
      )
    print(accuracy)
  }
}

# Save diagnostic results
cat("\n=== Saving diagnostic plots ===\n")
pdf("diagnostic_comparison.pdf", width=12, height=8)

par(mfrow=c(2,2))

# Mean diff comparison
hist(df_original_classified$mean_diff, breaks=30, col=rgb(0,0,1,0.5),
     main="Mean Diff Distribution", xlab="mean_diff", xlim=range(c(df_original_classified$mean_diff, df_presampled_classified$mean_diff)))
hist(df_presampled_classified$mean_diff, breaks=30, col=rgb(1,0,0,0.5), add=TRUE)
legend("topright", c("Original", "Presampled"), fill=c(rgb(0,0,1,0.5), rgb(1,0,0,0.5)))
abline(v=c(-tol_in, tol_in), lty=2, col="black", lwd=2)

# Cohen's d comparison
hist(abs(df_original_classified$cohens_d), breaks=30, col=rgb(0,0,1,0.5),
     main="Abs Cohen's d Distribution", xlab="|Cohen's d|", xlim=range(c(abs(df_original_classified$cohens_d), abs(df_presampled_classified$cohens_d))))
hist(abs(df_presampled_classified$cohens_d), breaks=30, col=rgb(1,0,0,0.5), add=TRUE)
legend("topright", c("Original", "Presampled"), fill=c(rgb(0,0,1,0.5), rgb(1,0,0,0.5)))
abline(v=c(0.2, 0.5, 0.8), lty=2, col="black")

# Outcome proportions
barplot(rbind(orig_summary$pct, pre_summary$pct),
        beside=TRUE, names.arg=orig_summary$outcome_class,
        col=c("blue", "red"), main="Outcome Proportions (%)",
        legend.text=c("Original", "Presampled"))

# Mean diff vs Cohen's d
plot(df_original_classified$mean_diff, abs(df_original_classified$cohens_d),
     col=rgb(0,0,1,0.3), pch=16, main="Mean Diff vs |Cohen's d|",
     xlab="mean_diff", ylab="|Cohen's d|",
     xlim=range(c(df_original_classified$mean_diff, df_presampled_classified$mean_diff)),
     ylim=range(c(abs(df_original_classified$cohens_d), abs(df_presampled_classified$cohens_d))))
points(df_presampled_classified$mean_diff, abs(df_presampled_classified$cohens_d),
       col=rgb(1,0,0,0.3), pch=16)
abline(h=0.5, lty=2, col="black")
abline(v=c(-tol_in, tol_in), lty=2, col="black")
legend("topright", c("Original", "Presampled"), fill=c(rgb(0,0,1,0.5), rgb(1,0,0,0.5)))

dev.off()

cat("\nDiagnostic complete! Check diagnostic_comparison.pdf\n")
