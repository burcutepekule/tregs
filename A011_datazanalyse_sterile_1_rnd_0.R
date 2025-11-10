rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(gridExtra)

# ============================================================================
# LOAD PARAMETERS
# ============================================================================
params_df = read.csv("./original_lhs_parameters.csv", stringsAsFactors = FALSE)

# ============================================================================
# LOAD SIMULATION RESULTS
# ============================================================================

# Function to classify outcomes
classify_outcome <- function(df) {
  # Get final timepoint for each replicate
  final_state <- df %>%
    filter(t == max(t)) %>%
    group_by(param_set_id, rep_id, tregs_on) %>%
    summarise(
      final_healthy = mean(epithelial_healthy),
      final_injured = mean(epithelial_inj_1 + epithelial_inj_2 + epithelial_inj_3 +
                             epithelial_inj_4 + epithelial_inj_5),
      mean_injured_last_100 = NA,  # Will calculate below
      .groups = 'drop'
    )

  # Calculate mean injury over last 100 timesteps
  last_100_state <- df %>%
    filter(t > (max(t) - 100)) %>%
    group_by(param_set_id, rep_id, tregs_on) %>%
    summarise(
      mean_injured_last_100 = mean(epithelial_inj_1 + epithelial_inj_2 + epithelial_inj_3 +
                                      epithelial_inj_4 + epithelial_inj_5),
      .groups = 'drop'
    )

  final_state$mean_injured_last_100 <- last_100_state$mean_injured_last_100

  # Classify: chronic if mean injury in last 100 timesteps > threshold
  final_state$outcome <- ifelse(final_state$mean_injured_last_100 > 50, "chronic", "resolved")

  return(final_state)
}

# Load results for all parameter sets
all_outcomes <- data.frame()

for(param_set_id in 0:9999){
  tryCatch({
    df_treg_0 <- readRDS(paste0('./mass_sim_results_R/longitudinal_df_param_set_id_',
                                 param_set_id,'_tregs_0.rds'))
    df_treg_1 <- readRDS(paste0('./mass_sim_results_R/longitudinal_df_param_set_id_',
                                 param_set_id,'_tregs_1.rds'))

    # Classify outcomes
    outcomes_0 <- classify_outcome(df_treg_0)
    outcomes_1 <- classify_outcome(df_treg_1)

    all_outcomes <- rbind(all_outcomes, outcomes_0, outcomes_1)

  }, error = function(e) {
    # Skip if file doesn't exist
  })
}

# Merge with parameters
results_with_params <- all_outcomes %>%
  left_join(params_df, by = "param_set_id")

# ============================================================================
# CALCULATE INFLAMMATORY INDICES
# ============================================================================

# 1. Signal potential (accumulation / activation threshold)
results_with_params <- results_with_params %>%
  mutate(
    # DAMP potential
    DAMP_potential = (add_DAMPs / DAMPs_decay) / activation_threshold_DAMPs,

    # SAMP potential
    SAMP_potential = (add_SAMPs / SAMPs_decay) / activation_threshold_SAMPs,

    # Signal ratio (primary metric)
    Signal_Ratio = DAMP_potential / SAMP_potential,

    # ROS damage potential
    ROS_potential = (activity_ROS_M1_baseline / ros_decay) / th_ROS_epith_recover,

    # Feedback strength
    Feedback_Strength = rate_leak_commensal_injury / epith_recovery_chance,

    # Combined chronicity risk
    Chronicity_Risk = Signal_Ratio * log(1 + Feedback_Strength),

    # Simple inflammatory index
    Inflammatory_Index = (add_DAMPs / DAMPs_decay / activation_threshold_DAMPs) /
                        (add_SAMPs / SAMPs_decay / activation_threshold_SAMPs),

    # Extended index with diffusion
    Inflammatory_Index_Extended = Inflammatory_Index *
                                  sqrt(diffusion_speed_SAMPs / diffusion_speed_DAMPs),

    # Comprehensive instability score
    Instability_Score = Signal_Ratio *
                       ROS_potential *
                       (1 / epith_recovery_chance) *
                       rate_leak_commensal_injury *
                       exp(-2 * treg_discrimination_efficiency)
  )

# ============================================================================
# VISUALIZATION 1: HEATMAP - Signal Ratio vs Feedback Strength
# ============================================================================

# Aggregate outcomes per parameter set
param_outcomes <- results_with_params %>%
  group_by(param_set_id, tregs_on) %>%
  summarise(
    chronic_fraction = mean(outcome == "chronic"),
    Signal_Ratio = mean(Signal_Ratio),
    Feedback_Strength = mean(Feedback_Strength),
    Chronicity_Risk = mean(Chronicity_Risk),
    Inflammatory_Index = mean(Inflammatory_Index),
    DAMP_potential = mean(DAMP_potential),
    SAMP_potential = mean(SAMP_potential),
    .groups = 'drop'
  )

# Create bins for heatmap
param_outcomes <- param_outcomes %>%
  mutate(
    Signal_Ratio_bin = cut(Signal_Ratio, breaks = 20),
    Feedback_Strength_bin = cut(Feedback_Strength, breaks = 20)
  )

# Heatmap aggregation
heatmap_data <- param_outcomes %>%
  group_by(Signal_Ratio_bin, Feedback_Strength_bin, tregs_on) %>%
  summarise(
    mean_chronic = mean(chronic_fraction),
    n = n(),
    .groups = 'drop'
  ) %>%
  filter(n >= 3)  # Only show bins with at least 3 parameter sets

# Plot heatmap
p1 <- ggplot(heatmap_data, aes(x = Signal_Ratio_bin, y = Feedback_Strength_bin, fill = mean_chronic)) +
  geom_tile() +
  scale_fill_viridis(name = "Chronic\nFraction", limits = c(0, 1)) +
  facet_wrap(~tregs_on, labeller = label_both) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6)) +
  labs(title = "Phase Diagram: Signal Ratio vs Feedback Strength",
       x = "Signal Ratio (DAMP/SAMP potential)",
       y = "Feedback Strength (leakage/healing)")

print(p1)
ggsave("./figures/A011_phase_diagram_heatmap.png", p1, width = 12, height = 6, dpi = 300)

# Scatter plot version for better visibility
p2 <- ggplot(param_outcomes, aes(x = Signal_Ratio, y = Feedback_Strength,
                                  color = chronic_fraction)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_viridis(name = "Chronic\nFraction", limits = c(0, 1)) +
  facet_wrap(~tregs_on, labeller = label_both) +
  theme_minimal() +
  labs(title = "Phase Diagram: Signal Ratio vs Feedback Strength (Scatter)",
       x = "Signal Ratio (DAMP/SAMP potential)",
       y = "Feedback Strength (leakage/healing)") +
  geom_hline(yintercept = median(param_outcomes$Feedback_Strength),
             linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = median(param_outcomes$Signal_Ratio),
             linetype = "dashed", alpha = 0.5)

print(p2)
ggsave("./figures/A011_phase_diagram_scatter.png", p2, width = 12, height = 6, dpi = 300)

# ============================================================================
# VISUALIZATION 2: Chronicity Risk Index
# ============================================================================

p3 <- ggplot(param_outcomes, aes(x = Chronicity_Risk, fill = chronic_fraction > 0.5)) +
  geom_histogram(bins = 50, alpha = 0.7) +
  facet_wrap(~tregs_on, labeller = label_both, scales = "free_y") +
  theme_minimal() +
  scale_fill_manual(values = c("FALSE" = "blue", "TRUE" = "red"),
                    name = "Majority\nChronic",
                    labels = c("FALSE" = "Resolved", "TRUE" = "Chronic")) +
  labs(title = "Distribution of Chronicity Risk Index",
       x = "Chronicity Risk = Signal_Ratio × log(1 + Feedback_Strength)",
       y = "Count")

print(p3)
ggsave("./figures/A011_chronicity_risk_distribution.png", p3, width = 10, height = 6, dpi = 300)

# ============================================================================
# VISUALIZATION 3: PCA ANALYSIS
# ============================================================================

# Prepare data for PCA (use unique parameter sets)
pca_data <- results_with_params %>%
  group_by(param_set_id, tregs_on) %>%
  summarise(
    chronic_fraction = mean(outcome == "chronic"),
    add_DAMPs_norm = mean(add_DAMPs / DAMPs_decay),
    add_SAMPs_norm = mean(add_SAMPs / SAMPs_decay),
    inv_threshold_DAMPs = mean(1 / activation_threshold_DAMPs),
    inv_threshold_SAMPs = mean(1 / activation_threshold_SAMPs),
    diffusion_speed_DAMPs = mean(diffusion_speed_DAMPs),
    diffusion_speed_SAMPs = mean(diffusion_speed_SAMPs),
    rate_leak_commensal_injury = mean(rate_leak_commensal_injury),
    inv_epith_recovery = mean(1 / epith_recovery_chance),
    ROS_production = mean(activity_ROS_M1_baseline / ros_decay),
    inv_ROS_threshold = mean(1 / th_ROS_epith_recover),
    DAMPs_decay = mean(DAMPs_decay),
    SAMPs_decay = mean(SAMPs_decay),
    treg_discrimination = mean(treg_discrimination_efficiency),
    .groups = 'drop'
  )

# Select parameters for PCA
pca_matrix <- pca_data %>%
  select(add_DAMPs_norm, add_SAMPs_norm, inv_threshold_DAMPs, inv_threshold_SAMPs,
         diffusion_speed_DAMPs, diffusion_speed_SAMPs,
         rate_leak_commensal_injury, inv_epith_recovery,
         ROS_production, inv_ROS_threshold, DAMPs_decay, SAMPs_decay,
         treg_discrimination) %>%
  as.matrix()

# Remove any rows with NA or Inf
valid_rows <- complete.cases(pca_matrix) & !apply(pca_matrix, 1, function(x) any(is.infinite(x)))
pca_matrix <- pca_matrix[valid_rows, ]
pca_data_valid <- pca_data[valid_rows, ]

# Perform PCA
pca_result <- prcomp(pca_matrix, scale. = TRUE, center = TRUE)

# Add PC scores to data
pca_data_valid$PC1 <- pca_result$x[, 1]
pca_data_valid$PC2 <- pca_result$x[, 2]
pca_data_valid$PC3 <- pca_result$x[, 3]

# Plot PCA biplot
p4 <- ggplot(pca_data_valid, aes(x = PC1, y = PC2, color = chronic_fraction)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_viridis(name = "Chronic\nFraction", limits = c(0, 1)) +
  facet_wrap(~tregs_on, labeller = label_both) +
  theme_minimal() +
  labs(title = "PCA: Parameter Space Colored by Chronic Outcome",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "% variance)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "% variance)"))

print(p4)
ggsave("./figures/A011_PCA_biplot.png", p4, width = 12, height = 6, dpi = 300)

# Plot PC1 vs PC3
p5 <- ggplot(pca_data_valid, aes(x = PC1, y = PC3, color = chronic_fraction)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_viridis(name = "Chronic\nFraction", limits = c(0, 1)) +
  facet_wrap(~tregs_on, labeller = label_both) +
  theme_minimal() +
  labs(title = "PCA: PC1 vs PC3",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "% variance)"),
       y = paste0("PC3 (", round(summary(pca_result)$importance[2, 3] * 100, 1), "% variance)"))

print(p5)
ggsave("./figures/A011_PCA_PC1_PC3.png", p5, width = 12, height = 6, dpi = 300)

# PCA loadings plot
loadings_df <- as.data.frame(pca_result$rotation[, 1:3])
loadings_df$variable <- rownames(loadings_df)

p6 <- ggplot(loadings_df, aes(x = PC1, y = PC2, label = variable)) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.3, "cm")), color = "red", alpha = 0.7) +
  geom_text(hjust = "outward", vjust = "outward", size = 3) +
  theme_minimal() +
  labs(title = "PCA Loadings: Parameter Contributions",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "% variance)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "% variance)"))

print(p6)
ggsave("./figures/A011_PCA_loadings.png", p6, width = 10, height = 8, dpi = 300)

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

cat("\n========================================\n")
cat("SUMMARY: Inflammatory Index Analysis\n")
cat("========================================\n\n")

cat("Variance explained by PCs:\n")
print(summary(pca_result)$importance[, 1:5])

cat("\n\nTop 5 parameters contributing to PC1 (inflammatory axis):\n")
pc1_loadings <- sort(abs(pca_result$rotation[, 1]), decreasing = TRUE)[1:5]
print(pc1_loadings)

cat("\n\nCorrelation of indices with chronic outcome:\n")
correlations <- pca_data_valid %>%
  summarise(
    cor_PC1 = cor(PC1, chronic_fraction),
    cor_Signal_Ratio = cor(add_DAMPs_norm / add_SAMPs_norm, chronic_fraction),
    cor_Feedback = cor(rate_leak_commensal_injury * inv_epith_recovery, chronic_fraction),
    cor_Treg = cor(treg_discrimination, chronic_fraction)
  )
print(correlations)

cat("\n\nMedian values separating chronic vs resolved:\n")
comparison <- pca_data_valid %>%
  mutate(outcome_binary = ifelse(chronic_fraction > 0.5, "chronic", "resolved")) %>%
  group_by(outcome_binary, tregs_on) %>%
  summarise(
    median_PC1 = median(PC1),
    median_DAMP_SAMP_ratio = median(add_DAMPs_norm / add_SAMPs_norm),
    median_leakage = median(rate_leak_commensal_injury),
    median_healing = median(1 / inv_epith_recovery),
    .groups = 'drop'
  )
print(comparison)

# ============================================================================
# SAVE RESULTS
# ============================================================================

saveRDS(results_with_params, "./analysis_results/A011_results_with_inflammatory_indices.rds")
saveRDS(pca_result, "./analysis_results/A011_PCA_result.rds")
saveRDS(pca_data_valid, "./analysis_results/A011_PCA_data.rds")

cat("\n\nAnalysis complete! Figures saved to ./figures/\n")
cat("Data saved to ./analysis_results/\n")
