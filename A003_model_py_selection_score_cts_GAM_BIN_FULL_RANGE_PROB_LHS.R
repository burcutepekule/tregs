rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)  # For read_csv
library(stringr)
library(zoo)
library(mgcv)

source("/Users/burcutepekule/Dropbox/Treg_problem_v2/MISC/PLOT_FUNCTIONS.R")
df_params  = read_csv('/Users/burcutepekule/Desktop/tregs/mass_sim_results_full_range_LHS/sampled_parameters.csv', show_col_types = FALSE)
df_raw     = readRDS('/Users/burcutepekule/Desktop/tregs/df_summary_selection_score_raw_full_range_LHS.rds')
df_raw$tol = 25*0.05

df_model    = df_raw %>% dplyr::filter(comparison=='Treg_OFF_ON' & injury_type=='sterile')
df_model    = df_model %>% dplyr::mutate(nonzero=ifelse(effect_size %in% c('Medium','Large') & abs(mean_diff)>tol, 1, 0))
# df_model    = df_model %>% dplyr::mutate(nonzero=ifelse(abs(mean_diff)>0, 1, 0))

df_model    = merge(df_model, df_params, by='param_set_id')

#-- to check
df_model_nz = df_model %>% dplyr::filter(nonzero==1)
df_model_z  = df_model %>% dplyr::filter(nonzero==0)
dim(df_model_nz)[1]/dim(df_model_z)[1] # ~1?

# All parameter names (exactly as in your model)
param_names <- c(
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

# Output folder for saving plots
outdir <- "/Users/burcutepekule/Desktop/tregs/mass_sim_results/gam_histograms"
dir.create(outdir, showWarnings = FALSE)

# Function to create one histogram plot for a given parameter
plot_param_distribution <- function(param) {
  ggplot(df_model, aes(x = .data[[param]], fill = factor(nonzero))) +
    geom_density(alpha = 0.4, color = NA) +
    scale_fill_manual(values = c("0" = "gray70", "1" = "#3182bd"),
                      labels = c("Zero outcome", "Nonzero outcome")) +
    labs(
      x = param,
      y = "Density",
      fill = "Outcome",
      title = paste("Distribution of", param, "by outcome (nonzero)")
    ) +
    theme_minimal(base_size = 13)
}

# # Loop through all parameters, plot and save
# walk(param_names, function(param) {
#   p <- plot_param_distribution(param)
#   ggsave(
#     filename = file.path(outdir, paste0(param, "_histogram.png")),
#     plot = p,
#     width = 6,
#     height = 4,
#     dpi = 300
#   )
#   message("✅ Saved: ", param)
# })

table(df_model$nonzero)

df_agg = df_model %>%
  group_by(param_set_id) %>%
  summarise(
    nonzero_n = sum(nonzero == 1),
    zero_n    = sum(nonzero == 0),
    across(where(is.numeric) & !any_of(c("nonzero")), ~ first(.x)),
    .groups="drop"
  )

# This naturally handles different numbers of replicates per param set.
# If there’s extra-binomial variation, try family = quasibinomial as a robustness check.

# models the probability that a given parameter set produces a nonzero outcome across its replicates.
gam_binom_counts <- gam(
  cbind(nonzero_n, zero_n) ~
    s(th_ROS_microbe, bs="cs") +
    s(th_ROS_epith_recover, bs="cs") +
    s(epith_recovery_chance, bs="cs") +
    s(rat_com_pat_threshold, bs="cs") +
    s(diffusion_speed_DAMPs, bs="cs") +
    s(diffusion_speed_SAMPs, bs="cs") +
    s(diffusion_speed_ROS, bs="cs") +
    s(add_ROS, bs="cs") + s(add_DAMPs, bs="cs") + s(add_SAMPs, bs="cs") +
    s(ros_decay, bs="cs") + s(DAMPs_decay, bs="cs") + s(SAMPs_decay, bs="cs") +
    s(activation_threshold_DAMPs, bs="cs") + s(activation_threshold_SAMPs, bs="cs") +
    s(activity_engulf_M0_baseline, bs="cs") + s(activity_engulf_M1_baseline, bs="cs") +
    s(activity_engulf_M2_baseline, bs="cs") + s(activity_ROS_M1_baseline, bs="cs") +
    s(rate_leak_commensal_injury, bs="cs") + s(rate_leak_pathogen_injury, bs="cs") +
    s(rate_leak_commensal_baseline, bs="cs") + s(active_age_limit, bs="cs") +
    s(treg_discrimination_efficiency, bs="cs"),
  data   = df_agg,
  family = binomial,
  method = "REML"
)

# gam_binom_counts = readRDS('gam_binom_counts_full_range_RND.rds')

summary(gam_binom_counts)      # see which predictors matter
plot(gam_binom_counts, pages=1, shade=TRUE)   # visualize smooth effects
gam.check(gam_binom_counts)    # residual diagnostics

dev.off()
plot(gam_binom_counts, select = 3, shade = TRUE, seWithMean = TRUE, rug = TRUE)

# 🧩 Goal: You want to sample new parameter sets more intelligently —
# not purely random, but biased toward regions where the predicted activation probability ≈ 0.5
# That way, you’ll: 
# Avoid wasting simulations in regions where everything collapses (p ≈ 0) or always activates (p ≈ 1).
# Generate a balanced dataset (≈ equal 0s and 1s).
# Explore the model’s critical thresholds more efficiently.

# find approximate thresholds where the model predicts a 50% 
# activation probability for each significant parameter.
# For instance, for epithelial recovery chance:
  
# Create a fine grid for one parameter while holding others at mean

# -------------------------------------------------------------------------------------
grid <- df_agg %>%
  summarise(across(where(is.numeric), \(x) median(x, na.rm = TRUE))) %>%
  slice(rep(1, 500)) %>%
  mutate(epith_recovery_chance = seq(0, 1, length.out = 500))

grid$pred_p <- predict(gam_binom_counts, grid, type = "response")

plot(grid$epith_recovery_chance, grid$pred_p, type = "l")
abline(h = 0.5, col = "red", lty = 2)


# -------------------------------------------------------------------------------------
# df_model has columns: nonzero (0/1) + numeric parameters
num_cols <- param_names
num_cols <- setdiff(num_cols, c("nonzero"))  # exclude outcome

# alpha   <- 0.10         # keep middle 80% of the nonzero distribution
alpha   <- 0.25         # keep middle 50% of the nonzero distribution
margin  <- 0.05         # widen each side by 5% of full range

full_min <- sapply(df_model[num_cols], min, na.rm=TRUE)
full_max <- sapply(df_model[num_cols], max, na.rm=TRUE)

nz <- df_model[df_model$nonzero==1, num_cols, drop=FALSE]

nz_q_lo <- sapply(nz, \(x) quantile(x, alpha,     na.rm=TRUE))
nz_q_hi <- sapply(nz, \(x) quantile(x, 1 - alpha, na.rm=TRUE))

# widen a bit, then clip to global min/max
new_lo <- pmax(full_min, nz_q_lo - margin*(full_max-full_min))
new_hi <- pmin(full_max, nz_q_hi + margin*(full_max-full_min))

new_bounds <- data.frame(var=num_cols, new_lo, new_hi, full_min, full_max)
new_bounds

# Loop through all parameters, plot and save
walk(param_names, function(param) {
  p = plot_param_distribution(param)
  new_bounds_var = new_bounds %>% filter(var==param)
  lo = new_bounds_var$new_lo
  hi = new_bounds_var$new_hi
  p  = p +
    geom_vline(xintercept = lo, color = "red", linetype = "dashed", linewidth = 0.6) +
    geom_vline(xintercept = hi, color = "red", linetype = "dashed", linewidth = 0.6)
  
  ggsave(
    filename = file.path(outdir, paste0(param, "_histogram.png")),
    plot = p,
    width = 6,
    height = 4,
    dpi = 300
  )
  message("✅ Saved with bounds: ", param, " [", round(lo,3), ", ", round(hi,3), "]")
})


new_bounds_use = new_bounds[,1:3]
new_bounds_use$new_hi = round(new_bounds_use$new_hi , 2)
new_bounds_use$new_lo = round(new_bounds_use$new_lo , 2)
new_bounds_use_constraint = new_bounds_use %>% filter(var %in% c('th_ROS_epith_recover',
                                                                 'epith_recovery_chance',
                                                                 'ros_decay','DAMPs_decay'))

############

library(dplyr)
library(purrr)

num_cols <- setdiff(param_names, "nonzero")

# function to find the best quarter of the range for nonzero==1
find_best_quarter <- function(param, outcome, bins = 50) {
  df <- tibble(param = param, nonzero = outcome) %>%
    mutate(bin = ntile(param, bins)) %>%
    group_by(bin) %>%
    summarise(
      prob_nonzero = mean(nonzero, na.rm = TRUE),
      lo = min(param, na.rm = TRUE),
      hi = max(param, na.rm = TRUE),
      .groups = "drop"
    )
  
  # slide over consecutive bins to get 25% of total bins
  win_size <- ceiling(bins * 0.25)
  if (nrow(df) <= win_size) return(c(min(param), max(param)))
  
  # compute rolling mean of nonzero probability
  roll_means <- zoo::rollapply(df$prob_nonzero, win_size, mean, align = "left", fill = NA)
  best_start <- which.max(roll_means)
  
  best_lo <- df$lo[best_start]
  best_hi <- df$hi[min(best_start + win_size - 1, nrow(df))]
  c(best_lo, best_hi)
}

# compute new bounds for each numeric parameter
new_bounds <- map_dfr(num_cols, function(v) {
  best_range <- find_best_quarter(df_model[[v]], df_model$nonzero)
  full_min <- min(df_model[[v]], na.rm = TRUE)
  full_max <- max(df_model[[v]], na.rm = TRUE)
  tibble(var = v, new_lo = best_range[1], new_hi = best_range[2], full_min, full_max)
})

new_bounds
# Loop through all parameters, plot and save
walk(param_names, function(param) {
  p = plot_param_distribution(param)
  new_bounds_var = new_bounds %>% filter(var==param)
  lo = new_bounds_var$new_lo
  hi = new_bounds_var$new_hi
  p  = p +
    geom_vline(xintercept = lo, color = "red", linetype = "dashed", linewidth = 0.6) +
    geom_vline(xintercept = hi, color = "red", linetype = "dashed", linewidth = 0.6)
  
  ggsave(
    filename = file.path(outdir, paste0(param, "_histogram.png")),
    plot = p,
    width = 6,
    height = 4,
    dpi = 300
  )
  message("✅ Saved with bounds: ", param, " [", round(lo,3), ", ", round(hi,3), "]")
})


new_bounds_use = new_bounds[,1:3]
new_bounds_use$new_hi = round(new_bounds_use$new_hi , 2)
new_bounds_use$new_lo = round(new_bounds_use$new_lo , 2)
new_bounds_use_constraint = new_bounds_use %>% filter(var %in% c(
                                                                 'epith_recovery_chance',
                                                                 'ros_decay','DAMPs_decay'))

