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
df_raw    = readRDS('/Users/burcutepekule/Desktop/tregs/all_comparison_results_0.rds')
df_params = read_csv('/Users/burcutepekule/Desktop/tregs/mass_sim_results/sampled_parameters.csv', show_col_types = FALSE)

t_max    = 500
inj_type = 'sterile'
df_model = df_raw %>% dplyr::filter(comparison=='Treg_OFF_ON' & injury_type==inj_type)
df_model = df_model %>% dplyr::mutate(log_integ_diff = sign(integ_diff)*log(1+abs(integ_diff))/(t_max-ss_start))
df_model = inner_join(df_model, df_params, by='param_set_id')

#----- filter based on ss_start, it cannot be too large otherwise not much to compare!
ss_start_threshold   = 450
param_id_all_below = df_model %>%
  dplyr::group_by(param_set_id) %>%
  dplyr::summarise(all_below = all(ss_start < ss_start_threshold), .groups = "drop") %>%
  dplyr::filter(all_below) %>%
  dplyr::pull(param_set_id)
df_model = df_model %>% dplyr::filter(param_set_id %in% param_id_all_below)

# #----- filter based on replicate_id, less than 10 means incomplete!
param_id_all_complete = df_model %>%
  dplyr::group_by(param_set_id, comparison, injury_type) %>%
  dplyr::summarise(all_complete = (n_distinct(replicate_id) == 10), .groups = "drop") %>%
  dplyr::group_by(param_set_id) %>%
  dplyr::summarise(all_complete = all(all_complete), .groups = "drop") %>%
  dplyr::filter(all_complete) %>%
  dplyr::pull(param_set_id)
df_model = df_model %>% dplyr::filter(param_set_id %in% param_id_all_complete)

#------- Effect size based on cohens d -----------------------------------------
df_model$tol = 25*0.25
df_model = df_model %>% dplyr::mutate(abs_cohens_d = abs(cohens_d))
# df_model = df_model %>% dplyr::mutate(effect_size = case_when(
#   abs_cohens_d < 0.2 ~ "Negligible",
#   abs_cohens_d < 0.5 & abs_cohens_d>= 0.2  ~ "Small",
#   abs_cohens_d < 0.8 & abs_cohens_d>= 0.5 ~ "Medium",
# ))

df_summary = df_model %>%
  dplyr::group_by(param_set_id, injury_type, comparison) %>%
  dplyr::summarise(
    n_better = sum(mean_diff > tol, na.rm = TRUE),
    n_drift  = sum((mean_diff <= tol & mean_diff >= -1*tol), na.rm = TRUE),
    n_worse  = sum(mean_diff < -1*tol, na.rm = TRUE),
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

# Create outcome proportions
df_clustering = df_summary %>%
  mutate(
    total = n_better + n_drift + n_worse,
    frac_better = n_better / total,
    frac_drift = n_drift / total,
    frac_worse = n_worse / total,
    frac_select = frac_better+frac_worse,
  )

# # Extract just the outcome proportions
# outcome_matrix = df_clustering %>%
#   select(frac_better, frac_drift, frac_worse) %>%
#   as.matrix()

# Extract just the outcome proportions
outcome_matrix = df_clustering %>%
  select(frac_select, frac_drift) %>%
  as.matrix()

# Compute distance matrix
dist_outcomes = dist(outcome_matrix, method = "euclidean")

# Hierarchical clustering
hc_outcomes = hclust(dist_outcomes, method = "ward.D2")

# Choose k based on the plots above (let's say k=4 for now)
k = 2 # ADJUST THIS based on your elbow/silhouette plots

# Cut the dendrogram
df_clustering$cluster = cutree(hc_outcomes, k = k)

# Examine cluster characteristics
cluster_summary = df_clustering %>%
  group_by(cluster) %>%
  summarise(
    n_param_sets = n(),
    
    # Mean fractions
    mean_frac_better = mean(frac_better),
    mean_frac_drift = mean(frac_drift),
    mean_frac_worse = mean(frac_worse),
    
    # Standard deviations
    sd_frac_better = sd(frac_better),
    sd_frac_worse = sd(frac_worse),
    
    # Count of any events
    n_with_any_better = sum(n_better > 0),
    n_with_any_worse = sum(n_worse > 0),
    
    .groups = "drop"
  ) %>%
  arrange(desc(mean_frac_better))

print("Cluster Summary:")
print(cluster_summary)

# 2D visualization: Better vs Worse
ggplot(df_clustering, aes(x = frac_better, y = frac_worse, 
                          color = factor(cluster))) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_brewer(palette = "Set1", name = "Cluster") +
  labs(title = "Outcome Clusters: Better vs Worse",
       subtitle = paste(k, "clusters identified"),
       x = "Fraction Better", 
       y = "Fraction Worse") +
  theme_minimal() +
  theme(legend.position = "right")

# Heatmap of cluster profiles
cluster_summary_long = cluster_summary %>%
  select(cluster, mean_frac_better, mean_frac_drift, mean_frac_worse) %>%
  pivot_longer(cols = starts_with("mean_frac_"), 
               names_to = "outcome", 
               values_to = "fraction") %>%
  mutate(outcome = gsub("mean_frac_", "", outcome))

ggplot(cluster_summary_long, aes(x = factor(cluster), y = outcome, fill = fraction)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f", fraction)), color = "white", size = 5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0.5, limits = c(0, 1)) +
  labs(title = "Cluster Outcome Profiles",
       subtitle = "Mean fraction for each outcome type",
       x = "Cluster", y = "Outcome Type") +
  theme_minimal()

#------- what parameters define cluster 2?

# Compare parameter distributions between clusters
param_comparison <- df_clustering %>%
  group_by(cluster) %>%
  summarise(
    across(all_of(param_names), 
           list(median = median, 
                q25 = ~quantile(., 0.25),
                q75 = ~quantile(., 0.75)),
           .names = "{.col}_{.fn}"),
    .groups = "drop"
  )

# Transpose for easier reading
param_comparison_t <- param_comparison %>%
  pivot_longer(-cluster, names_to = "param_stat", values_to = "value") %>%
  separate(param_stat, into = c("parameter", "statistic"), sep = "_(?=[^_]+$)") %>%
  pivot_wider(names_from = c(cluster, statistic), 
              values_from = value,
              names_sep = "_")

print("Parameter comparison between clusters:")
print(param_comparison_t)

# Statistical tests: which parameters differ significantly?
library(broom)

param_tests <- map_df(param_names, function(param) {
  cluster1_vals <- df_clustering %>% filter(cluster == 1) %>% pull(!!sym(param))
  cluster2_vals <- df_clustering %>% filter(cluster == 2) %>% pull(!!sym(param))
  
  # Wilcoxon test (non-parametric, good for your sample size imbalance)
  test_result <- wilcox.test(cluster2_vals, cluster1_vals)
  
  tibble(
    parameter = param,
    median_drift = median(cluster1_vals),
    median_selection = median(cluster2_vals),
    diff = median(cluster2_vals) - median(cluster1_vals),
    p_value = test_result$p.value
  )
}) %>%
  arrange(p_value)

print("\nParameters that distinguish 'Selection' from 'Drift' clusters:")
print(param_tests)

# Visualize top differences
top_params <- head(param_tests, 6)$parameter

plots <- lapply(top_params, function(param) {
  ggplot(df_clustering, aes(x = factor(cluster), y = !!sym(param), 
                            fill = factor(cluster))) +
    geom_boxplot() +
    scale_fill_manual(values = c("1" = "gray", "2" = "darkgreen"),
                      labels = c("Drift", "Selection")) +
    labs(title = param, x = "Cluster", y = "Value") +
    theme_minimal() +
    theme(legend.position = "none")
})

library(patchwork)
wrap_plots(plots, ncol = 3)