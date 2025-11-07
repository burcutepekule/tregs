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
ss_start_threshold   = 250
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

df_model$tol = 25*0.25
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

library(ggplot2)
library(ggdendro)
library(factoextra)
library(cluster)

# Create outcome proportions
df_clustering = df_summary %>%
  mutate(
    total = n_better + n_drift + n_worse,
    frac_better = div0(n_better,total),
    frac_drift = div0(n_drift,total),
    frac_worse = div0(n_worse,total)
  )

# Create a binary outcome: any selection vs pure drift
df_clustering = df_clustering %>%
  mutate(
    frac_select = frac_worse #or -, if truly looking for biological selection?
  )

# # ---------- Conditional Inference Trees (CTree): These use statistical tests for splitting and handle interactions better:

library(partykit)

ctree_model = ctree(frac_select ~ .,
                    data = df_clustering %>% select(frac_select, all_of(param_names)),
                    control = ctree_control(
                      # === Splitting criteria ===
                      mincriterion = 0.90,        # 1-p-value threshold (higher = stricter splits)
                      # 0.95 = p<0.05, 0.99 = p<0.01, 0.999 = p<0.001
                      
                      # === Node size controls ===
                      minsplit = 10,               # Min observations to attempt a split
                      minbucket = 10,              # Min observations in terminal node
                      maxdepth = 5,               # Max tree depth (0 = unlimited)
                      
                      # === Test type ===
                      testtype = "Univariate",    # Multiple testing correction
                      # Options: "Bonferroni", "MonteCarlo", "Univariate", "Teststatistic"
                      
                      # === If using MonteCarlo ===
                      # testtype = "MonteCarlo",
                      # nresample = 9999,          # Number of permutations
                      
                      # === Handling ties ===
                      splittry = 2L,              # Number of tries for binary splits (for ties)
                      
                      # === Multiway splits ===
                      splitflavour = "ctree"      # How to handle multiway splits
                      # Options: "ctree", "exhaustive"
                    ))

png('tree_worse.png', width=2800,height=1800,units = 'px')
plot(ctree_model, 
     main = "Conditional Inference Tree",
     type = "simple")
dev.off()

# ----------
df_clustering$ctree_node = predict(ctree_model, type = "node")

node_summary = df_clustering %>%
  group_by(ctree_node) %>%
  summarise(
    n = n(),
    mean_select = mean(frac_select),
    mean_better = mean(frac_better),
    mean_worse = mean(frac_worse),
    sd_select = sd(frac_select),
    .groupsx = "drop"
  ) %>%
  arrange(desc(mean_select))

opt_nodes     = as.numeric(node_summary[1:1,]$ctree_node)
df_clustering = df_clustering %>% dplyr::mutate(condition_label=ifelse(ctree_node %in% opt_nodes, 'optimal', 'not_optimal'))

# Alternative: Density plot (smoother)
ggplot(df_clustering, aes(x = frac_select, fill = condition_label, color = condition_label)) +
  geom_density(alpha = 0.5, linewidth = 1) +
  scale_fill_manual(values = c("optimal" = "darkgreen",
                               "not_optimal" = "gray"),
                    name = NULL) +
  scale_color_manual(values = c("optimal" = "darkgreen",
                                "not_optimal" = "gray"),
                     name = NULL) +
  labs(title = "Distribution of Selection Fraction",
       subtitle = "Optimal conditions",
       x = "Fraction Selected (Better)",
       y = "Density") +
  theme_minimal() +
  theme(legend.position = "top")

# ---- rules?

# Function to extract path/rules for a given node
get_node_rules = function(tree, node_id) {
  # Get the path from root to the target node
  path = partykit:::.list.rules.party(tree, i = node_id)
  return(path)
}
node_rules = get_node_rules(ctree_model, opt_nodes)
print(node_rules)

dim(df_clustering)[1]

print(opt_nodes)
df_clustering_opt = df_clustering %>% dplyr::filter(condition_label=='optimal')
100*mean(df_clustering_opt$frac_select)

df_clustering_nz = df_clustering %>% dplyr::filter(total>0)

ctree_model_worse=ctree_model
saveRDS(ctree_model_worse,'ctree_model_worse.rds')