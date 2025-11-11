# rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)  # For read_csv
library(stringr)
library(zoo)

source("/Users/burcutepekule/Dropbox/Treg_problem_v2/MISC/PLOT_FUNCTIONS.R")
inj_type             = 'sterile'
ss_start_threshold   = 4500
t_max                = 5000
tol_in               = 25*0.25

# df_raw    = readRDS('/Users/burcutepekule/Desktop/tregs/all_comparison_results_sterile_1_trnd_0_old.rds')
df_raw    = readRDS('/Users/burcutepekule/Desktop/tregs/all_comparison_results_sterile_1_trnd_0.rds')
length(unique(df_raw$param_set_id))
df_params = read_csv('/Users/burcutepekule/Desktop/tregs/original_lhs_parameters.csv', show_col_types = FALSE)

df_raw_keep = df_raw
df_raw = df_raw %>% dplyr::filter(comparison=='Treg_OFF_ON' & injury_type==inj_type)
hist(df_raw_keep$mean_diff)

#----- filter based on ss_start, it cannot be too large otherwise not much to compare!
param_id_all_below = df_raw %>%
  dplyr::group_by(param_set_id) %>%
  dplyr::summarise(all_below = all(ss_start < ss_start_threshold), .groups = "drop") %>%
  dplyr::filter(all_below) %>%
  dplyr::pull(param_set_id)
df_raw = df_raw %>% dplyr::filter(param_set_id %in% param_id_all_below)

#----- filter based on replicate_id, less than 10 means incomplete!
# param_id_all_complete = df_raw %>%
#   # dplyr::group_by(param_set_id, comparison, injury_type) %>%
#   dplyr::group_by(param_set_id, comparison) %>%
#   dplyr::summarise(all_complete = (n_distinct(replicate_id) == 10), .groups = "drop") %>%
#   dplyr::group_by(param_set_id) %>%
#   dplyr::summarise(all_complete = all(all_complete), .groups = "drop") %>%
#   dplyr::filter(all_complete) %>%
#   dplyr::pull(param_set_id)
# df_raw = df_raw %>% dplyr::filter(param_set_id %in% param_id_all_complete)
# length(unique(param_id_all_complete)) 

df_raw     = df_raw %>% dplyr::mutate(abs_cohens_d = abs(cohens_d))
df_raw     = df_raw %>% dplyr::mutate(effect_size = case_when(
  abs_cohens_d < 0.2 ~ "Negligible",
  abs_cohens_d < 0.5 & abs_cohens_d>= 0.2  ~ "Small",
  abs_cohens_d < 0.8 & abs_cohens_d>= 0.5 ~ "Medium",
  TRUE ~ "Large"
))

df_raw_plot_nz = df_raw %>% dplyr::filter(effect_size %in% c("Medium","Large") & (mean_diff>tol_in | mean_diff<(-1*tol_in)))
df_raw_plot_z  = df_raw %>% dplyr::filter((mean_diff<=tol_in & mean_diff>=(-1*tol_in))) # this is not true- what about effect_size small but abs(mean_diff)>tol_in?
df_raw_plot    = rbind(df_raw_plot_z, df_raw_plot_nz)
hist(df_raw_plot$mean_diff,30)

x   = df_raw_plot$mean_diff
round(100*sum(x>tol_in)/length(x),2)
round(100*sum(x<=tol_in & x>=(-1*tol_in))/length(x),2)
round(100*sum(x<(-1*tol_in))/length(x),2)
hist(x,30)

#-----------
var_df = df_raw_plot %>%
  group_by(param_set_id) %>%
  summarise(
    variance = var(mean_diff),
    sd = sd(mean_diff),
    mean = mean(mean_diff),
    min = min(mean_diff),
    max = max(mean_diff),
    abs_diff = abs(max-min),
    n_replicates = n()
  ) %>%
  arrange(desc(variance))


df_raw = df_raw %>% inner_join(var_df, by='param_set_id')
df_raw = df_raw %>% dplyr::mutate(high_var = ifelse(sd>1,1,0))
df_raw_use = distinct(df_raw[c('param_set_id','sd','mean','min','max','abs_diff','high_var')])
df_raw_use = df_raw_use %>% inner_join(df_params, by='param_set_id')

# ------------
df_raw_plot_nz = df_raw_plot_nz %>% inner_join(var_df, by='param_set_id')

plot(df_raw_plot_nz$sd, df_raw_plot_nz$mean)
abline(lm(mean ~ sd, data = df_raw_plot_nz), col = "red")

plot(df_raw_use$sd, df_raw_use$min)
abline(lm(min ~ sd, data = df_raw_use), col = "red")

# Filter for the two groups
df_plot = df_raw_use %>%
  filter(high_var %in% c(0, 1)) %>%
  mutate(high_var = factor(high_var, labels = c("Low Variance", "High Variance")))

ggplot(df_plot, aes(x = mean, fill = high_var)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Low Variance" = "blue", "High Variance" = "red")) +
  theme_minimal()

ggplot(df_plot, aes(x = activity_engulf_M1_baseline, fill = high_var)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Low Variance" = "blue", "High Variance" = "red")) +
  theme_minimal()

ggplot(df_plot, aes(x = activity_engulf_M2_baseline, fill = high_var)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Low Variance" = "blue", "High Variance" = "red")) +
  theme_minimal()

ggplot(df_plot, aes(x = SAMPs_decay, fill = high_var)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Low Variance" = "blue", "High Variance" = "red")) +
  theme_minimal()

ggplot(df_plot, aes(x = (SAMPs_decay-DAMPs_decay), fill = high_var)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Low Variance" = "blue", "High Variance" = "red")) +
  theme_minimal()

# Low treg_discrimination_efficiency - Tregs can't distinguish commensals
# This is actually more sophisticed and gives a bimodal dist. probably because 
# too low leads to bad outcomes consistently, 
# and too high leads to good outcomes consistently?

ggplot(df_plot, aes(x = treg_discrimination_efficiency, fill = high_var)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Low Variance" = "blue", "High Variance" = "red")) +
  theme_minimal()

df_plot_low = df_plot %>% dplyr::filter(high_var=='High Variance')
hist(df_plot_low$treg_discrimination_efficiency)
plot(df_plot_low$treg_discrimination_efficiency, df_plot_low$mean)

ggplot(df_plot, aes(x = (diffusion_speed_DAMPs), fill = high_var)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Low Variance" = "blue", "High Variance" = "red")) +
  theme_minimal()

ggplot(df_plot, aes(x = (activation_threshold_SAMPs), fill = high_var)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Low Variance" = "blue", "High Variance" = "red")) +
  theme_minimal()

ggplot(df_plot, aes(x = (SAMPs_decay), fill = high_var)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Low Variance" = "blue", "High Variance" = "red")) +
  theme_minimal()

# ------- VERY CLEAR?
ggplot(df_plot, aes(x = activation_threshold_DAMPs, fill = high_var)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Low Variance" = "blue", "High Variance" = "red")) +
  theme_minimal()

ggplot(df_plot, aes(x = activation_threshold_DAMPs-activation_threshold_SAMPs, fill = high_var)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Low Variance" = "blue", "High Variance" = "red")) +
  theme_minimal()

ggplot(df_plot, aes(x = SAMPs_decay, fill = high_var)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Low Variance" = "blue", "High Variance" = "red")) +
  theme_minimal()

# ------- VERY CLEAR?
ggplot(df_plot, aes(x = diffusion_speed_DAMPs, fill = high_var)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Low Variance" = "blue", "High Variance" = "red")) +
  theme_minimal()

ggplot(df_plot, aes(x = diffusion_speed_SAMPs, fill = high_var)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Low Variance" = "blue", "High Variance" = "red")) +
  theme_minimal()

ggplot(df_plot, aes(x = activation_threshold_DAMPs, fill = high_var)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Low Variance" = "blue", "High Variance" = "red")) +
  theme_minimal()


### epith_recovery_chance (line 769): If too low, injury never heals, commensal leakage persists
ggplot(df_plot, aes(x = epith_recovery_chance, fill = high_var)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Low Variance" = "blue", "High Variance" = "red")) +
  theme_minimal()

### recoevry in general?
ggplot(df_plot, aes(x = th_ROS_epith_recover+epith_recovery_chance, fill = high_var)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Low Variance" = "blue", "High Variance" = "red")) +
  theme_minimal()

###### MORE COMPOSITE METRICS
df_plot = df_plot %>% dplyr::mutate(DAMP_potential = (add_DAMPs / DAMPs_decay) / activation_threshold_DAMPs)
df_plot = df_plot %>% dplyr::mutate(SAMP_potential = (add_SAMPs / SAMPs_decay) / activation_threshold_SAMPs)
df_plot = df_plot %>% dplyr::mutate(Signal_Ratio = DAMP_potential / SAMP_potential)
df_plot = df_plot %>% dplyr::mutate(ROS_potential = (activity_ROS_M1_baseline / ros_decay) / th_ROS_epith_recover)
df_plot = df_plot %>% dplyr::mutate(Feedback_Strength = rate_leak_commensal_injury / epith_recovery_chance)
df_plot = df_plot %>% dplyr::mutate(Chronicity_Risk = Signal_Ratio * log(1 + Feedback_Strength))
df_plot = df_plot %>% dplyr::mutate(Inflammatory_Index = (add_DAMPs / DAMPs_decay / activation_threshold_DAMPs) /
                                      (add_SAMPs / SAMPs_decay / activation_threshold_SAMPs))
df_plot = df_plot %>% dplyr::mutate(Inflammatory_Index_Extended = Inflammatory_Index *
                                      sqrt(diffusion_speed_SAMPs / diffusion_speed_DAMPs))
df_plot = df_plot %>% dplyr::mutate(Instability_Score = Signal_Ratio *
                                      ROS_potential *
                                      (1 / epith_recovery_chance) *
                                      rate_leak_commensal_injury *
                                      exp(-2 * treg_discrimination_efficiency))



ggplot(df_plot, aes(x = log(Instability_Score), fill = high_var)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Low Variance" = "blue", "High Variance" = "red")) +
  theme_minimal()

#----------------------------------
pca_data = df_plot %>%
  group_by(param_set_id) %>%
  summarise(
    standard_dev = sd,
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
pca_matrix = pca_data %>%
  select(add_DAMPs_norm, add_SAMPs_norm, inv_threshold_DAMPs, inv_threshold_SAMPs,
         diffusion_speed_DAMPs, diffusion_speed_SAMPs,
         rate_leak_commensal_injury, inv_epith_recovery,
         ROS_production, inv_ROS_threshold, DAMPs_decay, SAMPs_decay,
         treg_discrimination) %>%
  as.matrix()

pca_result = prcomp(pca_matrix, scale. = TRUE, center = TRUE)

# Add PC scores to data
pca_data$PC1 = pca_result$x[, 1]
pca_data$PC2 = pca_result$x[, 2]
pca_data$PC3 = pca_result$x[, 3]

# Plot PCA biplot
library(viridis)

p_pca = ggplot(pca_data, aes(x = PC1, y = PC2, color = standard_dev)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_viridis(name = "Chronic\nFraction", limits = c(0, 1)) +
  theme_minimal() +
  labs(title = "PCA: Parameter Space Colored by Chronic Outcome",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "% variance)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "% variance)"))
p_pca
