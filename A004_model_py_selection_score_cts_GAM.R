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

x        = df_model$mean_diff
tol_in   = 25*0.25
round(100*sum(x<=tol_in & x>=(-1*tol_in))/length(x),2)
round(100*sum(x>tol_in)/length(x),2)
round(100*sum(x<(-1*tol_in))/length(x),2)
hist(x)


# Binary: did ANY "better" outcome occur?
df_summary$any_better <- df_summary$n_better > 0
df_summary$any_worse <- df_summary$n_worse > 0

# Prepare predictor matrix
X <- model.matrix(~ th_ROS_microbe + th_ROS_epith_recover + 
                    epith_recovery_chance + rat_com_pat_threshold + 
                    diffusion_speed_DAMPs + diffusion_speed_SAMPs + 
                    diffusion_speed_ROS + add_ROS + add_DAMPs + add_SAMPs + 
                    ros_decay + DAMPs_decay + SAMPs_decay + 
                    activation_threshold_DAMPs + activation_threshold_SAMPs + 
                    activity_engulf_M0_baseline + activity_engulf_M1_baseline + 
                    activity_engulf_M2_baseline + activity_ROS_M1_baseline + 
                    rate_leak_commensal_injury + rate_leak_pathogen_injury + 
                    rate_leak_commensal_baseline + active_age_limit + 
                    treg_discrimination_efficiency - 1, 
                  data = df_summary)
library(glmnet)

# Regularized logistic regression (prevents overfitting with rare events)
model_better_lasso <- cv.glmnet(X, df_summary$any_better, 
                                family = "binomial",
                                alpha = 1,  # Lasso
                                type.measure = "auc")

model_worse_lasso <- cv.glmnet(X, df_summary$any_worse,
                               family = "binomial", 
                               alpha = 1,
                               type.measure = "auc")

# Check performance
print(paste("Better AUC:", max(model_better_lasso$cvm)))
print(paste("Worse AUC:", max(model_worse_lasso$cvm)))

# Get important predictors
coef_better <- coef(model_better_lasso, s = "lambda.min")
coef_worse <- coef(model_worse_lasso, s = "lambda.min")

# Non-zero coefficients
important_better <- coef_better[coef_better[,1] != 0, , drop = FALSE]
important_worse <- coef_worse[coef_worse[,1] != 0, , drop = FALSE]

print("Parameters associated with 'better' outcomes:")
print(important_better)

print("Parameters associated with 'worse' outcomes:")
print(important_worse)

#-------------------------------------------------------------------------------

ggplot(df_summary, aes(x = th_ROS_epith_recover)) +
  geom_density(aes(fill = "drift", color = "drift"), alpha = 0.3) +
  geom_density(data = df_summary %>% filter(any_better), 
               aes(fill = "better", color = "better"), alpha = 0.5) +
  geom_density(data = df_summary %>% filter(any_worse),
               aes(fill = "worse", color = "worse"), alpha = 0.5) +
  scale_fill_manual(values = c("drift" = "gray", "better" = "green", "worse" = "red")) +
  scale_color_manual(values = c("drift" = "gray", "better" = "green", "worse" = "red")) +
  labs(title = "ROS threshold for epithelial recovery",
       subtitle = "Lower values → escape from drift (both directions)",
       x = "th_ROS_epith_recover", y = "Density") +
  theme_minimal()


#-------------------------------------------------------------------------------

# Visualize this "fork in the road"
df_summary %>%
  mutate(outcome = case_when(
    any_better ~ "Better",
    any_worse ~ "Worse",
    TRUE ~ "Drift"
  )) %>%
  ggplot(aes(x = th_ROS_epith_recover, fill = outcome)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Drift" = "gray", "Better" = "green", "Worse" = "red")) +
  geom_vline(xintercept = median(df_summary$th_ROS_epith_recover[df_summary$any_better | df_summary$any_worse]),
             linetype = "dashed") +
  labs(title = "ROS Threshold: The Drift Escape Parameter",
       subtitle = "Low values → escape drift, but direction depends on other parameters",
       x = "th_ROS_epith_recover") +
  theme_minimal()

#-------------------------------------------------------------------------------

# th_ROS_epith_recover: OR = 0.54
# 
# Going from min→max threshold reduces odds of "better" by 46%
# Low ROS recovery threshold nearly doubles the odds of better outcomes
# 
# th_ROS_epith_recover: OR = 0.51
# 
# Same direction as "better"!
# Low ROS threshold doubles the odds of BOTH better and worse
# This is your "instability" parameter

#-------------------------------------------------------------------------------

df_summary_th_ros_low  = df_summary %>% dplyr::filter(th_ROS_epith_recover<0.2)
df_summary_th_ros_high = df_summary %>% dplyr::filter(th_ROS_epith_recover>=0.2)
