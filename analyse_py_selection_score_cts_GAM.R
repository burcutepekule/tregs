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
df_params  = read_csv('/Users/burcutepekule/Desktop/tregs/mass_sim_results/sampled_parameters.csv', show_col_types = FALSE)
df_summary = readRDS('/Users/burcutepekule/Desktop/tregs/df_summary_selection_score_cts_full_range.rds')
df_summary_merged = merge(df_summary, df_params, by='param_set_id')

df_model = df_summary_merged %>% dplyr::filter(comparison=='Treg_OFF_ON' & injury_type=='sterile')

df_model$nonzero = as.numeric(df_model$outcome_mean != 0)

gam_binary = gam(nonzero ~ 
              s(th_ROS_microbe, bs = "cs") +
              s(th_ROS_epith_recover, bs = "cs") +
              s(epith_recovery_chance, bs = "cs") +
              s(rat_com_pat_threshold, bs = "cs") +
              s(diffusion_speed_DAMPs, bs = "cs") +
              s(diffusion_speed_SAMPs, bs = "cs") +
              s(diffusion_speed_ROS, bs = "cs") +
              s(add_ROS, bs = "cs") +
              s(add_DAMPs, bs = "cs") +
              s(add_SAMPs, bs = "cs") +
              s(ros_decay, bs = "cs") +
              s(DAMPs_decay, bs = "cs") +
              s(SAMPs_decay, bs = "cs") +
              s(activation_threshold_DAMPs, bs = "cs") +
              s(activation_threshold_SAMPs, bs = "cs") +
              s(activity_engulf_M0_baseline, bs = "cs") +
              s(activity_engulf_M1_baseline, bs = "cs") +
              s(activity_engulf_M2_baseline, bs = "cs") +
              s(activity_ROS_M1_baseline, bs = "cs") +
              s(rate_leak_commensal_injury, bs = "cs") +
              s(rate_leak_pathogen_injury, bs = "cs") +
              s(rate_leak_commensal_baseline, bs = "cs") +
              s(active_age_limit, bs = "cs") +
              s(treg_discrimination_efficiency, bs = "cs"), # + comparison + injury_type, # if also modelled as factors
            data = df_model,
            family = binomial)

summary(gam_binary)      # see which predictors matter
plot(gam_binary, pages=1, shade=TRUE)   # visualize smooth effects
gam.check(gam_binary)    # residual diagnostics

library(pROC)
roc(df_model$nonzero, fitted(gam_binary))$auc

plot(gam_binary, select = 3, pages = 24, shade = TRUE)

