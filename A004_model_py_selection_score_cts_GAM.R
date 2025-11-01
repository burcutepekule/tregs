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
df_summary = readRDS('/Users/burcutepekule/Desktop/tregs/df_summary_selection_score_cts.rds')
df_raw     = readRDS('/Users/burcutepekule/Desktop/tregs/df_summary_selection_score_raw.rds')

df_model = df_raw %>% dplyr::filter(comparison=='Treg_OFF_ON' & injury_type=='sterile')
df_model = merge(df_model, df_params, by='param_set_id')

df_model = df_model %>% dplyr::filter(effect_size %in% c('Medium','Large') & abs(mean_diff)>4*tol)

df_model_significant = df_model %>% dplyr::mutate(nonzero=ifelse(effect_size %in% c('Medium','Large') & abs(mean_diff)>4*tol, 1, 0))


df_model = df_model_significant %>% dplyr::filter(comparison=='Treg_OFF_ON' & injury_type=='sterile')


# df_summary_merged = df_summary_merged %>% dplyr::filter(outcome_sd<1)

# df_model = df_summary_merged %>% dplyr::filter(comparison=='Treg_OFF_ON' & injury_type=='pathogenic')

hist(df_model$outcome_mean)

gam_model = gam(outcome_mean ~ 
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
              s(treg_discrimination_efficiency, bs = "cs"),
              # + comparison + injury_type, # if also modelled as factors
            data = df_model,
            method = "REML")

summary(gam_model)      # see which predictors matter
plot(gam_model, pages=1, shade=TRUE)   # visualize smooth effects
gam.check(gam_model)    # residual diagnostics

xlibrary(pROC)
roc(df_model$outcome_mean, fitted(gam_model))$auc

dev.off()
plot(gam_model, select = 2, shade = TRUE, seWithMean = TRUE, rug = TRUE)
