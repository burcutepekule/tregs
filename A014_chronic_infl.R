rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)  # For read_csv
library(stringr)
library(zoo)
library(mgcv)
library(factoextra)
library(cluster)

source("/Users/burcutepekule/Dropbox/Treg_problem_v2/MISC/PLOT_FUNCTIONS.R")
# source("~/Desktop/tregs/A013_datazanalyse_sterile_1_rnd_0.R")
df_plot      = readRDS('df_plot_A13.rds') # comes from 
df_plot_keep = df_plot
# vector of parameter names
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

df_plot = df_plot[c('param_set_id','replicate_id',
                    'ss_start_tregs_off','ss_start_tregs_on',
                    'mean_treg_on','mean_treg_off')]


df_summary = df_plot %>%
  group_by(param_set_id) %>%
  summarise(
    mean_treg_on_mean  = mean(mean_treg_on,  na.rm = TRUE),
    mean_treg_off_mean = mean(mean_treg_off, na.rm = TRUE),
    mean_treg_on_off_diff_mean = mean_treg_on_mean-mean_treg_off_mean
  )


tol_in = 25*0.25
df_summary = df_summary %>% dplyr::mutate(rounded_on_off_diff_mean = ifelse(
  mean_treg_on_off_diff_mean < tol_in & mean_treg_on_off_diff_mean > -1*tol_in, 
  0, mean_treg_on_off_diff_mean
))


df_plot = merge(df_plot, df_summary, by='param_set_id')

