rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)  # For read_csv
library(stringr)
library(zoo)

source("/Users/burcutepekule/Dropbox/Treg_problem_v2/MISC/PLOT_FUNCTIONS.R")
inj_type             = 'sterile'
ss_start_threshold   = 450
t_max                = 500
tol_in               = 25*0.25

df_raw = readRDS('/Users/burcutepekule/Desktop/tregs/all_comparison_results_0_pre.rds')
df_params_pre = read_csv('/Users/burcutepekule/Desktop/tregs/mass_sim_results_presampled/loaded_parameters.csv', show_col_types = FALSE)

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

# df_raw_plot_nz = df_raw %>% dplyr::filter(effect_size %in% c("Medium","Large") & (mean_diff>tol_in | mean_diff<(-1*tol_in)))
# df_raw_plot_z  = df_raw %>% dplyr::filter((mean_diff<=tol_in | mean_diff>=(-1*tol_in))) # this is not true- what about effect_size small but abs(mean_diff)>tol_in?
# df_raw_plot    = rbind(df_raw_plot_z, df_raw_plot_nz)
# hist(df_raw_plot$mean_diff,30)

x   = df_raw$mean_diff
round(100*sum(x>tol_in)/length(x),2)
round(100*sum(x<=tol_in & x>=(-1*tol_in))/length(x),2)
round(100*sum(x<(-1*tol_in))/length(x),2)
hist(x,30)

