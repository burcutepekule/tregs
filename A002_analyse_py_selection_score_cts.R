rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)  # For read_csv
library(stringr)
library(zoo)

source("/Users/burcutepekule/Dropbox/Treg_problem_v2/MISC/PLOT_FUNCTIONS.R")
# df_raw      = readRDS('/Users/burcutepekule/Desktop/tregs/all_comparison_results_0_full_range_RND.rds')
# df_raw      = readRDS('/Users/burcutepekule/Desktop/tregs/all_comparison_results_0_full_range_LHS.rds')
t_max  = 500
df_raw = readRDS('/Users/burcutepekule/Desktop/tregs/all_comparison_results_0.rds')
df_raw = df_raw %>% dplyr::mutate(log_integ_diff = sign(integ_diff)*log(1+abs(integ_diff))/(t_max-ss_start))
hist(df_raw$log_integ_diff, 15)
length(unique(df_raw$param_set_id))


df_raw_keep = df_raw

#----- filter based on ss_start, it cannot be too large otherwise not much to compare!
ss_start_threshold   = 150
param_id_all_below = df_raw %>%
  dplyr::group_by(param_set_id) %>%
  dplyr::summarise(all_below = all(ss_start < ss_start_threshold), .groups = "drop") %>%
  dplyr::filter(all_below) %>%
  dplyr::pull(param_set_id)
df_raw = df_raw %>% dplyr::filter(param_set_id %in% param_id_all_below)

#----- filter based on replicate_id, less than 10 means incomplete!
param_id_all_complete = df_raw %>%
  dplyr::group_by(param_set_id, comparison, injury_type) %>%
  dplyr::summarise(all_complete = (n_distinct(replicate_id) == 10), .groups = "drop") %>%
  dplyr::group_by(param_set_id) %>%
  dplyr::summarise(all_complete = all(all_complete), .groups = "drop") %>%
  dplyr::filter(all_complete) %>%
  dplyr::pull(param_set_id)
df_raw = df_raw %>% dplyr::filter(param_set_id %in% param_id_all_complete)

length(unique(df_raw$param_set_id)) 
hist(df_raw$log_integ_diff, 15)

df_raw     = df_raw %>% dplyr::mutate(abs_cohens_d = abs(cohens_d))
df_raw     = df_raw %>% dplyr::mutate(effect_size = case_when(
  abs_cohens_d < 0.2 ~ "Negligible",
  abs_cohens_d < 0.5 & abs_cohens_d>= 0.2  ~ "Small",
  abs_cohens_d < 0.8 & abs_cohens_d>= 0.5 ~ "Medium",
  TRUE ~ "Large"
))

# meaningful_effect_size = 25*0.05 # at least 5%
#
# df_raw = df_raw %>% dplyr::mutate(better = ifelse(((effect_size=='Large' | effect_size=='Medium') & mean_diff>(+1)*meaningful_effect_size),1,0))
# df_raw = df_raw %>% dplyr::mutate(worse  = ifelse(((effect_size=='Large' | effect_size=='Medium') & mean_diff<(-1)*meaningful_effect_size),1,0))
# 
# df_raw = df_raw %>% dplyr::mutate(outcome = ifelse(better==0 & worse==0, 'drift',
#                                                    ifelse(better == 0 & worse == 1, 'worse',
#                                                           ifelse( better == 1 & worse == 0, 'better', NA))))

df_raw = df_raw %>% dplyr::mutate(outcome = sign(mean_diff)*log(1+abs_cohens_d*abs(mean_diff)))
# df_raw = df_raw %>% dplyr::mutate(outcome = sign(mean_diff)*log10(1+abs_cohens_d*abs(mean_diff)))

hist(df_raw$outcome,30)
hist(df_raw$mean_diff,30)
# hist(sign(df_raw$mean_diff)*log10(abs(1+df_raw$mean_diff)),30)
table(round(df_raw$mean_diff,3))

tol = 0 # at least 5%
hist(x)
x   = df_raw$log_integ_diff
round(100*sum(x<=tol & x>=(-1*tol))/length(x),2)
round(100*sum(x>tol)/length(x),2)
round(100*sum(x<(-1*tol))/length(x),2)

# df_raw_non_negligible = df_raw %>% dplyr::filter(effect_size %in% c("Medium", "Large"))


df_summary = df_raw %>%
  dplyr::group_by(param_set_id, injury_type, comparison) %>%
  dplyr::summarise(
    outcome_mean = mean(outcome, na.rm = TRUE),
    n_better = sum(effect_size %in% c("Medium", "Large") & mean_diff > tol, na.rm = TRUE),
    n_drift  = sum(effect_size %in% c("Medium", "Large") & (mean_diff <= tol & mean_diff >= -1*tol), na.rm = TRUE),
    n_worse  = sum(effect_size %in% c("Medium", "Large") & mean_diff < -1*tol, na.rm = TRUE),
    n_negligible = sum(effect_size %in% c("Negligible", "Small")),
    .groups = "drop"
  )


df_summary = df_summary %>% dplyr::mutate(selection_score = outcome_mean)
df_summary$tol = tol
df_raw$tol = tol

# saveRDS(df_summary, '/Users/burcutepekule/Desktop/tregs/df_summary_selection_score_cts_full_range.rds')
saveRDS(df_summary, paste0('/Users/burcutepekule/Desktop/tregs/df_summary_selection_score_cts_',ss_start_threshold,'.rds'))
saveRDS(df_raw, paste0('/Users/burcutepekule/Desktop/tregs/df_summary_selection_score_raw_',ss_start_threshold,'.rds'))

