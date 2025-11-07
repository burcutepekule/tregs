rm(list=ls())
library(dplyr)
library(tidyr)
library(readr)
source("/Users/burcutepekule/Dropbox/Treg_problem_v2/MISC/PLOT_FUNCTIONS.R")

ss_start_threshold = 250
df_params = read_csv('/Users/burcutepekule/Desktop/tregs/mass_sim_results/sampled_parameters.csv', show_col_types = FALSE)
df_raw    = readRDS(paste0('/Users/burcutepekule/Desktop/tregs/df_summary_selection_score_raw_',ss_start_threshold,'.rds'))
df_model  = df_raw %>% dplyr::filter(comparison=='Treg_OFF_ON' & injury_type=='sterile')
# df_model  = df_model %>% dplyr::mutate(nonzero=ifelse(abs(mean_diff)>tol, 1, 0))
df_model  = df_model %>% dplyr::mutate(nonzero=ifelse(effect_size %in% c('Medium','Large') & abs(mean_diff)>tol, 1, 0))
df_model  = df_model %>% dplyr::mutate(mean_diff_effect=mean_diff*nonzero)
df_model  = merge(df_model, df_params, by='param_set_id')
df_model  = df_model %>% dplyr::mutate(group_th = ifelse(th_ROS_epith_recover<0.05 & epith_recovery_chance<0.05,1,0))

df_model_1 = df_model %>% dplyr::filter(group_th==1)
df_model_0 = df_model %>% dplyr::filter(group_th==0)

ggplot(df_model, aes(x = mean_diff_effect, fill = factor(group_th))) +
  geom_density(alpha = 0.4, color = NA) +
  scale_fill_manual(values = c("0" = "gray70", "1" = "#3182bd"),
                    labels = c("Zero outcome", "Nonzero outcome")) +
  labs(
    x = 'Epith. score difference',
    y = "Density",
    fill = "Outcome",
    title = ""
  ) +
  theme_minimal(base_size = 13)

# table(round(df_model_1$mean_diff_effect))
100*round(length(which(df_model_1$mean_diff_effect==0))/dim(df_model_1)[1],3)
# table(round(df_model_0$mean_diff_effect))
100*round(length(which(df_model_0$mean_diff_effect==0))/dim(df_model_0)[1],3)