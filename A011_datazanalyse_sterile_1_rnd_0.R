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
ss_start_threshold   = 450
t_max                = 500
tol_in               = 25*0.25

df_raw    = readRDS('/Users/burcutepekule/Desktop/tregs/all_comparison_results_sterile_1_trnd_0_old.rds')
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

table(df_raw_plot_nz$param_set_id)

x   = df_raw$mean_diff
round(100*sum(x>tol_in)/length(x),2)
round(100*sum(x<=tol_in & x>=(-1*tol_in))/length(x),2)
round(100*sum(x<(-1*tol_in))/length(x),2)
hist(x,30)

#-----------
var_df = df_raw %>%
  group_by(param_set_id) %>%
  summarise(
    variance = var(mean_diff),
    sd = sd(mean_diff),
    mean = mean(mean_diff),
    min = min(mean_diff),
    max = max(mean_diff),
    n_replicates = n()
  ) %>%
  arrange(desc(variance))


df_raw = df_raw %>% inner_join(var_df, by='param_set_id')
df_raw = df_raw %>% dplyr::mutate(high_var = ifelse(sd>1,1,0))
df_raw_use = distinct(df_raw[c('param_set_id','sd','mean','min','max','high_var')])
df_raw_use = df_raw_use %>% inner_join(df_params, by='param_set_id')

# Filter for the two groups
df_plot = df_raw_use %>%
  filter(high_var %in% c(0, 1)) %>%
  mutate(high_var = factor(high_var, labels = c("Low Variance", "High Variance")))

ggplot(df_plot, aes(x = SAMPs_decay, fill = high_var)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Low Variance" = "blue", "High Variance" = "red")) +
  theme_minimal()

ggplot(df_plot, aes(x = log10(SAMPs_decay/DAMPs_decay), fill = high_var)) +
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
df_plot_low = df_plot %>% dplyr::filter(high_var=='Low Variance')
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

ggplot(df_plot, aes(x = (activation_threshold_DAMPs), fill = high_var)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Low Variance" = "blue", "High Variance" = "red")) +
  theme_minimal()

ggplot(df_plot, aes(x = (activation_threshold_DAMPs-activation_threshold_SAMPs), fill = high_var)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Low Variance" = "blue", "High Variance" = "red")) +
  theme_minimal()

ggplot(df_plot, aes(x = (DAMPs_decay-SAMPs_decay), fill = high_var)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Low Variance" = "blue", "High Variance" = "red")) +
  theme_minimal()

ggplot(df_plot, aes(x = (diffusion_speed_DAMPs-diffusion_speed_SAMPs), fill = high_var)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Low Variance" = "blue", "High Variance" = "red")) +
  theme_minimal()

ggplot(df_plot, aes(x = (diffusion_speed_SAMPs), fill = high_var)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Low Variance" = "blue", "High Variance" = "red")) +
  theme_minimal()

ggplot(df_plot, aes(x = activation_threshold_DAMPs-(add_DAMPs/DAMPs_decay), fill = high_var)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Low Variance" = "blue", "High Variance" = "red")) +
  theme_minimal()

ggplot(df_plot, aes(x = activation_threshold_DAMPs-(add_DAMPs-diffusion_speed_DAMPs-DAMPs_decay), fill = high_var)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Low Variance" = "blue", "High Variance" = "red")) +
  theme_minimal()

### epith_recovery_chance (line 769): If too low, injury never heals, commensal leakage persists
ggplot(df_plot, aes(x = (epith_recovery_chance), fill = high_var)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Low Variance" = "blue", "High Variance" = "red")) +
  theme_minimal()

